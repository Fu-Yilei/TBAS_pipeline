from __future__ import annotations

import csv
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

BED_FILE_MAP = {
    "Arthrogryposis": "adaptive_regions/annotated/merged/Arthro_2582_padded_top11_manual_annotated.sorted.merged.bed",
    "Epilepsy": "adaptive_regions/annotated/merged/Epilepsy_V2_1662_padded_top8_manual_annotated.sorted.merged.bed",
    "Microcephaly": "adaptive_regions/annotated/merged/Microcephaly_1769_padded_top3_manual_annotated.sorted.merged.bed",
    "Neurodegeneration": "adaptive_regions/annotated/merged/Neurodegen_V2_1704_padded_top8_manual_annotated.sorted.merged.bed",
    "CMRG": "adaptive_regions/annotated/merged/GRCh38_CMRG_benchmark_gene_coordinates_padding10K.bed",
    "WGS": None,
}

TR_BED_FILE_MAP = {
    "Arthrogryposis": "analysis/tr_regions/adotto_catalog.hg38.lite.Arthro.bed",
    "Epilepsy": "analysis/tr_regions/adotto_catalog.hg38.lite.Epilepsy.bed",
    "Microcephaly": "analysis/tr_regions/adotto_catalog.hg38.lite.Microcephaly.bed",
    "Neurodegeneration": "analysis/tr_regions/adotto_catalog.hg38.lite.Neurodegen.bed",
    "CMRG": "analysis/tr_regions/adotto_catalog.hg38.lite.CMRG.bed",
    "WGS": "adotto_catalog.hg38.lite.bed",
}

STAGE_ORDER = [
    "demultiplex",
    "fastq_extract",
    "minimap2",
    "bam_sort",
    "sniffles_global",
    "bam_mosdepth",
    "clair3_local",
    "kanpig_pileup",
    "kanpig_gt",
    "kanpig_trio",
    "whatshap_single_sample_local_phasing",
    "whatshap_haplotag",
    "medaka_local",
    "medaka_patho",
    "modkit",
    "modkit_nohp",
    "tdb",
]


@dataclass(frozen=True)
class PipelineSettings:
    manifest: Path
    output_folder: Path = Path("GREGoR_adaptive_sampling")
    reference: Path = Path("GRCh38-2.1.0_genome_mainchrs.fa")
    kit_name: str = "SQK-NBD114-96"
    clair3_model: Path = Path("clair3/bin/models/r1041_e82_400bps_sup_v520")
    adotto_pheno: Path = Path("adaptive_regions/patho_tr/adotto_pheno.bed")


class CommandRunner:
    def __init__(self, dry_run: bool = False, print_commands: bool = True) -> None:
        self.dry_run = dry_run
        self.print_commands = print_commands

    def run(self, command: str) -> None:
        if self.print_commands:
            print(command)
        if self.dry_run:
            return
        subprocess.run(command, shell=True, executable="/bin/bash", check=True)


def _quote(value: object) -> str:
    return shlex.quote(str(value))


def _detect_delimiter(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".tsv":
        return "\t"
    if suffix == ".csv":
        return ","
    if suffix in {".txt", ".manifest"}:
        first_line = path.read_text().splitlines()[0] if path.exists() else ""
        return "\t" if "\t" in first_line else ","
    raise ValueError(
        f"Unsupported manifest extension: {path.suffix}. Use .csv or .tsv."
    )


def load_manifest(path: Path) -> list[dict[str, str]]:
    delimiter = _detect_delimiter(path)
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")
        for line_no, row in enumerate(reader, start=2):
            normalized = {
                (key or "").strip(): (value or "").strip()
                for key, value in row.items()
            }
            sample_id = normalized.get("sample_id", "")
            if not sample_id:
                raise ValueError(f"Missing sample_id at manifest line {line_no}")
            rows.append(normalized)
    if not rows:
        raise ValueError(f"Manifest is empty: {path}")
    return rows


def parse_sample_barcodes(sample_id: str) -> list[str]:
    sample_text = str(sample_id)
    match = re.search(r"(\d+)-(\d+)", sample_text)
    if match:
        start, end = int(match.group(1)), int(match.group(2))
    else:
        match = re.search(r"(\d+)_(\d+)", sample_text)
        if not match:
            raise ValueError(
                f"Cannot parse barcode range from sample_id={sample_id!r}. "
                "Expected format like '10_12_*' or '43-45'."
            )
        start, end = int(match.group(1)), int(match.group(2))

    if start > end:
        raise ValueError(
            f"Invalid barcode range in sample_id={sample_id!r}: start > end"
        )
    return [f"barcode{i:02d}" for i in range(start, end + 1)] + ["unassigned"]


def parse_readgroup_prefix_from_header(header_text: str, kit_name: str) -> str:
    pattern = re.compile(rf"(.+)_{re.escape(kit_name)}_barcode\d{{2}}$")
    for line in header_text.splitlines():
        if not line.startswith("@RG"):
            continue
        for field in line.split("\t"):
            if not field.startswith("ID:"):
                continue
            rg_id = field[3:]
            match = pattern.search(rg_id)
            if match:
                return match.group(1)
    raise ValueError("Could not infer read-group prefix from BAM header")


def get_readgroup_prefix_from_bam(bam_path: Path, kit_name: str) -> str:
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found for RG parsing: {bam_path}")
    proc = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        text=True,
        capture_output=True,
        check=True,
    )
    return parse_readgroup_prefix_from_header(proc.stdout, kit_name)


def _resolve_bed_path(value: str) -> str | None:
    if not value:
        return None
    if value.lower() in {"none", "nan", "na"}:
        return None
    return BED_FILE_MAP.get(value, value)


def _resolve_tr_bed_path(value: str) -> str | None:
    if not value:
        return None
    if value.lower() in {"none", "nan", "na"}:
        return None
    return TR_BED_FILE_MAP.get(value, value)


def _find_calls_bam(sample_dir: Path) -> Path:
    matches = sorted(sample_dir.glob("calls*.bam"))
    if not matches:
        raise FileNotFoundError(f"No calls*.bam found in {sample_dir}")
    return matches[0]


class TBASPipeline:
    def __init__(
        self,
        manifest_rows: Sequence[dict[str, str]],
        settings: PipelineSettings,
        runner: CommandRunner | None = None,
    ) -> None:
        self.rows = list(manifest_rows)
        self.settings = settings
        self.runner = runner or CommandRunner()
        self._stage_handlers = {
            "demultiplex": self._stage_demultiplex,
            "fastq_extract": self._stage_fastq_extract,
            "minimap2": self._stage_minimap2,
            "bam_sort": self._stage_bam_sort,
            "sniffles_global": self._stage_sniffles_global,
            "bam_mosdepth": self._stage_bam_mosdepth,
            "clair3_local": self._stage_clair3_local,
            "kanpig_pileup": self._stage_kanpig_pileup,
            "kanpig_gt": self._stage_kanpig_gt,
            "kanpig_trio": self._stage_kanpig_trio,
            "whatshap_single_sample_local_phasing": self._stage_whatshap_phase,
            "whatshap_haplotag": self._stage_whatshap_haplotag,
            "medaka_local": self._stage_medaka_local,
            "medaka_patho": self._stage_medaka_patho,
            "modkit": self._stage_modkit,
            "modkit_nohp": self._stage_modkit_nohp,
            "tdb": self._stage_tdb,
        }

    @classmethod
    def from_manifest(
        cls,
        settings: PipelineSettings,
        runner: CommandRunner | None = None,
    ) -> "TBASPipeline":
        return cls(load_manifest(settings.manifest), settings, runner=runner)

    def run(self, stages: Sequence[str] | None = None) -> None:
        selected = list(stages or STAGE_ORDER)
        unknown = [stage for stage in selected if stage not in self._stage_handlers]
        if unknown:
            raise ValueError(f"Unknown stages: {', '.join(unknown)}")
        for stage in selected:
            self._stage_handlers[stage]()

    def _require(self, row: dict[str, str], key: str, stage: str) -> str:
        value = row.get(key, "").strip()
        if not value:
            raise ValueError(
                f"Manifest column '{key}' is required for stage '{stage}' "
                f"(sample_id={row.get('sample_id', '<missing>')})"
            )
        return value

    def _sample_dir(self, sample_id: str) -> Path:
        return self.settings.output_folder / sample_id

    def _trio_barcodes(self, sample_id: str) -> list[str]:
        barcodes = parse_sample_barcodes(sample_id)[:3]
        if len(barcodes) < 3:
            raise ValueError(
                f"Expected trio-style sample_id with at least 3 barcodes: {sample_id}"
            )
        return barcodes

    def _sample_tag(self, sample_id: str, barcode: str) -> str:
        return f"{sample_id}_{barcode}"

    def _calls_bam(self, row: dict[str, str], sample_dir: Path) -> Path:
        calls_bam = row.get("calls_bam", "").strip()
        if calls_bam:
            return Path(calls_bam)
        return _find_calls_bam(sample_dir)

    def _demux_prefix(self, row: dict[str, str], calls_bam: Path) -> str:
        prefix = row.get("read_group_prefix", "").strip()
        if prefix:
            return prefix
        return get_readgroup_prefix_from_bam(calls_bam, self.settings.kit_name)

    def _resolve_medaka_tr_bed(self, row: dict[str, str], stage: str) -> str:
        tr_bed_raw = row.get("tr_bed_file", "").strip()
        if tr_bed_raw:
            tr_bed = _resolve_tr_bed_path(tr_bed_raw)
            if tr_bed is None:
                raise ValueError(
                    f"Invalid tr_bed_file for sample_id={row.get('sample_id')}: "
                    f"{tr_bed_raw!r}"
                )
            return tr_bed

        bed_file = row.get("bed_file", "").strip()
        if not bed_file:
            raise ValueError(
                f"Manifest column 'tr_bed_file' or 'bed_file' is required for "
                f"stage '{stage}' (sample_id={row.get('sample_id', '<missing>')})"
            )

        tr_bed = _resolve_tr_bed_path(bed_file)
        if tr_bed is None:
            raise ValueError(
                f"No TR BED mapping for sample_id={row.get('sample_id')} "
                f"bed_file={bed_file}"
            )
        return tr_bed

    def _stage_demultiplex(self) -> None:
        stage = "demultiplex"
        threads = 5
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            sample_dir.mkdir(parents=True, exist_ok=True)
            calls_bam = self._calls_bam(row, sample_dir)
            rg_prefix = self._demux_prefix(row, calls_bam)
            for barcode in self._trio_barcodes(sample_id):
                full_rg = f"{rg_prefix}_{self.settings.kit_name}_{barcode}"
                output_bam = sample_dir / f"{sample_id}_{barcode}.bam"
                command = (
                    f"samtools view -b@ {threads} {_quote(calls_bam)} -r "
                    f"{_quote(full_rg)} -o {_quote(output_bam)}"
                )
                self.runner.run(command)

    def _stage_fastq_extract(self) -> None:
        stage = "fastq_extract"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_bam = sample_dir / f"{sample_id}_{barcode}.bam"
                output_fastq = sample_dir / f"{sample_id}_{barcode}.fastq"
                command = (
                    f"samtools fastq -@ {threads} -T RG,Mm,Ml,MM,ML -0 "
                    f"{_quote(output_fastq)} {_quote(input_bam)}"
                )
                self.runner.run(command)

    def _stage_minimap2(self) -> None:
        stage = "minimap2"
        threads = 20
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_fastq = sample_dir / f"{sample_id}_{barcode}.fastq"
                output_sam = sample_dir / f"{sample_id}_{barcode}.sam"
                command = (
                    f"minimap2 -t {threads} -a -Y -y -x map-ont -o "
                    f"{_quote(output_sam)} {_quote(self.settings.reference)} "
                    f"{_quote(input_fastq)}"
                )
                self.runner.run(command)

    def _stage_bam_sort(self) -> None:
        stage = "bam_sort"
        threads = 5
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_sam = sample_dir / f"{sample_id}_{barcode}.sam"
                output_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                command = (
                    f"samtools sort -@ {threads} -o {_quote(output_bam)} "
                    f"--write-index {_quote(input_sam)}"
                )
                self.runner.run(command)

    def _stage_sniffles_global(self) -> None:
        stage = "sniffles_global"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_vcf = sample_dir / f"{sample_id}_{barcode}_global.germline.vcf.gz"
                output_snf = sample_dir / f"{sample_id}_{barcode}_global.germline.snf"
                command = (
                    f"sniffles -i {_quote(input_bam)} --reference "
                    f"{_quote(self.settings.reference)} -v {_quote(output_vcf)} "
                    f"--allow-overwrite --snf {_quote(output_snf)} --sample-id "
                    f"{_quote(sample_tag)} -t {threads} --output-rnames"
                )
                self.runner.run(command)

    def _stage_bam_mosdepth(self) -> None:
        stage = "bam_mosdepth"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            bed_file_loc = _resolve_bed_path(self._require(row, "bed_file", stage))
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_dir = sample_dir / f"{sample_id}_{barcode}_mosdepth"
                output_dir.mkdir(parents=True, exist_ok=True)
                prefix = output_dir / barcode
                if bed_file_loc is None:
                    command = (
                        f"mosdepth -t {threads} -x {_quote(prefix)} "
                        f"{_quote(input_bam)}"
                    )
                else:
                    command = (
                        f"mosdepth -t {threads} -b {_quote(bed_file_loc)} -x "
                        f"{_quote(prefix)} {_quote(input_bam)}"
                    )
                self.runner.run(command)

    def _stage_clair3_local(self) -> None:
        stage = "clair3_local"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            bed_file_loc = _resolve_bed_path(self._require(row, "bed_file", stage))
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_dir = sample_dir / f"{sample_id}_{barcode}_clair3_local"
                output_dir.mkdir(parents=True, exist_ok=True)
                command = (
                    f"run_clair3.sh -b {_quote(input_bam)} -f "
                    f"{_quote(self.settings.reference)} -m "
                    f"{_quote(self.settings.clair3_model)} --sample_name="
                    f"{_quote(sample_tag)} -t {threads} -p ont -o "
                    f"{_quote(output_dir)}"
                )
                if bed_file_loc is not None:
                    command = f"{command} --bed_fn={_quote(bed_file_loc)}"
                self.runner.run(command)

    def _stage_kanpig_pileup(self) -> None:
        stage = "kanpig_pileup"
        threads = 4
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_plup = sample_dir / f"{sample_id}_{barcode}_kanpig.plup.gz"
                command = (
                    f"kanpig plup -b {_quote(input_bam)} -t {threads} | "
                    f"bedtools sort -header | bgzip > {_quote(output_plup)}"
                )
                self.runner.run(command)
                self.runner.run(
                    f"tabix -f -s 1 -b 2 -e 2 {_quote(output_plup)}"
                )

    def _stage_kanpig_gt(self) -> None:
        stage = "kanpig_gt"
        threads = 4
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            bed_file_loc = _resolve_bed_path(self._require(row, "bed_file", stage))
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_reads = sample_dir / f"{sample_id}_{barcode}_kanpig.plup.gz"
                input_vcf = sample_dir / f"{sample_id}_{barcode}_global.germline.vcf.gz"
                output_vcf = (
                    sample_dir
                    / f"{sample_id}_{barcode}_global.germline.kanpig.genotyped.vcf.gz"
                )
                command = (
                    f"kanpig gt --input {_quote(input_vcf)} -t {threads} --reads "
                    f"{_quote(input_reads)} --sample {_quote(sample_tag)} "
                    f"--reference {_quote(self.settings.reference)}"
                )
                if bed_file_loc is not None:
                    command = f"{command} --bed {_quote(bed_file_loc)}"
                command = f"{command} | bcftools sort -Oz -W -o {_quote(output_vcf)}"
                self.runner.run(command)

    def _stage_kanpig_trio(self) -> None:
        stage = "kanpig_trio"
        threads = 4
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            bed_file_loc = _resolve_bed_path(self._require(row, "bed_file", stage))
            sample_dir = self._sample_dir(sample_id)
            trio = self._trio_barcodes(sample_id)
            proband, mother, father = trio[0], trio[1], trio[2]
            proband_reads = sample_dir / f"{sample_id}_{proband}_kanpig.plup.gz"
            mother_reads = sample_dir / f"{sample_id}_{mother}_kanpig.plup.gz"
            father_reads = sample_dir / f"{sample_id}_{father}_kanpig.plup.gz"
            input_vcf = sample_dir / f"{sample_id}_{proband}_global.germline.vcf.gz"
            output_vcf = (
                sample_dir
                / f"{sample_id}_{proband}_global.germline.kanpig.trio.genotyped.vcf.gz"
            )
            command = (
                f"kanpig trio --input {_quote(input_vcf)} --proband "
                f"{_quote(proband_reads)} --mother {_quote(mother_reads)} "
                f"--father {_quote(father_reads)} -t {threads} --reference "
                f"{_quote(self.settings.reference)}"
            )
            if bed_file_loc is not None:
                command = f"{command} --bed {_quote(bed_file_loc)}"
            command = f"{command} | bcftools sort -Oz -W -o {_quote(output_vcf)}"
            self.runner.run(command)

    def _stage_whatshap_phase(self) -> None:
        stage = "whatshap_single_sample_local_phasing"
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_vcf = sample_dir / f"{sample_id}_{barcode}_clair3_local/merge_output.vcf.gz"
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_vcf = sample_dir / f"{sample_tag}.local.phased.vcf.gz"
                command = (
                    f"whatshap phase --reference {_quote(self.settings.reference)} "
                    f"-o {_quote(output_vcf)} --ignore-read-groups "
                    f"{_quote(input_vcf)} {_quote(input_bam)}"
                )
                self.runner.run(command)
                self.runner.run(f"tabix -f -p vcf {_quote(output_vcf)}")

    def _stage_whatshap_haplotag(self) -> None:
        stage = "whatshap_haplotag"
        threads = 3
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_vcf = sample_dir / f"{sample_tag}.local.phased.vcf.gz"
                input_bam = sample_dir / f"{sample_id}_{barcode}.sorted.bam"
                output_bam = sample_dir / f"{sample_id}_{barcode}.HP.bam"
                output_hp = sample_dir / f"{sample_id}_{barcode}.read.HP.gz"
                command = (
                    f"whatshap haplotag --ignore-read-groups --output-haplotag-list "
                    f"{_quote(output_hp)} --reference {_quote(self.settings.reference)} "
                    f"--output-threads {threads} -o {_quote(output_bam)} "
                    f"{_quote(input_vcf)} {_quote(input_bam)}"
                )
                self.runner.run(command)
                self.runner.run(f"samtools index -@ {threads} {_quote(output_bam)}")

    def _stage_medaka_local(self) -> None:
        stage = "medaka_local"
        threads = 10
        model = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0:consensus"
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            proband_gender = self._require(row, "proband_gender", stage).lower()
            tr_bed = self._resolve_medaka_tr_bed(row, stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_bam = sample_dir / f"{sample_id}_{barcode}.HP.bam"
                output_tr = sample_dir / f"{sample_id}_{barcode}_tr"
                output_tr.mkdir(parents=True, exist_ok=True)
                command = (
                    f"medaka tandem --workers {threads} --ignore_read_groups "
                    f"--model {_quote(model)} --phasing prephased --sample_name "
                    f"{_quote(sample_tag)} {_quote(input_bam)} "
                    f"{_quote(self.settings.reference)} {_quote(tr_bed)} "
                    f"{_quote(proband_gender)} {_quote(output_tr)}"
                )
                self.runner.run(command)

    def _stage_medaka_patho(self) -> None:
        stage = "medaka_patho"
        threads = 10
        model = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0:consensus"
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            proband_gender = self._require(row, "proband_gender", stage).lower()
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_bam = sample_dir / f"{sample_id}_{barcode}.HP.bam"
                output_tr = sample_dir / f"{sample_id}_{barcode}_pheno_tr"
                output_tr.mkdir(parents=True, exist_ok=True)
                command = (
                    f"medaka tandem --workers {threads} --ignore_read_groups "
                    f"--model {_quote(model)} --phasing prephased --sample_name "
                    f"{_quote(sample_tag)} {_quote(input_bam)} "
                    f"{_quote(self.settings.reference)} {_quote(self.settings.adotto_pheno)} "
                    f"{_quote(proband_gender)} {_quote(output_tr)}"
                )
                self.runner.run(command)

    def _stage_modkit(self) -> None:
        stage = "modkit"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_bam = sample_dir / f"{sample_id}_{barcode}.HP.bam"
                output_dir = sample_dir / f"{sample_id}_{barcode}_modkit"
                output_dir.mkdir(parents=True, exist_ok=True)
                prefix = f"{sample_id}_{barcode}.HP"
                output_file = output_dir / f"{sample_id}_{barcode}.HP.bed"
                command = (
                    f"modkit pileup --ref {_quote(self.settings.reference)} "
                    f"--phased --prefix {_quote(prefix)} -t {threads} "
                    f"{_quote(input_bam)} {_quote(output_file)}"
                )
                self.runner.run(command)

    def _stage_modkit_nohp(self) -> None:
        stage = "modkit_nohp"
        threads = 10
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                input_bam = sample_dir / f"{sample_id}_{barcode}.HP.bam"
                output_dir = sample_dir / f"{sample_id}_{barcode}_modkit"
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / f"{sample_id}_{barcode}.bed"
                command = (
                    f"modkit pileup --ref {_quote(self.settings.reference)} -t {threads} "
                    f"{_quote(input_bam)} {_quote(output_file)}"
                )
                self.runner.run(command)

    def _stage_tdb(self) -> None:
        stage = "tdb"
        for row in self.rows:
            sample_id = self._require(row, "sample_id", stage)
            sample_dir = self._sample_dir(sample_id)
            for barcode in self._trio_barcodes(sample_id):
                sample_tag = self._sample_tag(sample_id, barcode)
                input_vcf = sample_dir / f"{sample_id}_{barcode}_tr/medaka_to_ref.TR.vcf"
                output_tdb = sample_dir / f"{sample_id}_{barcode}_tr/{sample_id}_{barcode}.tdb"
                command = (
                    f"tdb create {_quote(input_vcf)} -o {_quote(output_tdb)} "
                    f"-s {_quote(sample_tag)}"
                )
                self.runner.run(command)
