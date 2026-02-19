from __future__ import annotations

import csv
from pathlib import Path

import pytest

from tbas_pipeline.pipeline import (
    PipelineSettings,
    TBASPipeline,
    parse_readgroup_prefix_from_header,
    parse_sample_barcodes,
)


class RecordingRunner:
    def __init__(self) -> None:
        self.commands: list[str] = []

    def run(self, command: str) -> None:
        self.commands.append(command)


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_parse_sample_barcodes() -> None:
    assert parse_sample_barcodes("4_6_Gregor_Trio") == [
        "barcode04",
        "barcode05",
        "barcode06",
        "unassigned",
    ]
    assert parse_sample_barcodes("43-45") == [
        "barcode43",
        "barcode44",
        "barcode45",
        "unassigned",
    ]
    with pytest.raises(ValueError):
        parse_sample_barcodes("invalid_sample")


def test_parse_readgroup_prefix_from_header() -> None:
    header = (
        "@HD\tVN:1.6\n"
        "@RG\tID:runprefix_dna_r10_sup_SQK-NBD114-96_barcode04\n"
        "@RG\tID:runprefix_dna_r10_sup\n"
    )
    prefix = parse_readgroup_prefix_from_header(header, "SQK-NBD114-96")
    assert prefix == "runprefix_dna_r10_sup"


def test_pipeline_dry_run_starts_with_demultiplex(tmp_path: Path) -> None:
    sample_id = "4_6_Gregor_Trio"
    output_root = tmp_path / "GREGoR_adaptive_sampling"
    sample_dir = output_root / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    calls_bam = sample_dir / "calls_2025-07-25_T21-05-41.bam"
    calls_bam.touch()

    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": sample_id,
                "bed_file": "Epilepsy",
                "proband_gender": "female",
                "calls_bam": str(calls_bam),
                "read_group_prefix": "2bd2bfe2-0259-4d7e-89fc-1819065db57c_dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
            }
        ],
    )

    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest, output_folder=output_root)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["demultiplex", "fastq_extract", "minimap2", "bam_sort"])

    assert len(runner.commands) == 12
    assert runner.commands[0].startswith("samtools view -b@")
    assert "barcode04" in runner.commands[0]
    assert "samtools fastq" in runner.commands[3]
    assert "minimap2" in runner.commands[6]
    assert "samtools sort" in runner.commands[9]


def test_unknown_stage_raises(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "Epilepsy",
                "proband_gender": "female",
            }
        ],
    )
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=RecordingRunner())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        pipeline.run(["not_a_stage"])


def test_kanpig_pileup_adds_tabix_index(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "CMRG",
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["kanpig_pileup"])

    assert len(runner.commands) == 6
    assert runner.commands[0].startswith("kanpig plup -b ")
    assert runner.commands[1].startswith("tabix -f -s 1 -b 2 -e 2 ")
    assert sum(cmd.startswith("tabix -f -s 1 -b 2 -e 2 ") for cmd in runner.commands) == 3


def test_whatshap_phase_adds_tabix_index(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "CMRG",
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["whatshap_single_sample_local_phasing"])

    assert len(runner.commands) == 6
    assert runner.commands[0].startswith("whatshap phase --reference ")
    assert runner.commands[1].startswith("tabix -f -p vcf ")
    assert sum(cmd.startswith("tabix -f -p vcf ") for cmd in runner.commands) == 3


def test_whatshap_haplotag_adds_bam_index(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "CMRG",
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["whatshap_haplotag"])

    assert len(runner.commands) == 6
    assert runner.commands[0].startswith("whatshap haplotag --ignore-read-groups ")
    assert runner.commands[1].startswith("samtools index -@ 3 ")
    assert sum(cmd.startswith("samtools index -@ 3 ") for cmd in runner.commands) == 3


def test_modkit_commands_match_installed_cli(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "CMRG",
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["modkit", "modkit_nohp"])

    assert len(runner.commands) == 6
    modkit_cmds = runner.commands[:3]
    modkit_nohp_cmds = runner.commands[3:]

    for cmd in modkit_cmds:
        assert cmd.startswith("modkit pileup --ref ")
        assert "--phased" in cmd
        assert "--prefix " in cmd
        assert "--preset " not in cmd
        assert "--partition-tag " not in cmd

    for cmd in modkit_nohp_cmds:
        assert cmd.startswith("modkit pileup --ref ")
        assert "--preset " not in cmd


def test_medaka_local_prefers_tr_bed_file(tmp_path: Path) -> None:
    output_root = tmp_path / "out"
    manifest = tmp_path / "manifest.csv"
    custom_tr_bed = tmp_path / "custom_tr.bed"
    custom_tr_bed.write_text("chr1\t1\t2\tTR1\n")
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "bed_file": "CMRG",
                "tr_bed_file": str(custom_tr_bed),
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest, output_folder=output_root)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    pipeline.run(["medaka_local"])

    assert len(runner.commands) == 3
    assert str(custom_tr_bed) in runner.commands[0]
    assert "analysis/tr_regions/adotto_catalog.hg38.lite.CMRG.bed" not in runner.commands[0]


def test_medaka_local_requires_tr_bed_or_bed(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.csv"
    _write_manifest(
        manifest,
        [
            {
                "sample_id": "4_6_Gregor_Trio",
                "proband_gender": "female",
            }
        ],
    )
    runner = RecordingRunner()
    settings = PipelineSettings(manifest=manifest)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="tr_bed_file.*bed_file"):
        pipeline.run(["medaka_local"])
