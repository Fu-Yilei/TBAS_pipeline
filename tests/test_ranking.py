from __future__ import annotations

from pathlib import Path

from tbas_pipeline import ranking

# --------------------------------------------------------------------------- #
# Genotype helpers
# --------------------------------------------------------------------------- #
def test_genotype_helpers() -> None:
    assert ranking.extract_gt("1/1:110:33:0,33") == "1/1"
    assert ranking.is_hom_alt("1|1")
    assert ranking.is_hom_alt("2/2")
    assert not ranking.is_hom_alt("0/1")
    assert ranking.is_het("0|1")
    assert ranking.is_het("1/2")
    assert not ranking.is_het("1/1")
    assert ranking.is_hom_ref("0/0")
    assert ranking.carries_alt("0/1") and not ranking.carries_alt("0/0")
    assert ranking.is_missing("./.")
    assert ranking.parse_float(".") is None
    assert ranking.parse_float("0.75") == 0.75


# --------------------------------------------------------------------------- #
# SNV tiering
# --------------------------------------------------------------------------- #
_MULTIANNO_COLS = [
    "Chr", "Start", "End", "Ref", "Alt", "Gene.refGeneWithVer",
    "REVEL_score", "CADD_phred", "ALL.sites.2015_08", "CLNSIG",
    "Otherinfo13", "Otherinfo14", "Otherinfo15",
]


def _write_multianno(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w") as fh:
        fh.write("\t".join(_MULTIANNO_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, ".")) for c in _MULTIANNO_COLS) + "\n")


def _variant(chrom, pos, gene, revel, cadd, af, clnsig, pro, mom, dad) -> dict:
    return {
        "Chr": chrom, "Start": pos, "End": pos, "Ref": "A", "Alt": "T",
        "Gene.refGeneWithVer": gene, "REVEL_score": revel, "CADD_phred": cadd,
        "ALL.sites.2015_08": af, "CLNSIG": clnsig,
        "Otherinfo13": pro, "Otherinfo14": mom, "Otherinfo15": dad,
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def test_snv_tiering_assigns_expected_tiers(tmp_path: Path) -> None:
    manno = tmp_path / "trio.hg38_multianno.txt"
    _write_multianno(manno, [
        # Tier 1: de novo, high CADD
        _variant("chr1", 100, "GENE_DN", ".", "25", ".", ".", "0/1", "0/0", "0/0"),
        # Tier 1: homozygous recessive, high REVEL
        _variant("chr1", 200, "GENE_HOM", "0.80", "5", ".", ".", "1/1", "0/1", "0/1"),
        # Tier 1: ClinVar pathogenic, low scores
        _variant("chr1", 300, "GENE_CV", "0.1", "5", ".", "Pathogenic", "0/1", "0/0", "0/1"),
        # Tier 2: compound het (two vars in GENE_CH, one per parent), high CADD
        _variant("chr2", 400, "GENE_CH", ".", "25", ".", ".", "0/1", "0/1", "0/0"),
        _variant("chr2", 500, "GENE_CH", ".", "25", ".", ".", "0/1", "0/0", "0/1"),
        # Tier 3: ultra-rare homozygous, no ClinVar, low scores
        _variant("chr3", 600, "GENE_UR", "0.1", "5", ".", ".", "1/1", "0/1", "0/1"),
        # Tier 4: heterozygous high-impact (both parents carriers => not de novo/compound)
        _variant("chr4", 700, "GENE_HET", ".", "25", ".", ".", "0/1", "0/1", "0/1"),
        # Dropped: too common (AF >= 0.001)
        _variant("chr5", 800, "GENE_COMMON", ".", "40", "0.5", "Pathogenic", "0/1", "0/0", "0/0"),
        # Dropped: benign ClinVar, low impact, heterozygous
        _variant("chr6", 900, "GENE_BEN", "0.1", "5", ".", "Benign", "0/1", "0/1", "0/1"),
    ])
    out = tmp_path / "ranked.tsv"
    n = ranking.rank_snvs(manno, out)
    rows = _read_tsv(out)
    by_gene = {r["Gene"]: r["Tier"] for r in rows}

    assert by_gene.get("GENE_DN") == "1"
    assert by_gene.get("GENE_HOM") == "1"
    assert by_gene.get("GENE_CV") == "1"
    assert by_gene.get("GENE_CH") == "2"
    assert by_gene.get("GENE_UR") == "3"
    assert by_gene.get("GENE_HET") == "4"
    assert "GENE_COMMON" not in by_gene   # AF filter
    assert "GENE_BEN" not in by_gene      # benign ClinVar, not high-impact
    # output is sorted by tier ascending
    tiers = [int(r["Tier"]) for r in rows]
    assert tiers == sorted(tiers)
    assert n == len(rows)


def test_snv_common_variant_dropped_even_if_pathogenic(tmp_path: Path) -> None:
    manno = tmp_path / "m.hg38_multianno.txt"
    _write_multianno(manno, [
        _variant("chr1", 100, "G", ".", "40", "0.02", "Pathogenic", "1/1", "0/1", "0/1"),
    ])
    out = tmp_path / "o.tsv"
    assert ranking.rank_snvs(manno, out) == 0


# --------------------------------------------------------------------------- #
# SV filtering
# --------------------------------------------------------------------------- #
def test_filter_svs_keeps_pathogenic(tmp_path: Path) -> None:
    tsv = tmp_path / "annotsv.tsv"
    tsv.write_text(
        "SV_chrom\tSV_start\tACMG_class\n"
        "chr1\t100\tfull=5\n"        # pathogenic
        "chr2\t200\t3\n"             # VUS
        "chr3\t300\tlikely_pathogenic\n"
        "chr4\t400\tNA\n"
    )
    out = tmp_path / "path.tsv"
    assert ranking.filter_svs(tsv, out) == 2
    kept = _read_tsv(out)
    assert {r["SV_chrom"] for r in kept} == {"chr1", "chr3"}


def test_filter_svs_prefers_norm_column(tmp_path: Path) -> None:
    tsv = tmp_path / "a.tsv"
    tsv.write_text(
        "SV_chrom\tACMG_class\tACMG_class_norm\n"
        "chr1\t3\tpathogenic\n"      # norm says pathogenic -> keep
        "chr2\t5\tbenign\n"          # norm says benign -> drop
    )
    out = tmp_path / "o.tsv"
    assert ranking.filter_svs(tsv, out) == 1
    assert _read_tsv(out)[0]["SV_chrom"] == "chr1"


# --------------------------------------------------------------------------- #
# Trio tandem-repeat comparison
# --------------------------------------------------------------------------- #
def _write_tr_vcf(path: Path, records: list[tuple[str, int, str, str, str]]) -> None:
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt, gt in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n")


def test_tr_reports_proband_distinct_allele(tmp_path: Path) -> None:
    # Proband carries a 6bp allele neither parent has; a second locus is inherited.
    pro = tmp_path / "pro.vcf"
    mom = tmp_path / "mom.vcf"
    dad = tmp_path / "dad.vcf"
    _write_tr_vcf(pro, [
        ("chr1", 100, "AAAA", "AAAAAA,AA", "1/2"),   # lengths 6,2
        ("chr2", 500, "AAAA", "AA", "0/1"),          # lengths 4,2 (inherited)
    ])
    _write_tr_vcf(mom, [
        ("chr1", 100, "AAAA", ".", "0/0"),           # lengths 4,4
        ("chr2", 500, "AAAA", "AA", "0/1"),          # lengths 4,2
    ])
    _write_tr_vcf(dad, [
        ("chr1", 100, "AAAA", "AA", "0/1"),          # lengths 4,2
        ("chr2", 500, "AAAA", "AA", "0/1"),          # lengths 4,2
    ])
    out = tmp_path / "tr.tsv"
    n = ranking.rank_tandem_repeats(pro, mom, dad, out)
    rows = _read_tsv(out)
    assert n == 1
    assert rows[0]["Chr"] == "chr1" and rows[0]["Pos"] == "100"
    assert rows[0]["Proband_distinct_lengths"] == "6"
    assert rows[0]["Expansion_vs_parents"] == "yes"
    assert rows[0]["Parents_compared"] == "mother,father"


def test_tr_strchive_restriction(tmp_path: Path) -> None:
    pro = tmp_path / "pro.vcf"
    mom = tmp_path / "mom.vcf"
    dad = tmp_path / "dad.vcf"
    for p, gt in ((pro, "1/1"), (mom, "0/0"), (dad, "0/0")):
        _write_tr_vcf(p, [("chr1", 100, "AAAA", "AAAAAAAA", gt)])
    # Proband hom-long, parents hom-ref -> distinct, but only report if in STRchive.
    strchive = tmp_path / "strchive.bed"
    strchive.write_text("chr1\t90\t120\tDISEASE_LOCUS\n")
    out = tmp_path / "tr.tsv"
    assert ranking.rank_tandem_repeats(pro, mom, dad, out, strchive_bed=strchive) == 1
    assert _read_tsv(out)[0]["STRchive_locus"] == "DISEASE_LOCUS"
    # A STRchive BED that does not cover the locus -> nothing reported.
    strchive2 = tmp_path / "s2.bed"
    strchive2.write_text("chr9\t1\t2\tOTHER\n")
    out2 = tmp_path / "tr2.tsv"
    assert ranking.rank_tandem_repeats(pro, mom, dad, out2, strchive_bed=strchive2) == 0
