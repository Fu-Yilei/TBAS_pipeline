"""Variant / SV / tandem-repeat ranking for TBAS trios.

This module turns annotation output into ranked candidate tables. It is pure
standard library (no pandas) so the package stays dependency-free and installs
anywhere; the annotation itself is produced by external tools (ANNOVAR,
AnnotSV) in the pipeline stages that call these functions.

The three entry points mirror the three analyses:
  - rank_snvs(): the SNV tiering system on an ANNOVAR multianno table.
  - filter_svs(): pathogenic SV filtering on an AnnotSV table.
  - rank_tandem_repeats(): trio-aware tandem-repeat comparison across the
    proband and both parents, restricted to pathogenic (STRchive) loci.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

# --- SNV tiering thresholds (from the study methods) ---
REVEL_THRESHOLD = 0.75          # high-impact if REVEL > 0.75
CADD_THRESHOLD = 20.0           # high-impact if CADD Phred > 20
AF_REPORTING_MAX = 0.001        # only report variants with 1000G ALL AF < 0.001

# A ClinVar clinical-significance label counts toward Tier 1 only when it
# reports (likely) pathogenicity — a plain "Benign" label should never promote a
# variant. This matches the pathogenic filter used in the trio analysis.
CLINVAR_PATHOGENIC = re.compile(r"pathogenic", re.IGNORECASE)

_ALLELE_SPLIT = re.compile(r"[/|]")


# --------------------------------------------------------------------------- #
# Genotype helpers
# --------------------------------------------------------------------------- #
def extract_gt(sample_field: str | None) -> str:
    """Return the GT token (e.g. ``1/1``) from a VCF sample field.

    ANNOVAR carries the raw VCF sample column through as ``GT:GQ:DP:...``; the
    genotype is everything before the first colon.
    """
    if not sample_field:
        return ""
    return sample_field.split(":", 1)[0].strip()


def _alleles(gt: str) -> list[str]:
    return [a for a in _ALLELE_SPLIT.split(gt) if a != ""]


def is_missing(gt: str) -> bool:
    al = _alleles(gt)
    return (not al) or any(a == "." for a in al)


def carries_alt(gt: str) -> bool:
    return any(a.isdigit() and int(a) >= 1 for a in _alleles(gt))


def is_hom_alt(gt: str) -> bool:
    al = _alleles(gt)
    return (
        len(al) >= 2
        and all(a.isdigit() and int(a) >= 1 for a in al)
        and len(set(al)) == 1
    )


def is_het(gt: str) -> bool:
    al = _alleles(gt)
    if len(al) < 2 or any(a == "." for a in al):
        return False
    return carries_alt(gt) and len(set(al)) > 1


def is_hom_ref(gt: str) -> bool:
    al = _alleles(gt)
    return len(al) >= 2 and all(a == "0" for a in al)


def parse_float(value: str | None) -> float | None:
    """Parse an ANNOVAR numeric cell; ``.``/empty (absent) become ``None``."""
    if value is None:
        return None
    v = value.strip()
    if v in ("", ".", "NA", "nan", "None"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parent_origin(pro: str, mom: str, dad: str) -> str | None:
    """For a proband-het variant, which parent transmitted the ALT (in trans)."""
    if not is_het(pro):
        return None
    if carries_alt(mom) and is_hom_ref(dad):
        return "mother"
    if carries_alt(dad) and is_hom_ref(mom):
        return "father"
    return None


# --------------------------------------------------------------------------- #
# SNV tiering
# --------------------------------------------------------------------------- #
SNV_OUTPUT_COLUMNS = [
    "Tier", "Chr", "Start", "End", "Ref", "Alt", "Gene",
    "Inheritance", "REVEL_score", "CADD_phred", "AF_1000g_ALL", "CLNSIG",
    "Proband_GT", "Mother_GT", "Father_GT",
]


def _snv_candidate(fields: Sequence[str], idx: dict[str, int]) -> dict | None:
    """Build a candidate record from one multianno line, or None to skip."""
    def get(name: str) -> str:
        i = idx.get(name)
        return fields[i] if i is not None and i < len(fields) else ""

    pro = extract_gt(get("proband"))
    if not carries_alt(pro):
        return None

    af = parse_float(get("af"))
    if af is not None and af >= AF_REPORTING_MAX:  # too common to report
        return None

    revel = parse_float(get("revel"))
    cadd = parse_float(get("cadd"))
    high_impact = (cadd is not None and cadd > CADD_THRESHOLD) or (
        revel is not None and revel > REVEL_THRESHOLD
    )
    clnsig = get("clnsig").strip()
    clinvar_pathogenic = bool(CLINVAR_PATHOGENIC.search(clnsig))
    mom = extract_gt(get("mother"))
    dad = extract_gt(get("father"))
    proband_hom = is_hom_alt(pro)

    # Keep only rows that can reach some tier; this bounds memory on WGS inputs.
    if not (high_impact or clinvar_pathogenic or proband_hom):
        return None

    flags = set()
    if is_hom_ref(mom) and is_hom_ref(dad):
        flags.add("de_novo")
    if proband_hom and not is_hom_alt(mom) and not is_hom_alt(dad):
        flags.add("hom_recessive")
    if is_het(pro):
        flags.add("het")

    return {
        "chr": get("chr"), "start": get("start"), "end": get("end"),
        "ref": get("ref"), "alt": get("alt"), "gene": get("gene").strip(),
        "pro": pro, "mom": mom, "dad": dad,
        "revel": revel, "cadd": cadd, "af": af, "clnsig": clnsig,
        "high_impact": high_impact, "clinvar_pathogenic": clinvar_pathogenic,
        "proband_hom": proband_hom, "flags": flags,
        "origin": parent_origin(pro, mom, dad),
    }


def _assign_tier(c: dict, compound_genes: set[str]) -> int | None:
    """Return the highest (lowest-numbered) tier a candidate qualifies for."""
    ultra_rare = c["af"] is None  # absent from 1000G entirely
    is_compound = (
        c["gene"] in compound_genes and is_het(c["pro"]) and c["origin"] is not None
    )
    tiers: list[int] = []
    # Tier 1: homozygous or de novo high-impact, or ClinVar (likely) pathogenic.
    if (("hom_recessive" in c["flags"] or "de_novo" in c["flags"]) and c["high_impact"]) \
            or c["clinvar_pathogenic"]:
        tiers.append(1)
    # Tier 2: compound heterozygous, same score thresholds.
    if is_compound and c["high_impact"]:
        tiers.append(2)
    # Tier 3: ultra-rare homozygous with no pathogenic ClinVar annotation.
    if c["proband_hom"] and ultra_rare and not c["clinvar_pathogenic"]:
        tiers.append(3)
    # Tier 4: heterozygous meeting the Tier-1 high-impact criteria.
    if is_het(c["pro"]) and c["high_impact"]:
        tiers.append(4)
    return min(tiers) if tiers else None


def _inheritance_label(c: dict, tier: int) -> str:
    parts = sorted(c["flags"])
    if tier == 2:
        parts.append(f"compound_het({c['origin']})")
    return ",".join(parts) if parts else "carrier"


def rank_snvs(
    multianno_path: str | Path,
    out_path: str | Path,
    *,
    proband_field: str = "Otherinfo13",
    mother_field: str = "Otherinfo14",
    father_field: str = "Otherinfo15",
    af_field: str = "ALL.sites.2015_08",
    revel_field: str = "REVEL_score",
    cadd_field: str = "CADD_phred",
    clnsig_field: str = "CLNSIG",
    gene_field: str = "Gene.refGeneWithVer",
) -> int:
    """Tier and rank trio SNVs from an ANNOVAR multianno table.

    Streams the (potentially multi-GB) table, keeps only tier-eligible
    candidates, detects compound-heterozygous pairs per gene, assigns each
    candidate the highest tier it qualifies for, and writes a ranked TSV sorted
    by tier then REVEL then CADD. Returns the number of ranked variants.
    """
    multianno_path = Path(multianno_path)
    with multianno_path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wanted = {
            "chr": "Chr", "start": "Start", "end": "End", "ref": "Ref",
            "alt": "Alt", "gene": gene_field, "proband": proband_field,
            "mother": mother_field, "father": father_field, "af": af_field,
            "revel": revel_field, "cadd": cadd_field, "clnsig": clnsig_field,
        }
        idx = {key: header.index(col) for key, col in wanted.items() if col in header}
        for required in ("chr", "start", "proband", "mother", "father"):
            if required not in idx:
                raise ValueError(
                    f"Multianno file missing expected column for '{required}' "
                    f"({wanted[required]}): {multianno_path}"
                )
        candidates: list[dict] = []
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            cand = _snv_candidate(fields, idx)
            if cand is not None:
                candidates.append(cand)

    # Compound-het detection: a gene with >=2 proband-het candidates transmitted
    # from different parents (one in trans from mother, one from father).
    by_gene: dict[str, set[str]] = defaultdict(set)
    for c in candidates:
        if is_het(c["pro"]) and c["gene"] and c["origin"] in ("mother", "father"):
            by_gene[c["gene"]].add(c["origin"])
    compound_genes = {g for g, origins in by_gene.items() if {"mother", "father"} <= origins}

    ranked: list[tuple[int, dict]] = []
    for c in candidates:
        tier = _assign_tier(c, compound_genes)
        if tier is not None:
            ranked.append((tier, c))

    # Sort by tier, then by descending REVEL then CADD (missing scores last).
    ranked.sort(key=lambda tc: (
        tc[0],
        -(tc[1]["revel"] if tc[1]["revel"] is not None else -1.0),
        -(tc[1]["cadd"] if tc[1]["cadd"] is not None else -1.0),
    ))

    out_path = Path(out_path)
    with out_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(SNV_OUTPUT_COLUMNS)
        for tier, c in ranked:
            writer.writerow([
                tier, c["chr"], c["start"], c["end"], c["ref"], c["alt"],
                c["gene"], _inheritance_label(c, tier),
                "" if c["revel"] is None else c["revel"],
                "" if c["cadd"] is None else c["cadd"],
                "." if c["af"] is None else c["af"], c["clnsig"],
                c["pro"], c["mom"], c["dad"],
            ])
    return len(ranked)


# --------------------------------------------------------------------------- #
# SV pathogenic filtering (AnnotSV)
# --------------------------------------------------------------------------- #
_ACMG_INT = re.compile(r"\d+")


def _is_pathogenic_acmg(value: str) -> bool:
    """True if an AnnotSV ACMG class means (likely) pathogenic.

    Handles the normalized text form (``pathogenic``/``likely_pathogenic``) and
    the numeric ACMG class (4 = likely pathogenic, 5 = pathogenic), including
    AnnotSV's ``full=N`` prefixed form.
    """
    v = (value or "").strip().lower()
    if not v or v in (".", "na"):
        return False
    if "pathogenic" in v:  # covers pathogenic and likely_pathogenic
        return True
    return any(int(n) in (4, 5) for n in _ACMG_INT.findall(v))


def filter_svs(annotsv_tsv: str | Path, out_path: str | Path) -> int:
    """Keep AnnotSV rows classified (likely) pathogenic; write them out.

    Prefers the ``ACMG_class_norm`` column and falls back to ``ACMG_class``.
    Returns the number of pathogenic SVs kept.
    """
    annotsv_tsv = Path(annotsv_tsv)
    with annotsv_tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty AnnotSV file: {annotsv_tsv}")
        acmg_col = next(
            (c for c in ("ACMG_class_norm", "ACMG_class") if c in reader.fieldnames),
            None,
        )
        if acmg_col is None:
            raise ValueError(
                f"AnnotSV file has no ACMG_class(_norm) column: {annotsv_tsv}"
            )
        kept = [row for row in reader if _is_pathogenic_acmg(row.get(acmg_col, ""))]
        fieldnames = reader.fieldnames

    out_path = Path(out_path)
    with out_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)
    return len(kept)


# --------------------------------------------------------------------------- #
# Trio tandem-repeat comparison
# --------------------------------------------------------------------------- #
def _load_bed_intervals(bed_path: str | Path) -> dict[str, list[tuple[int, int, str]]]:
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with Path(bed_path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            name = f[3] if len(f) > 3 else ""
            intervals[f[0]].append((int(f[1]), int(f[2]), name))
    return intervals


def _overlapping_name(
    intervals: dict[str, list[tuple[int, int, str]]],
    chrom: str, qstart: int, qend: int,
) -> str | None:
    """Name of the first interval on ``chrom`` that overlaps ``[qstart, qend]``.

    Overlap (not point containment) matters because the tandem-repeat catalog
    used by medaka and STRchive anchor the same locus at slightly different
    coordinates.
    """
    for start, end, name in intervals.get(chrom, ()):  # STRchive is small; linear is fine
        if start <= qend and qstart <= end:
            return name or f"{chrom}:{start}-{end}"
    return None


def _parse_tr_vcf(vcf_path: str | Path):
    """Yield (chrom, pos, ref_end, allele_lengths) per record of a sample VCF.

    Allele length is the length of the called allele sequence (REF for index 0,
    the matching ALT otherwise) — what a repeat expansion changes. ``ref_end``
    is the last reference base the record spans, used for interval overlap.
    """
    with Path(vcf_path).open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            gt = extract_gt(f[9])
            seqs = [ref] + (alt.split(",") if alt not in (".", "") else [])
            lengths: list[int] = []
            for a in _alleles(gt):
                if a.isdigit() and int(a) < len(seqs):
                    seq = seqs[int(a)]
                    lengths.append(0 if seq in (".", "*") else len(seq))
            if lengths:
                yield chrom, pos, pos + len(ref) - 1, lengths


def sample_allele_lengths(vcf_path: str | Path) -> dict[tuple[str, int], list[int]]:
    """Map (chrom, pos) -> the sample's called tandem-repeat allele lengths."""
    return {(c, p): lengths for c, p, _end, lengths in _parse_tr_vcf(vcf_path)}


TR_OUTPUT_COLUMNS = [
    "Chr", "Pos", "STRchive_locus", "Proband_allele_lengths",
    "Mother_allele_lengths", "Father_allele_lengths",
    "Proband_distinct_lengths", "Parents_compared", "Expansion_vs_parents",
]


def _distinct_lengths(proband: list[int], parents: list[list[int]],
                      tolerance: int) -> list[int]:
    """Proband allele lengths differing from every allele of every parent."""
    return [
        p for p in proband
        if all(all(abs(p - x) > tolerance for x in lens) for lens in parents)
    ]


def rank_tandem_repeats(
    proband_vcf: str | Path,
    mother_vcf: str | Path,
    father_vcf: str | Path,
    out_path: str | Path,
    *,
    strchive_bed: str | Path | None = None,
    length_tolerance: int = 0,
) -> int:
    """Report pathogenic tandem-repeat loci where the proband has an allele
    distinct from both parents.

    Restricts to STRchive (pathogenic) loci when a BED is given, then at each
    locus compares the proband's tandem-repeat allele lengths against the
    parents. A locus is reported when the proband carries at least one allele
    length that matches neither parent (within ``length_tolerance``). A parent
    with no confident call at the locus cannot be compared, so the report notes
    which parents were actually compared (``Parents_compared``); a fully
    confident "distinct from both parents" call is one where that column reads
    ``mother,father``. ``Expansion_vs_parents`` flags when the proband's longest
    allele exceeds every compared parent's. Returns the number of reported loci.
    """
    pro_end: dict[tuple[str, int], int] = {}
    pro: dict[tuple[str, int], list[int]] = {}
    for chrom, pos, ref_end, lengths in _parse_tr_vcf(proband_vcf):
        pro[(chrom, pos)] = lengths
        pro_end[(chrom, pos)] = ref_end
    mom = sample_allele_lengths(mother_vcf)
    dad = sample_allele_lengths(father_vcf)
    strchive = _load_bed_intervals(strchive_bed) if strchive_bed else None

    rows = []
    for (chrom, pos), p_lens in sorted(pro.items()):
        m_lens = mom.get((chrom, pos))
        f_lens = dad.get((chrom, pos))
        compared: list[tuple[str, list[int]]] = []
        if m_lens:
            compared.append(("mother", m_lens))
        if f_lens:
            compared.append(("father", f_lens))
        if not compared:
            continue  # no parent called -> nothing to compare against
        locus = None
        if strchive is not None:
            locus = _overlapping_name(strchive, chrom, pos, pro_end[(chrom, pos)])
            if locus is None:
                continue  # not a pathogenic (STRchive) locus
        distinct = _distinct_lengths(p_lens, [lens for _, lens in compared],
                                     length_tolerance)
        if not distinct:
            continue
        expansion = all(max(p_lens) > max(lens) for _, lens in compared)
        rows.append([
            chrom, pos, locus or "",
            ",".join(map(str, p_lens)),
            ",".join(map(str, m_lens)) if m_lens else ".",
            ",".join(map(str, f_lens)) if f_lens else ".",
            ",".join(map(str, distinct)),
            ",".join(name for name, _ in compared),
            "yes" if expansion else "no",
        ])

    out_path = Path(out_path)
    with out_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(TR_OUTPUT_COLUMNS)
        writer.writerows(rows)
    return len(rows)
