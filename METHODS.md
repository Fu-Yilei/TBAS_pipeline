# Annotation and ranking methods

This documents the annotation databases and the variant/SV/tandem-repeat ranking
implemented by the `annovar`, `rank_snv`, `annotsv`, `rank_sv`, and `rank_tr`
pipeline stages (`tbas_pipeline/ranking.py`).

## Installing the annotation resources

The ranking stages need three external resources beyond the base toolchain:

- **AnnotSV** — installed by `environment.yml` (bioconda). Its annotation files
  ship with the conda package; if a lookup complains, run the package's
  `INSTALL_annotations.sh` once.
- **STRchive** BED — already shipped in this repo
  (`analysis/tr_regions/strchive/STRchive-disease-loci.hg38.bed`); nothing to
  install.
- **ANNOVAR + its databases** — ANNOVAR is registration-gated and not on conda.
  Register and download it from
  <https://www.openbioinformatics.org/annovar/annovar_download_form.php>, put its
  scripts on PATH, then fetch the exact GRCh38 databases below with the helper:

  ```bash
  scripts/download_annovar_db.sh annovar/humandb
  # then run the pipeline with:  --annovar-humandb annovar/humandb
  ```

  (`gnomad41_genome` is tens of GB — allow disk and time.) The `rank_snv`
  scoring uses only 1000G AF, REVEL, CADD, and ClinVar, so the SNV tiers are
  unaffected if gnomAD is temporarily absent, but the full protocol is what the
  annotation stage runs.

## Small-variant annotation (ANNOVAR)

ANNOVAR (version 2025-03-02) is run on the family (trio) small-variant VCF. The
following GRCh38 reference annotation databases are used:

- **refGene** (20211019) — all annotated transcripts in RefSeq gene.
- **avsnp151** (20170929) — dbSNP151.
- **dbnsfp47a** (20240525) — predictive scores from dbNSFP (v3.0a), incl. REVEL
  and CADD.
- **clinvar_20250721** (20250721) — ClinVar variant clinical significance (CLNSIG).
- **1000g2015aug** (20150824) — 1000 Genomes allele frequency, as
  `ALL.sites.2015_08`, `AFR.sites.2015_08`, `AMR.sites.2015_08`,
  `EAS.sites.2015_08`, `EUR.sites.2015_08`, `SAS.sites.2015_08` (overall and per
  continental-ancestry cluster).
- **gnomad41_genome** (20240602) — gnomAD genome collection (v4.1).

The exact `-protocol`/`-operation` strings are defined as `ANNOVAR_PROTOCOL` /
`ANNOVAR_OPERATION` in `tbas_pipeline/pipeline.py`. The trio VCF is produced by
the `merge_trio_snv` stage (bcftools merge of the three per-member Clair3 VCFs,
proband first so the proband/mother/father genotypes land in the
`Otherinfo13/14/15` columns).

## SNV tiering (`rank_snv`)

A tiering system filters variants by allele frequency, REVEL score, and CADD
Phred. Minimum reporting threshold: **1000 Genomes ALL AF < 0.001** (variants
absent from 1000G pass, being rarer than any threshold). Within the reported
set, SNVs are ranked by REVEL and CADD. A variant is **high-impact** if
`CADD_phred > 20` **or** `REVEL_score > 0.75`.

- **Tier 1** — homozygous and de novo SNVs that are high-impact, plus any SNV
  carrying a ClinVar (likely) pathogenic significance label.
- **Tier 2** — compound-heterozygous SNVs meeting the high-impact thresholds
  (two proband-het variants in one gene, transmitted in trans — one from each
  parent).
- **Tier 3** — ultra-rare homozygous SNVs (absent from 1000G) lacking a
  pathogenic ClinVar annotation.
- **Tier 4** — heterozygous SNVs meeting the Tier-1 high-impact criteria.

Each variant is reported at the highest (lowest-numbered) tier it qualifies for.
Note on ClinVar: the Tier-1 ClinVar clause matches **(likely) pathogenic**
significance (`CLINVAR_PATHOGENIC` in `ranking.py`), not any label — a plain
`Benign` annotation does not promote a variant. Output:
`<sample>.snv.ranked.tsv`.

## Structural-variant filtering (AnnotSV, `annotsv` + `rank_sv`)

AnnotSV v3.5 annotates the structural variants (`AnnotSV -SVinputFile
{input.vcf}`; the pipeline also passes `-outputDir`/`-genomeBuild`). Pathogenic
candidate filtering keeps rows whose `ACMG_class_norm` (falling back to
`ACMG_class`) is **pathogenic** or **likely_pathogenic** (numeric ACMG classes 5
and 4, including AnnotSV's `full=N` form). Output:
`<sample>_<barcode>.sv.pathogenic.tsv`.

## Tandem-repeat trio comparison (STRchive, `rank_tr`)

STRchive is used to intersect the tandem-repeat loci and obtain the pathogenic
repeats (overlap is computed against each medaka TR record's reference span, so
catalog/STRchive coordinate offsets do not drop loci). Loci where the **proband
carries at least one tandem-repeat allele distinct from both parents** are
reported, with the compared parents noted and an expansion flag when the
proband's longest allele exceeds every compared parent's. A parent with no
confident call at a locus is annotated (`Parents_compared`) rather than assumed.
Output: `<sample>.tr.trio_distinct.tsv`.
