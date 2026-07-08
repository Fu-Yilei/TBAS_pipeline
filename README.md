# TBAS_pipeline
We present Trio-barcoded ONT Adaptive Sampling (TBAS), a cost-efficient long-read sequencing strategy combining sample barcoding and adaptive enrichment to sequence rare-disease trios on a single PromethION/P2 flow cell. TBAS achieved near-complete variant phasing and detection of small variants, structural variants, and tandem repeats with high accuracy and 77% potential solve rate. This scalable approach retains methylation data and enables clinically relevant, phenotype-guided long-read diagnostics at a fraction of current costs.


## Overview

- The software pipeline entrypoint is `tbas-pipeline`, implemented in `tbas_pipeline/`.
- The software version runs tools directly (no `sbatch` generation).
- The notebook `TBAS_pipeline_slurm.ipynb` is kept as-is for historical/reference workflow.
- The `analysis/` folder contains downstream analysis notebooks organized by topic (SNV/Trio calling, coverage, methylation, tandem repeats, and variant counting).
- Benchmarking materials and example HG002 results are provided under `benchmarking/`, with data available at the Zenodo record below.

## Analytical novelty

TBAS is not only a cheaper way to sequence trios — its analytical approach is
co-designed with long-read sequencing (LRS) to extract signal that short-read
trio pipelines cannot. Three points are novel:

1. **Automated, phenotype-guided trio variant ranking.** Small and structural
   variants are called per trio member and genotyped *trio-aware* (`kanpig
   trio`), then annotated (ANNOVAR functional/population annotation) and ranked
   by combining that annotation with the trio inheritance pattern — de novo,
   recessive/compound-heterozygous, and phenotype (HPO) fit — so candidate
   diagnoses surface automatically rather than by manual review. Because reads
   are natively phased on a single flow cell, compound-heterozygous and
   parent-of-origin calls are made directly instead of inferred.

2. **DNA methylation from LRS as an independent marker of pathogenic
   expansions — including from adaptive-sampling *rejected* reads.** TBAS
   preserves native 5mC/5hmC modification tags end to end (methylation tags are
   carried through read extraction and emitted as phased and unphased pileups by
   `modkit`). Locus hypermethylation is then used as an orthogonal signal for
   repeat-expansion disorders even where the expansion is not fully spanned. To
   our knowledge this is the first approach to combine LRS methylation with
   adaptive sampling to flag pathogenic expansions such as **NAXE** and
   **DIP2B**, notably by recovering methylation signal from reads that adaptive
   sampling *rejected* (unblocked/off-target) rather than discarding them
   (worked examples in `analysis/methylation/`).

3. **Trio-aware tandem-repeat length comparison.** Long reads span full repeat
   alleles, so `medaka tandem` resolves per-sample allele lengths and `tdb`
   collates the proband and both parents into one database. This lets the
   pipeline compare repeat length *across the trio* directly and flag a proband
   allele that is expanded relative to the parents (de novo expansion or further
   somatic/germline expansion) — a comparison prior tandem-repeat pipelines,
   built for single samples, did not make. See the `patho_ATN1_FGF14` example
   (ATN1/DRPLA `CAG`, FGF14/SCA27B `GAA`).

Together these turn TBAS's LRS strengths — native methylation, read phasing,
full-length repeat spanning, and the trio structure itself — into diagnostic
evidence, at a fraction of the cost of running trios on separate flow cells.

## Requirements

- Python 3.10+.
- Standard long-read analysis command-line tools used by the pipeline:
  `samtools`, `minimap2`, `sniffles`, `mosdepth`, `run_clair3.sh`,
  `kanpig`, `bedtools`, `bgzip`, `bcftools`, `whatshap`, `medaka`, `modkit`, `tdb`.

The easiest way to get all of these is the bundled conda environment (see
below), which installs the whole toolchain in one command.

## Getting started

1. Install the external tools and the package into a single self-contained
   environment. With micromamba (or swap `micromamba` for `mamba`/`conda`):

```bash
micromamba env create -f environment.yml
micromamba activate tbas
pip install -e .
```

This one environment provides the entire toolchain (samtools, minimap2,
sniffles, mosdepth, Clair3, kanpig, bedtools, bcftools, whatshap, medaka,
modkit, tdb) — no separate per-tool envs or from-source builds. `tdb` installs
from GitHub via pip as declared in `environment.yml`. For a byte-for-byte
reproducible environment, create from the pinned `environment.lock.yml` instead.

Or, if the command-line tools are already on your `PATH`, just:

```bash
pip install -e .
```

2. Verify the toolchain is complete before running:

```bash
tbas-pipeline --manifest example_data/test_subset_chr22/manifest_example.csv --check-deps
```

This lists the tools required by the selected stages and reports any that are
missing from `PATH`. A real run also performs this check automatically and
aborts early if a tool is missing (pass `--skip-dep-check` to override).

3. Create a manifest CSV/TSV with at least:
- `sample_id` (example: `4_6_Gregor_Trio`)
- `bed_file` (example: `Epilepsy`, `CMRG`, `WGS`, or a direct BED path)
- `proband_gender` (used by medaka stages, example: `female`)

Optional columns:
- `calls_bam` (explicit path to `calls*.bam`; if omitted, the pipeline searches under `<output_folder>/<sample_id>/calls*.bam`)
- `read_group_prefix` (skip BAM-header inference during demultiplex stage)
- `tr_bed_file` (TR catalog BED for `medaka_local`; if omitted, the pipeline derives this from `bed_file` via built-in mapping)

4. Run a dry run starting from demultiplexing:

```bash
tbas-pipeline \
  --manifest example_data/test_subset_chr22/manifest_example.csv \
  --output-folder demo_output \
  --stages demultiplex \
  --dry-run
```

5. Run full pipeline:

```bash
tbas-pipeline \
  --manifest example_data/test_subset_chr22/manifest_example.csv \
  --output-folder demo_output
```

See `OPTIMIZATION_PLAN.md` for a roadmap toward a fully push-button, clinically
deployable install (containers, pinned environments, resumability).

## Example data

Two real-data bundles under `example_data/` (subsets of the `4_6_Gregor_Trio`
sample) let you smoke-test the pipeline end to end:

- `example_data/test_subset_chr22/` — CMRG genes on `chr22:42.6-42.74 Mb`.
  Exercises stages 1-16 fully; `medaka_patho` is skipped because this window
  contains no known pathogenic STR loci.
- `example_data/patho_ATN1_FGF14/` — windows around **ATN1** (DRPLA, `CAG`) and
  **FGF14** (SCA27B, `GAA`). These loci are in the pathogenic catalog, so this
  bundle also exercises `medaka_patho` and `tdb` on real disease repeats.

Both runs need a GRCh38 reference, a Clair3 model directory, and the pathogenic
TR BED, passed via `--reference`, `--clair3-model`, and `--adotto-pheno`.

## Pipeline stages

Default stage order:

1. `demultiplex`
2. `fastq_extract`
3. `minimap2`
4. `bam_sort`
5. `sniffles_global`
6. `bam_mosdepth`
7. `clair3_local`
8. `kanpig_pileup`
9. `kanpig_gt`
10. `kanpig_trio`
11. `whatshap_single_sample_local_phasing`
12. `whatshap_haplotag`
13. `medaka_local`
14. `medaka_patho`
15. `modkit`
16. `modkit_nohp`
17. `tdb`

You can run a subset with `--stages stage1,stage2,...`.

## Notebook reference

- `TBAS_pipeline_slurm.ipynb` is preserved and unchanged.
- Use it as reference for the original SLURM-oriented workflow.

## Downstream analysis notebooks

The `analysis/` directory contains topic-focused notebooks and resources:

- SNV/Trio calling
	- `analysis/snv_trio_analysis.ipynb`: downstream analysis notebook for small variants and trio evaluation.

- Coverage
	- `analysis/covearge/mosdepth_analysis.ipynb`: coverage analysis with mosdepth.
	- `analysis/covearge/mosdepth_summary_totals.csv`: example summary output.

- Methylation
	- `analysis/methylation/methylation_analysis.ipynb`: methylation analysis notebook.
	- `analysis/methylation/*/`: per-sample directories with supporting files.

- Tandem repeats (TR)
	- `analysis/tr/adaptive_sampling_tr_analysis.ipynb`: TR analysis on adaptive sampling data.
	- `analysis/tr/tr_regions/`: TR region catalogs used by the analysis, including:
		- `adotto_catalog.hg38.lite.*.bed`
		- `strchive/STRchive-disease-loci.hg38.bed`

- Variant counts
	- `analysis/variant_counts/count_variants_by_barcode.sh`
	- `analysis/variant_counts/count_variants_by_barcode_parallel.sh`
	- `analysis/variant_counts/variant_counts.tsv`

## Benchmarking

- Folder: `benchmarking/`
	- `adaptive_methylation_benchmark.ipynb`: methylation benchmarking notebook.
	- `benchmarking/hg002_variant_calls/`: example HG002 variant call files and indices for reference.

Benchmarking data DOI: https://doi.org/10.5281/zenodo.17398577


## Data availability and citation

If you use these materials, please reference the Zenodo record:

- TBAS benchmarking data: https://doi.org/10.5281/zenodo.17398577

## License

This project is distributed under the terms of the license in `LICENSE`.
