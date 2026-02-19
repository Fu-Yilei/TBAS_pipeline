# TBAS_pipeline
We present Trio-barcoded ONT Adaptive Sampling (TBAS), a cost-efficient long-read sequencing strategy combining sample barcoding and adaptive enrichment to sequence rare-disease trios on a single PromethION/P2 flow cell. TBAS achieved near-complete variant phasing and detection of small variants, structural variants, and tandem repeats with high accuracy and 77% potential solve rate. This scalable approach retains methylation data and enables clinically relevant, phenotype-guided long-read diagnostics at a fraction of current costs.


## Overview

- The software pipeline entrypoint is `tbas-pipeline`, implemented in `tbas_pipeline/`.
- The software version runs tools directly (no `sbatch` generation).
- The notebook `TBAS_pipeline_slurm.ipynb` is kept as-is for historical/reference workflow.
- The `analysis/` folder contains downstream analysis notebooks organized by topic (SNV/Trio calling, coverage, methylation, tandem repeats, and variant counting).
- Benchmarking materials and example HG002 results are provided under `benchmarking/`, with data available at the Zenodo record below.

## Requirements

- Python 3.9+.
- Standard long-read analysis command-line tools used by the pipeline:
  `samtools`, `minimap2`, `sniffles`, `mosdepth`, `run_clair3.sh`,
  `kanpig`, `bedtools`, `bgzip`, `bcftools`, `whatshap`, `medaka`, `modkit`, `tdb`.

## Getting started

1. Install the local package:

```bash
pip install -e .
```

2. Create a manifest CSV/TSV with at least:
- `sample_id` (example: `4_6_Gregor_Trio`)
- `bed_file` (example: `Epilepsy`, `CMRG`, `WGS`, or a direct BED path)
- `proband_gender` (used by medaka stages, example: `female`)

Optional columns:
- `calls_bam` (explicit path to `calls*.bam`; if omitted, the pipeline searches under `<output_folder>/<sample_id>/calls*.bam`)
- `read_group_prefix` (skip BAM-header inference during demultiplex stage)
- `tr_bed_file` (TR catalog BED for `medaka_local`; if omitted, the pipeline derives this from `bed_file` via built-in mapping)

3. Run a dry run starting from demultiplexing:

```bash
tbas-pipeline \
  --manifest manifest.csv \
  --output-folder /stornext/snfs170/next-gen/scratch/Yilei/projects/adaptive_sampling/GREGoR_adaptive_sampling \
  --stages demultiplex \
  --dry-run
```

4. Run full pipeline:

```bash
tbas-pipeline \
  --manifest manifest.csv \
  --output-folder /stornext/snfs170/next-gen/scratch/Yilei/projects/adaptive_sampling/GREGoR_adaptive_sampling
```

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
