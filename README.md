# TBAS_pipeline
We present Trio-barcoded ONT Adaptive Sampling (TBAS), a cost-efficient long-read sequencing strategy combining sample barcoding and adaptive enrichment to sequence rare-disease trios on a single PromethION/P2 flow cell. TBAS achieved near-complete variant phasing and detection of small variants, structural variants, and tandem repeats with high accuracy and 77% potential solve rate. This scalable approach retains methylation data and enables clinically relevant, phenotype-guided long-read diagnostics at a fraction of current costs.


## Overview

- The notebook `TBAS_pipeline_slurm.ipynb` walks through running each step of the TBAS pipeline on a SLURM cluster.
- The `analysis/` folder contains downstream analysis notebooks organized by topic (SNV/Trio calling, coverage, methylation, tandem repeats, and variant counting).
- Benchmarking materials and example HG002 results are provided under `benchmarking/`, with data available at the Zenodo record below.

## Requirements

- Access to a SLURM cluster (e.g., via `sbatch`, `squeue`, etc.).
- A Python environment with Jupyter (JupyterLab or classic notebook).
- Standard long-read analysis command-line tools as used in the notebooks (e.g., samtools, bcftools, bedtools, mosdepth). Consult the individual notebooks for any tool-specific steps and versions used in your environment.

## Getting started

1. Open this repository in JupyterLab or the classic Jupyter Notebook interface.
2. Launch and step through `TBAS_pipeline_slurm.ipynb` to run the pipeline on a SLURM cluster.
	 - Configure any paths to your input data and reference files as indicated in the notebook cells.
	 - Submit jobs or run commands as instructed by the notebook; monitor them with your normal SLURM tooling.
3. Use the analysis notebooks described below to perform downstream analyses and summarize results.

## Run on a SLURM cluster

- Primary entry point: `TBAS_pipeline_slurm.ipynb`.
- Assumptions: you have SLURM available on your compute environment, and the necessary tools are installed and visible in your PATH or modules.
- Tip: If your site uses environment modules, load the appropriate modules at the top of the notebook (or in your kernel startup) before executing the cells.

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

