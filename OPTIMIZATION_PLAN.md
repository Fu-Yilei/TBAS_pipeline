# TBAS Pipeline — Codebase Analysis & Optimization Plan

This document explains what the pipeline does, assesses it for clinical /
push-button use, and lays out a phased plan to make it easy to install and run
without deep bioinformatics expertise.

## 1. What the software does

TBAS (Trio-Barcoded ONT Adaptive Sampling) analyzes Oxford Nanopore long-read
sequencing of rare-disease **trios** (proband + mother + father) that were
barcoded onto a single flow cell with adaptive sampling. The Python package
(`tbas_pipeline/`) is a **thin orchestration layer**: it reads a manifest,
builds shell command strings, and runs external bioinformatics tools in order.
It carries no scientific Python dependencies of its own (`dependencies = []`);
all real work is done by command-line tools.

The input is a single basecalled/demultiplexed `calls*.bam` containing one
read-group per barcode. `sample_id` (e.g. `4_6_Gregor_Trio`) encodes the trio's
barcode range (`4_6` → barcode04/05/06). The 17 stages run per barcode:

| # | Stage | Tool | Purpose |
|---|-------|------|---------|
| 1 | demultiplex | samtools view | split BAM by read-group into per-barcode BAMs |
| 2 | fastq_extract | samtools fastq | extract reads + methylation tags |
| 3 | minimap2 | minimap2 | align to GRCh38 |
| 4 | bam_sort | samtools sort | sort + index |
| 5 | sniffles_global | sniffles | structural-variant calling |
| 6 | bam_mosdepth | mosdepth | coverage over target BED |
| 7 | clair3_local | run_clair3.sh | small-variant (SNV/indel) calling |
| 8–10 | kanpig_pileup/gt/trio | kanpig | SV genotyping incl. trio-aware |
| 11–12 | whatshap phase/haplotag | whatshap | phasing + haplotagging reads |
| 13–14 | medaka_local/patho | medaka tandem | tandem-repeat consensus (catalog + pathogenic loci) |
| 15–16 | modkit/modkit_nohp | modkit | methylation pileups (phased + unphased) |
| 17 | tdb | tdb | build tandem-repeat database |

The `analysis/` and `benchmarking/` notebooks are downstream/manual and are not
part of the automated pipeline.

## 2. Why it is hard to use clinically today

The orchestration code is clean, but the *operational surface* is not
push-button. Concrete blockers:

1. **Dependency sprawl with no environment spec.** The run needs 13 external
   tools. In the current dev machine they come from four different places:
   a base micromamba env, a separate `clair3` conda env, a hand-built Rust
   binary (`kanpig`), and a loose binary in a home directory (`modkit`).
   Nothing declares this, so a new user must reconstruct it by trial and error.
2. **GitHub / from-source installs.** `clair3`, `kanpig`, `modkit`, `tdb` are
   commonly installed from GitHub or built from source. That is exactly the
   friction that blocks a clinical lab.
3. **Hardcoded, machine-specific default paths.** Defaults such as
   `GRCh38-2.1.0_genome_mainchrs.fa`, `clair3/bin/models/...`, and the
   `BED_FILE_MAP` entries (`adaptive_regions/annotated/merged/...`) assume a
   specific working directory. They are not shipped with the package and not
   validated, so a fresh clone cannot run the defaults.
4. **No fail-fast validation (partly addressed in this change).** Until now, a
   missing tool or reference surfaced only when a stage crashed — potentially
   after hours. (This PR adds a preflight tool check; input-file validation is
   still open.)
5. **Not resumable / not idempotent.** Every invocation re-runs every selected
   stage from scratch. A crash in stage 14 means re-running 1–13. There is no
   "skip if output exists" and no checkpointing.
6. **No parallelism, no per-sample isolation.** Stages loop over barcodes
   serially, and any single failure aborts the whole run (`check=True` with
   `shell=True`). One bad sample takes down the batch.
7. **No provenance / logging.** Commands print to stdout but nothing is written
   to a per-run log, no tool versions are captured, no structured summary is
   produced — all of which a clinical/accredited setting expects.
8. **Reproducibility gaps / silent version breakage.** Model strings are
   hardcoded (`clair3 ...v520`, `medaka ...@v5.2.0`) and nothing pins tool
   versions to match them. This is not hypothetical: on a real end-to-end run
   of the shipped `example_data`, `medaka 2.1.1` (against NumPy 2.2) crashed on
   every tandem-repeat locus with `` `np.compat` was removed in the NumPy 2.0
   release `` yet **still exited 0**, producing an empty TR VCF. The failure
   only surfaced two stages later when `tdb` aborted with "input does not
   exist." `medaka >=2.2` fixes the NumPy-2 incompatibility. Two lessons: (a)
   pin tool versions (now noted in `environment.yml`), and (b) stages that can
   exit 0 on partial failure (medaka tandem) need explicit output validation.
9. **Efficiency.** `minimap2` is given the raw FASTA, so it rebuilds the whole
   genome index (~3 min, ~11 GB RAM) once *per barcode*. A prebuilt `.mmi`
   index would remove most of that cost.

## 3. Optimization plan (phased)

### Phase 0 — Quick wins (this PR, low risk)
- [x] Bump minimum Python to 3.10; add metadata/classifiers.
- [x] Preflight dependency check (`--check-deps`, auto fail-fast before real
      runs, `--skip-dep-check` escape hatch).
- [x] `environment.yml` so the whole toolchain installs into one self-contained
      env with a single `micromamba env create`.
- [x] `medaka_patho` no longer crashes the run when the target panel contains no
      pathogenic STR loci: the catalog is restricted to the target regions and
      the stage skips cleanly when nothing overlaps.
- [ ] Preflight **input** validation: reference FASTA (+ `.fai`), clair3 model
      dir, and BED files exist before running; build the `.fai` if missing.

### Phase 1 — One-command install (packaging)
- [x] Verified the entire toolchain installs into a single self-contained
      `tbas` env (no separate clair3/medaka envs, no from-source builds). Two
      package-name gotchas resolved: modkit ships as `ont-modkit`, and `tdb` is
      not on bioconda/PyPI so it installs from GitHub via pip.
- [x] Committed a pinned `environment.lock.yml` (from `micromamba env export`)
      alongside the loose `environment.yml` for reproducible installs.
- Ship a **container** (Docker + an Apptainer/Singularity build for HPC/clinical
  environments) that bakes in the toolchain and models. A container is the
  strongest reproducibility guarantee and is the format clinical pipelines
  usually require.
- Optionally submit the Python package to **bioconda / PyPI** so
  `pip install tbas-pipeline` or `mamba install tbas-pipeline` works.

### Phase 2 — Configuration over hardcoding
- Replace the hardcoded `PipelineSettings` defaults and `BED_FILE_MAP` /
  `TR_BED_FILE_MAP` with a **YAML/TOML config file** (`--config tbas.yaml`)
  covering reference, models, kit, thread counts, and region catalogs.
- **Bundle the region BEDs as package data** (or fetch them by version) so the
  panel names (`Epilepsy`, `CMRG`, …) resolve without a specific CWD.
- Make thread counts configurable instead of per-stage constants.

### Phase 3 — Robustness & resumability
- Add **idempotency**: skip a stage when its expected outputs already exist
  (with a `--force` override). This alone makes long runs recoverable.
- Write a **per-run log file** and a **run manifest** capturing tool versions,
  command lines, and exit status for provenance.
- Isolate per-sample failures (continue-on-error with a summary) instead of
  aborting the whole batch.
- Build the `minimap2` `.mmi` index once and reuse it across barcodes.

### Phase 4 — Workflow engine (for scale/clinical)
- For production/clinical throughput, consider re-expressing the DAG in
  **Snakemake or Nextflow**. These give resumability, parallel scheduling,
  cluster/cloud execution, and provenance essentially for free, and are the
  de-facto standard for accredited genomics pipelines. The current Python layer
  can remain as the reference/spec for the command lines.

### Phase 5 — Validation & CI
- Keep the current fast command-string unit tests, and add a **real-data
  integration test** using the shipped `example_data/test_subset_chr22` bundle
  (chr22:42.6–42.74 Mb, ~3k reads) that runs the full pipeline in CI against a
  containerized toolchain. This is what proves "it doesn't break on real data"
  on every change.

## 4. Suggested near-term priority

1. Phase 0 input validation (finish fail-fast).
2. Phase 1 container + pinned lock file (removes the biggest install barrier).
3. Phase 3 idempotency (makes real runs practical to operate).

These three, on top of what this PR already lands, convert the tool from
"requires the original author's machine" to "clone, `mamba env create`, run."
