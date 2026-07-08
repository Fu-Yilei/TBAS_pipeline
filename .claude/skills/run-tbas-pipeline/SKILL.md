---
name: run-tbas-pipeline
description: >-
  Guide a user end-to-end through running the TBAS (Trio-Barcoded ONT Adaptive
  Sampling) long-read pipeline in this repository: checking the software
  environment, building the sample manifest, choosing which stages to run,
  doing a dry run, executing, and interpreting the outputs (small variants,
  structural variants, tandem repeats including pathogenic repeat expansions,
  and methylation). Use this skill whenever someone wants to run TBAS, analyze
  an ONT/nanopore adaptive-sampling trio, process a `calls*.bam`, call
  variants/SVs/tandem-repeats/methylation from long reads with this tool, or
  asks "how do I run this pipeline" while in the TBAS_pipeline repo — even if
  they don't name the stages explicitly. Assume the person may be a clinician
  rather than a bioinformatician and may need plain-language guidance.
---

# Running the TBAS pipeline

The reader of this skill is you, Claude. The person you are helping may be a
clinician or lab scientist, not a bioinformatician, so translate everything into
plain language, never assume they know the jargon, and **confirm before any step
that writes files, is slow, or touches patient data**. Your job is to do the
mechanical work (paths, flags, manifest formatting) so they only make the
decisions that genuinely require a human.

TBAS sequences a rare-disease **trio** (proband + mother + father), barcoded
onto one flow cell with adaptive sampling. This pipeline takes the basecalled,
still-unaligned `calls*.bam` and produces, per trio member: aligned reads, small
variants, structural variants, tandem-repeat genotypes (including pathogenic
expansion loci), and methylation. It is a thin wrapper that runs standard
command-line tools in order — nothing here calls external services or moves
data off the machine.

Work through the steps below in order. Do not skip the dry run.

## Step 0 — Confirm you are set up

1. Confirm you are in the repository: there should be a `pyproject.toml` naming
   `tbas-pipeline` and a `tbas_pipeline/` package. If not, ask where they cloned
   the software and `cd` there.
2. Confirm the toolchain is installed. The repo ships a single self-contained
   conda environment. If `tbas-pipeline` is not already on PATH, guide them:
   ```bash
   micromamba env create -f environment.yml   # or mamba/conda
   micromamba activate tbas
   pip install -e .
   ```
   For a byte-identical environment, create from `environment.lock.yml` instead.
3. **Always** verify the external tools before doing anything else:
   ```bash
   tbas-pipeline --manifest <manifest> --check-deps
   ```
   This lists the tools the selected stages need and flags any missing from
   PATH. If something is missing, fix the environment before continuing — a real
   run can take a long time, and you do not want it to die halfway. If you do not
   have a manifest yet, you can point `--check-deps` at any manifest, including an
   example one under `example_data/`.

## Step 1 — Gather the inputs (conversationally)

Collect these from the user. Ask for what is missing; do not invent paths. If a
file they name does not exist, say so and ask again rather than guessing.

Required:
- **The `calls*.bam`** — the basecalled, demultiplexed reads for the trio (one
  read group per barcode). Often large.
- **A reference genome FASTA** (GRCh38) with its `.fai` index alongside. If the
  `.fai` is missing you can create it with `samtools faidx <ref.fa>`.
- **A Clair3 model directory** — for small-variant calling. In the bundled env
  it is at `$CONDA_PREFIX/bin/models/r1041_e82_400bps_sup_v520` (this must match
  the sequencing chemistry; the default assumes R10.4.1 SUP v5.2.0).
- **`sample_id`** — encodes the trio's barcode range, e.g. `4_6_Gregor_Trio`
  means barcodes 04/05/06, where the **first barcode is the proband**, then
  mother, then father. Confirm that ordering with the user; it matters for the
  trio and pathogenic-expansion comparisons.
- **`proband_gender`** (`male`/`female`) — needed by the tandem-repeat stages
  for sex-chromosome handling.
- **The target regions** (`bed_file`): either a built-in panel name
  (`Epilepsy`, `CMRG`, `Microcephaly`, `Neurodegeneration`, `Arthrogryposis`),
  `WGS` for whole genome, or a direct path to a BED. This is what was enriched.

Optional but common:
- **`adotto_pheno.bed`** — the pathogenic tandem-repeat catalog, needed only if
  they want the `medaka_patho` stage (disease repeat loci such as ATN1, FGF14).
- **`tr_bed_file`** — a tandem-repeat catalog BED for `medaka_local`; if omitted
  it is derived from the panel name.

## Step 2 — Build the manifest

Write a CSV (one row per trio) with a header. Minimum columns:
`sample_id,bed_file,proband_gender`. Add `calls_bam`, `tr_bed_file`, and
`read_group_prefix` when available. Example:

```csv
sample_id,bed_file,tr_bed_file,proband_gender,calls_bam
4_6_Gregor_Trio,CMRG,,female,/data/run1/calls_2025-05-20.bam
```

Show the finished manifest to the user and have them confirm the paths and the
proband/parent barcode order before running anything. Two ready-made real-data
examples live under `example_data/` (`test_subset_chr22` and
`patho_ATN1_FGF14`) — use them to demonstrate a run if the user just wants to
see it work.

## Step 3 — Choose which stages to run

The pipeline has 17 stages that build on each other (each consumes the previous
one's outputs), grouped into analyses. Ask the user what they want using
checkbox-style options via the AskUserQuestion tool. Offer these choices and let
them multi-select; recommend "Everything" as the default:

- **Everything** (recommended) — the full 17-stage run.
- **Alignment & coverage** — demultiplex → align → sort, plus SV calling
  (sniffles) and coverage (mosdepth).
- **Small variants** — Clair3 + kanpig genotyping (incl. trio-aware).
- **Tandem repeats** — medaka local catalog + pathogenic loci + tdb database.
- **Methylation** — modkit phased and unphased pileups.

Important: because stages depend on earlier outputs, running a later group
requires the earlier ones to have completed already (their files must exist in
the output folder). If the user picks only a downstream group on a fresh run,
warn them it will fail without the upstream outputs, and offer to run the full
pipeline or the needed prerequisites first. Translate their choice into the
concrete `--stages a,b,c` list (the stage names are in `references/stages.md`).

## Step 4 — Dry run first (do not skip)

Always show the exact commands before executing:

```bash
tbas-pipeline \
  --manifest <manifest.csv> \
  --output-folder <output_dir> \
  --reference <ref.fa> \
  --clair3-model <clair3_model_dir> \
  --adotto-pheno <adotto_pheno.bed> \
  --stages <selected or omit for all> \
  --dry-run
```

Read the printed commands back to the user in plain language ("this step splits
the reads by family member", "this one calls small variants", …). This is the
moment to catch a wrong path or panel before spending compute. Get explicit
confirmation to proceed.

## Step 5 — Execute

Run the same command without `--dry-run`. This is long-running: alignment loads
the whole genome index per family member, and the tandem-repeat/methylation
stages are compute-heavy. Run it in the background and monitor progress rather
than blocking. If a real run is missing a tool, the pipeline stops early with a
clear message (that is the preflight check doing its job) — fix and rerun.

## Step 6 — Interpret the outputs

Outputs land under `<output_folder>/<sample_id>/`, per barcode. Point the user to
what they care about:
- Small variants: `*_clair3_local/merge_output.vcf.gz`, and kanpig-genotyped /
  trio VCFs (`*.kanpig.genotyped.vcf.gz`, `*.kanpig.trio.genotyped.vcf.gz`).
- Structural variants: `*_global.germline.vcf.gz` (sniffles).
- Tandem repeats: `*_tr/medaka_to_ref.TR.vcf` and the `*.tdb` database;
  pathogenic loci in `*_pheno_tr/`. For a trio, compare the proband's repeat
  length against both parents to spot an expansion.
- Methylation: `*_modkit/*.bed` (phased `.HP.bed` and unphased).
- Coverage: `*_mosdepth/`.

If the user is chasing a specific diagnosis, tie it back: e.g. an expanded
`CAG`/`GAA` at a known disease gene in the proband but not the parents suggests a
de novo or inherited pathogenic expansion — but frame findings as candidates for
their own clinical review, not conclusions.

## When things go wrong

Read `references/troubleshooting.md` for the common, real failure modes and their
fixes (medaka/NumPy version crash, `medaka_patho` skipping when a panel has no
pathogenic loci, missing Clair3 model or `.fai`, wrong barcode order). Consult
`references/stages.md` for the full ordered stage list, the exact stage names to
pass to `--stages`, and what each stage consumes and produces.

## Guardrails

- This is clinical / patient data. Do not copy inputs or outputs anywhere the
  user did not ask for, and do not send anything to external services.
- Never fabricate file paths, sample IDs, or results. If unsure, ask.
- Confirm before long runs and before overwriting an existing output folder.
- Report failures honestly with the actual error; do not claim success you did
  not verify.
