# TBAS troubleshooting

Real failure modes seen on real data, with fixes. Match the symptom, apply the
fix, and re-run only the affected stages.

## `tdb` fails with "input ... does not exist" (empty tandem-repeat VCF)

Almost always caused by an **incompatible medaka version**. `medaka 2.1.x`
crashes under NumPy ≥ 2.0 with `` `np.compat` was removed in the NumPy 2.0
release ``, but **exits 0** — so every tandem-repeat locus silently fails,
`medaka_to_ref.TR.vcf` is never written, and the failure only surfaces two
stages later at `tdb`.

Fix: use `medaka >= 2.2` (the bundled `environment.yml` pins this). Check with
`medaka --version`. Do **not** downgrade NumPy — moving medaka forward is the
right fix. Then re-run `medaka_local` and `tdb`.

## `medaka_patho` errors with "failed to generate a consensus sequence"

The pathogenic-repeat catalog (`adotto_pheno.bed`) is genome-wide, but a
targeted panel only has coverage over its own regions. If none of the pathogenic
loci fall inside the target regions, medaka has nothing to call.

The pipeline already handles this: it intersects the pathogenic catalog with the
sample's `bed_file` and **skips the stage cleanly** when nothing overlaps,
printing `[medaka_patho] No pathogenic TR loci ... skipping stage`. That is
expected for panels like CMRG that contain no disease-repeat loci — it is not an
error. If you genuinely expect pathogenic loci (e.g. a neurodegeneration panel),
check that `bed_file` actually covers them.

## A tool is "not found" / run stops immediately

The preflight check (`--check-deps`, and the automatic check before a real run)
found a tool missing from PATH. Activate the environment
(`micromamba activate tbas`) or fix the install. To bypass intentionally, pass
`--skip-dep-check`, but that only defers the failure.

## Clair3 stage fails

- The `--clair3-model` path must point at a real model directory matching the
  sequencing chemistry (default assumes R10.4.1 SUP v5.2.0:
  `r1041_e82_400bps_sup_v520`). In the bundled env it is under
  `$CONDA_PREFIX/bin/models/`.
- Clair3 needs the reference `.fai`; create it with `samtools faidx <ref.fa>`.

## minimap2 is slow / uses lots of memory

Expected: it loads the whole genome index (~11 GB RAM, a few minutes) once per
family member. Not a bug. If this is a blocker, that is an optimization item
(prebuilt `.mmi` index) noted in `OPTIMIZATION_PLAN.md`, not a runtime fix.

## Trio results look swapped (proband vs parents)

`sample_id` encodes the barcode order and the **first** barcode is the proband,
then mother, then father (e.g. `4_6_...` → proband=04, mother=05, father=06).
If comparisons look inverted, confirm this ordering with the user and re-check
the manifest.

## The manifest won't load

Use `.csv` (comma) or `.tsv` (tab). The header must include `sample_id`, and
every row needs a non-empty `sample_id`. `bed_file` and `proband_gender` are
required for most stages.
