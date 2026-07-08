#!/usr/bin/env bash
# Download the ANNOVAR GRCh38 annotation databases used by the `rank_snv` stage.
#
# ANNOVAR itself is NOT redistributable and is not on conda: register and
# download it from
#   https://www.openbioinformatics.org/annovar/annovar_download_form.php
# then put its scripts on PATH (so `annotate_variation.pl` and
# `table_annovar.pl` are callable) before running this.
#
# Usage:
#   scripts/download_annovar_db.sh [HUMANDB_DIR]
#
# Default HUMANDB_DIR is ./annovar/humandb (matches the pipeline default; pass a
# different path here and to `--annovar-humandb`). gnomad41_genome is large
# (tens of GB) — ensure you have the space and time.
set -euo pipefail

HUMANDB="${1:-annovar/humandb}"
BUILD=hg38

if ! command -v annotate_variation.pl >/dev/null 2>&1; then
    echo "error: annotate_variation.pl not on PATH — install ANNOVAR first." >&2
    exit 1
fi

mkdir -p "$HUMANDB"

# The exact databases (and versions) used by the study, all served by ANNOVAR.
# If a pinned version is no longer hosted (e.g. a newer ClinVar), substitute the
# closest available build and update METHODS.md accordingly.
DBS=(
    refGeneWithVer        # RefSeq genes (versioned transcripts)
    avsnp151              # dbSNP 151
    dbnsfp47a             # dbNSFP v4.7a — provides REVEL_score, CADD_phred
    clinvar_20250721      # ClinVar (CLNSIG)
    1000g2015aug          # 1000 Genomes AF (ALL + AFR/AMR/EAS/EUR/SAS)
    gnomad41_genome       # gnomAD genome v4.1
)

for db in "${DBS[@]}"; do
    echo "=== downloading ${db} into ${HUMANDB} ==="
    annotate_variation.pl -buildver "$BUILD" -downdb -webfrom annovar "$db" "$HUMANDB/"
done

echo
echo "Done. Run the pipeline with:  --annovar-humandb ${HUMANDB}"
