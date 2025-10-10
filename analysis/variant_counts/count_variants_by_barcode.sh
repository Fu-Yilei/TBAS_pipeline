#!/usr/bin/env bash

# Usage: ./count_variants_by_samplekey.sh <directory>
ROOT="$1"
[[ -z "$ROOT" ]] && { echo "Usage: $0 <directory>" >&2; exit 1; }

# One line per sample key:
# sample_key = "<dirpath>|<basename_prefix_up_to__barcodeNN>"
# We also keep the barcode separately for readability in the output.

declare -A SNP_PATH SNP_CNT SV_PATH SV_CNT BARCODE_OF_KEY

# helper: extract sample_key and barcode from a path
extract_key_and_barcode() {
  local f="$1"
  local dir base barcode prefix key
  dir="$(dirname "$f")"
  base="$(basename "$f")"
  # barcodeNN
  if [[ "$base" =~ (.*_barcode)([0-9]+) ]]; then
    prefix="${BASH_REMATCH[1]}${BASH_REMATCH[2]}"  # everything through _barcodeNN
    barcode="${BASH_REMATCH[2]}"
    key="${dir}|${prefix}"
    printf "%s\t%s\n" "$key" "$barcode"
  else
    # no barcode pattern; skip by returning empty
    printf "\t\n"
  fi
}

# scan files and count on the fly (no indexing by barcode)
# SNP files
while IFS= read -r -d '' f; do
  kv="$(extract_key_and_barcode "$f")"
  key="${kv%%$'\t'*}"
  bc="${kv##*$'\t'}"
  [[ -z "$key" ]] && continue

  # Store barcode for this key (last write wins, but same barcode expected per key)
  BARCODE_OF_KEY["$key"]="$bc"
  SNP_PATH["$key"]="$f"

  # Count true SNPs
  SNP_CNT["$key"]=$(bcftools view -H -v snps "$f" | wc -l)
done < <(find "$ROOT" -type f -name '*_barcode*.local.phased.vcf.gz' -print0)

# SV files
while IFS= read -r -d '' f; do
  kv="$(extract_key_and_barcode "$f")"
  key="${kv%%$'\t'*}"
  bc="${kv##*$'\t'}"
  [[ -z "$key" ]] && continue

  BARCODE_OF_KEY["$key"]="$bc"
  SV_PATH["$key"]="$f"

  # Count records with SVTYPE= in INFO
  SV_CNT["$key"]=$(bcftools view -H "$f" | awk -F'\t' 'index($8,"SVTYPE=")>0{c++} END{print c+0}')
done < <(find "$ROOT" -type f -name '*_barcode*_global.germline.kanpig.genotyped.vcf.gz' -print0)

# Output
# Columns:
# sample_key  barcode  snp_vcf_path  snp_count  sv_vcf_path  sv_count
printf "sample_key\tbarcode\tsnp_vcf_path\tsnp_count\tsv_vcf_path\tsv_count\n"

# sort by sample_key for stable output
mapfile -t KEYS < <(printf "%s\n" "${!BARCODE_OF_KEY[@]}" | sort)

for key in "${KEYS[@]}"; do
  bc="${BARCODE_OF_KEY[$key]}"
  snp_p="${SNP_PATH[$key]:-}"
  sv_p="${SV_PATH[$key]:-}"
  snp_c="${SNP_CNT[$key]:-0}"
  sv_c="${SV_CNT[$key]:-0}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$key" "$bc" "${snp_p}" "${snp_c}" "${sv_p}" "${sv_c}"
done

