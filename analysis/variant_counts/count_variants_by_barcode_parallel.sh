#!/usr/bin/env bash
# Minimal checks, parallel per-file
# Usage:
#   NPROC=12 ./count_variants_by_barcode_parallel.sh /path/to/10_12_Gregor_Trio > variant_counts.tsv
# Defaults:
#   NPROC (concurrent processes) defaults to 8

set -euo pipefail

ROOT="${1:-}"
[[ -z "$ROOT" ]] && { echo "Usage: $0 <directory>" >&2; exit 1; }

NPROC="${NPROC:-8}"

# temp files
tmpdir="$(mktemp -d)"
snp_list="${tmpdir}/snp_files.list0"
sv_list="${tmpdir}/sv_files.list0"
snp_out="${tmpdir}/snp_counts.tsv"
sv_out="${tmpdir}/sv_counts.tsv"

# collect file lists (null-delimited)
find "$ROOT" -type f -name '*_barcode*.local.phased.vcf.gz' -print0 > "$snp_list"
find "$ROOT" -type f -name '*_barcode*_global.germline.kanpig.genotyped.vcf.gz' -print0 > "$sv_list"

# SNP counting worker: emits "barcode<TAB>path<TAB>count"
xargs -0 -I{} -P "$NPROC" bash -c '
  f="$1"
  bn="$(basename "$f")"
  if [[ "$bn" =~ barcode([0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
    c=$(bcftools view -H -v snps "$f" | wc -l)
    printf "%s\t%s\t%s\n" "$bc" "$f" "$c"
  fi
' _ < "$snp_list" > "$snp_out"

# SV counting worker: emits "barcode<TAB>path<TAB>count"
xargs -0 -I{} -P "$NPROC" bash -c '
  f="$1"
  bn="$(basename "$f")"
  if [[ "$bn" =~ barcode([0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
    c=$(bcftools view -H "$f" | awk -F"\t" '\''index($8,"SVTYPE=")>0{c++} END{print c+0}'\'')
    printf "%s\t%s\t%s\n" "$bc" "$f" "$c"
  fi
' _ < "$sv_list" > "$sv_out"

# Merge by barcode (full outer join), numeric sort by barcode
# Output header: barcode  snp_vcf_path  snp_count  sv_vcf_path  sv_count
{
  printf "barcode\tsnp_vcf_path\tsnp_count\tsv_vcf_path\tsv_count\n"
  awk -F'\t' '
    FNR==NR { snp_path[$1]=$2; snp_cnt[$1]=$3; seen[$1]=1; next }
             {  sv_path[$1]=$2;  sv_cnt[$1]=$3; seen[$1]=1      }
    END {
      for (b in seen) {
        sp = (b in snp_path) ? snp_path[b] : ""
        sc = (b in snp_cnt ) ? snp_cnt[b]  : 0
        vp = (b in sv_path ) ? sv_path[b]  : ""
        vc = (b in sv_cnt  ) ? sv_cnt[b]   : 0
        printf "%s\t%s\t%s\t%s\t%s\n", b, sp, sc, vp, vc
      }
    }
  ' "$snp_out" "$sv_out" | sort -k1,1n
} 

# cleanup
rm -rf "$tmpdir"

