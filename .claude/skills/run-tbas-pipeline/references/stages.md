# TBAS pipeline stages

The 17 stages run in this fixed order; each consumes the previous stages'
outputs. Pass a subset with `--stages name1,name2,...` (comma-separated, no
spaces), or omit `--stages` to run all. The stage **names** below are exactly
what `--stages` expects.

| # | Stage name | Plain meaning | Key output |
|---|------------|---------------|------------|
| 1 | `demultiplex` | Split the reads by family member (barcode) | per-barcode `.bam` |
| 2 | `fastq_extract` | Pull reads to FASTQ, keeping methylation tags | `.fastq` |
| 3 | `minimap2` | Align reads to the genome | `.sam` |
| 4 | `bam_sort` | Sort + index the alignments | `.sorted.bam` |
| 5 | `sniffles_global` | Call structural variants | `*_global.germline.vcf.gz` |
| 6 | `bam_mosdepth` | Coverage over the target regions | `*_mosdepth/` |
| 7 | `clair3_local` | Call small variants (SNVs/indels) | `*_clair3_local/merge_output.vcf.gz` |
| 8 | `kanpig_pileup` | Prepare SV genotyping input | `*_kanpig.plup.gz` |
| 9 | `kanpig_gt` | Genotype SVs per sample | `*.kanpig.genotyped.vcf.gz` |
| 10 | `kanpig_trio` | Trio-aware SV genotyping | `*.kanpig.trio.genotyped.vcf.gz` |
| 11 | `whatshap_single_sample_local_phasing` | Phase small variants | `*.local.phased.vcf.gz` |
| 12 | `whatshap_haplotag` | Tag reads by haplotype | `*.HP.bam` |
| 13 | `medaka_local` | Genotype tandem repeats (catalog) | `*_tr/medaka_to_ref.TR.vcf` |
| 14 | `medaka_patho` | Genotype pathogenic repeat loci | `*_pheno_tr/medaka_to_ref.TR.vcf` |
| 15 | `modkit` | Methylation pileup, phased | `*_modkit/*.HP.bed` |
| 16 | `modkit_nohp` | Methylation pileup, unphased | `*_modkit/*.bed` |
| 17 | `tdb` | Build tandem-repeat database | `*_tr/*.tdb` |
| 18 | `merge_trio_snv` | Merge the trio's small variants (proband first) | `*.snv.merged.vcf.gz` |
| 19 | `annovar` | Annotate small variants (ANNOVAR) | `*.snv.merged.hg38_multianno.txt` |
| 20 | `rank_snv` | Tier & rank SNVs (Tier 1–4) | `*.snv.ranked.tsv` |
| 21 | `annotsv` | Annotate structural variants (AnnotSV) | `*.annotsv.tsv` |
| 22 | `rank_sv` | Keep (likely) pathogenic SVs | `*.sv.pathogenic.tsv` |
| 23 | `rank_tr` | Trio tandem-repeat comparison (proband vs parents) | `*.tr.trio_distinct.tsv` |

Stages 18–23 are the annotation/ranking layer (see `METHODS.md`). They need
extra tools/resources: ANNOVAR (`table_annovar.pl` + a humandb directory,
`--annovar-humandb`), AnnotSV, and a STRchive BED (`--strchive-bed`). If a site
does not have these, run `--stages` up to `tdb` and do ranking separately.

## Dependency notes

- Stages 8–17 all need the sorted alignments (stage 4).
- Phasing/haplotag (11–12) need Clair3 output (7).
- The tandem-repeat and methylation stages (13–16) need the haplotagged BAM (12).
- `tdb` (17) needs `medaka_local` (13).
- So on a **fresh** output folder, do not start midway. Run the full pipeline,
  or start from `demultiplex`. You can safely *resume* from a later stage only
  when the earlier outputs already exist in the output folder.

## Analysis groupings (for the checkbox question)

- Alignment & coverage → `demultiplex,fastq_extract,minimap2,bam_sort,sniffles_global,bam_mosdepth`
- Small variants → `clair3_local,kanpig_pileup,kanpig_gt,kanpig_trio,whatshap_single_sample_local_phasing,whatshap_haplotag`
- Tandem repeats → `medaka_local,medaka_patho,tdb`
- Methylation → `modkit,modkit_nohp`
