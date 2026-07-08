# TBAS Example Test Data — pathogenic tandem repeats (ATN1 + FGF14)

This bundle exercises the **full** pipeline end-to-end, including the
`medaka_patho` and `tdb` stages, by targeting two clinically relevant
pathogenic tandem-repeat loci that have real coverage in this trio:

- **ATN1** (DRPLA), `chr12:6,936,716-6,936,773`, `CAG` repeat
- **FGF14** (SCA27B), `chr13:102,161,577-102,161,726`, `GAA` repeat

Unlike the `test_subset_chr22` bundle (CMRG genes, which contain no known
pathogenic STR loci and therefore make `medaka_patho` skip), both loci here are
present in the pathogenic catalog `adotto_pheno.bed`, so `medaka_patho`
genotypes them and `tdb` builds a database from the results.

- Sample: `4_6_Gregor_Trio` (barcodes 04/05/06)
- Regions: `chr12:6,900,000-6,980,000` and `chr13:102,120,000-102,200,000`
  (~80 kb around each locus, enough flanking sequence for local phasing)

Files:
- `4_6_Gregor_Trio/calls_patho_ATN1_FGF14.bam` — unaligned dorado reads
  (subset of the full trio BAM by read name), one read group per barcode
- `input_regions_patho_ATN1_FGF14.bed` — target regions (`bed_file`)
- `tr_regions_adotto_patho_ATN1_FGF14.bed` — Adotto TR catalog loci in-window
  for `medaka_local` (includes the ATN1 `CAG` and FGF14 `GAA` loci)
- `manifest_example.csv`

Run example (full pipeline):

```bash
tbas-pipeline \
  --manifest example_data/patho_ATN1_FGF14/manifest_example.csv \
  --output-folder /path/to/output_root \
  --reference /path/to/GRCh38.fa \
  --clair3-model /path/to/clair3_model_dir \
  --adotto-pheno /path/to/adotto_pheno.bed
```
