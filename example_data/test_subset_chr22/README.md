# TBAS Example Test Data (chr22 subset)

This bundle contains the exact subset BAM used in pipeline testing.

- Sample: `4_6_Gregor_Trio`
- Region: `chr22:42600000-42740000`

Files:
- `example_data/test_subset_chr22/4_6_Gregor_Trio/calls_chr22_42600000_42740000.bam`
- `example_data/test_subset_chr22/4_6_Gregor_Trio/calls_chr22_42600000_42740000.bam.bai`
- `example_data/test_subset_chr22/input_regions_cmrg_chr22_42600000_42740000.bed`
- `example_data/test_subset_chr22/tr_regions_adotto_cmrg_chr22_42600000_42740000.bed`
- `example_data/test_subset_chr22/manifest_example.csv`

Run example:

```bash
tbas-pipeline \
  --manifest example_data/test_subset_chr22/manifest_example.csv \
  --output-folder /path/to/output_root
```
