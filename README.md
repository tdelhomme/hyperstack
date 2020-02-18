# needlestack-ITC

## Description
Nextflow pipeline to run the intercept method on needlestack calling with a step of FDR computation with ML

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io).

2. External software:
- Rscript
- bcftools
- bgzip/tabix


## Input
  | Type      | Description     |
  |-----------|---------------|
  | train_vcf    | Input VCF used to train the random forest algorithm (bgzip/tabix format). |
  | apply_vcf    | Input VCF on which the random forest algorithm will be applied (bgzip/tabix format). |


## Parameters

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --output_folder   |      needlestack-ITC_result | Output folder  |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run tdelhomme/needlestack-ITC --train_vcf gold.vcf.bgz --apply_vcf test.vcf.bgz
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | ITC | Folder containing a .Rdata object per sample corresponding to the list of estimated ITC values |
  | PLOTS | Folder containing the output plots  |
  | RDATA | Folder containing the trained random forest in .Rdata format |
  | VCF | Folder containing annotated train and apply vcfs with a geno statistic: FPRF (false positive rate from the random forest) |
