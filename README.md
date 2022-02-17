# hyperstack

<img align="center" src="https://github.com/tdelhomme/needlestack-ITC/blob/master/method.png" width="600">

## Description
Nextflow pipeline to run the hyperstack method on needlestack calling with a step of FDR computation with ML.

This pipeline first trains a random forest algorithm based a set of true and false positve variant, defined with one of the following method:
  * on ethnicity of the samples in the VCF: TP are variants reported in the same ethnicity and FP are variants reported in another ethnicity, based on 1000 genome project database
  * on recurrence of mutations: variant is TP if observed >= x times in other samples (e.g. normal single cells) and FP is observed only once

Then it applies this random forest on a second VCF (new genotype statistic: FPRF).

Finally it computes the proportion of each base subtitution (96-profile) estimated form our intercept method.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io).

2. External software:
- Rscript
- bcftools
- bgzip/tabix


## Input
  | Type      | Description     |
  |-----------|---------------|
  | train_vcf    | Input VCF used to train the random forest algorithm (bgzip/tabix format and **annotated with annovar**). |
  | apply_vcf    | Input VCF on which the random forest algorithm will be applied (bgzip/tabix format). |


## Parameters

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --output_folder   |      needlestack-ITC_result | Output folder  |
| --germline_file   |      txt file containing one germline per line to be removed from the analysis (format: chr1:1245716_C/T)| TXT file  |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run tdelhomme/hyperstack --train_vcf gold.vcf.bgz --apply_vcf test.vcf.bgz
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | ITC | Folder containing a .Rdata object per sample corresponding to the list of estimated ITC values |
  | PLOTS | Folder containing the output plots  |
  | RDATA | Folder containing the trained random forest in .Rdata format |
  | VCF | Folder containing annotated train and apply vcfs with a geno statistic: FPRF (false positive rate from the random forest) |
