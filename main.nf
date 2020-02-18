#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2020 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  needlestack-ITC-nf: nextflow pipeline to run the intercept method "
log.info "  on needlestack variant calling with a step of ML to compute FDRs       "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --train_vcf train.vcf.bgz --apply_vcf apply.vcf.bgz'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --train_vcf                    FILE         input vcf for training.'
    log.info '    --apply_vcf                    FILE         input vcf to apply the model.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                FOLDER         Output folder (default: needlestack-ITC_result).'
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.train_vcf = null
params.apply_vcf = null
params.output_folder = "needlestack-ITC_result"

if(params.train_vcf == null | params.apply_vcf == null ){
  exit 1, "Please specify each of the following parameters: --train_vcf, --apply_vcf and --gold_list"
}

train = Channel.fromPath(params.train_vcf).view { "value: $it" }
train_tbi = file(params.train_vcf + '.tbi')
apply = Channel.fromPath(params.apply_vcf).view { "value: $it" }
apply_tbi = file(params.apply_vcf + '.tbi')

process training {

  publishDir params.output_folder+"/PLOTS/", mode: 'copy', pattern: "*.pdf"
  publishDir params.output_folder+"/RDATA/", mode: 'copy', pattern: "*.Rdata"
  publishDir params.output_folder+"/VCF/", mode: 'copy', pattern: "*.vcf"

  input:
  file t from train
  file train_tbi

  output:
  file "*.pdf" into plots
  file "*.Rdata" into rdata
  file "*.vcf" into vcftrain

  shell:
  '''
  Rscript !{baseDir}/bin/FDR_RF_train.r --vcf=!{t}
  '''
}

process application {

  publishDir params.output_folder+"/VCF/", mode: 'copy', pattern: "*.vcf"

  input:
  file a from apply
  file apply_tbi
  file model from rdata

  output:
  file "*.vcf" into vcfapply

  shell:
  '''
  Rscript !{baseDir}/bin/FDR_RF_apply.r --vcf=!{a} --model=!{model}
  '''
}

process extract_samples {

  input:
  file v from vcfapply

  output:
  file "*vcf" into sample_vcf mode flatten

  shell:
  '''
  ln -s `ls *RF_needlestack.vcf` input.vcf
  for sample in `bcftools query -l input.vcf`
  do
    bcftools view -s ${sample} input.vcf > ${sample}.vcf
  done

  '''
}

process ITC {

  tag {sm}

  publishDir params.output_folder+"/ITC/", mode: 'copy', pattern: "*.Rdata"
  publishDir params.output_folder+"/PLOT/", mode: 'copy', pattern: "*.pdf"

  input:
  file v from sample_vcf

  output:
  file "*Rdata" into itc
  file "*pdf" into outpdf

  shell:
  sm = v.baseName
  '''
  bgzip -c !{v} > !{v}.bgz
  tabix -p vcf !{v}.bgz
  Rscript !{baseDir}/bin/intercept_method.r --vcf=!{v}.bgz --bin_path=!{baseDir}/bin --sm=!{sm} --fdr_range=0.2,0.8
  '''

}
