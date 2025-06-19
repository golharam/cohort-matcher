#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  job_uuid: string
  bam_files:
    type: File[]
    secondaryFiles: [^.bai]
  reference:
    type: File
    secondaryFiles: [.fai, ^.dict]
  bed_file: File
  thread_count: int
  number_of_chunks: int
  output_prefix: string
  cromwell_engine: boolean

outputs:
  time_metrics_from_freebayes:
    type: File
    outputSource: genomel_pdc_freebayes/time_metrics
  time_metrics_from_picard_sortvcf:
    type: File
    outputSource: picard_sortvcf/time_metrics
  time_metrics_from_selectvariants:
    type: File
    outputSource: gatk3_selectvariants/time_metrics
  freebayes_vcf:
    type: File
    outputSource: gatk3_selectvariants/output_vcf

steps:
  genomel_pdc_freebayes:
    run: ../../tools/variant_calling/genomel_pdc_freebayes.cwl
    in:
      job_uuid: job_uuid
      bam_files: bam_files
      reference: reference
      bed_file: bed_file
      thread_count: thread_count
      number_of_chunks: number_of_chunks
      cromwell_engine: cromwell_engine
    out: [vcf_list, time_metrics]
