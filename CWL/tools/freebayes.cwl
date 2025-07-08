cwlVersion: v1.2
class: CommandLineTool
label: freebayes

inputs:
    sample_name:
        type: string
    bamfile:
        type: File
    targets_bed_file:
        type: File
    reference_fasta:
        type: File

baseCommand: [freebayes, --no-indels, --min-coverage, 15, --report-all-haplotype-alleles, --report-monomorphic]

arguments:
    - position: 0
      prefix: --fasta-reference
      valueFrom: $(inputs.reference_fasta.path)
    - position: 0
      prefix: --targets
      valueFrom: $(inputs.targets_bed_file.path)    
    - position: 0
      prefix: --vcf
      valueFrom: $(inputs.sample_name).vcf
    - position: 1
      valueFrom: $(inputs.bamfile.path)

outputs:
    vcf_file:
      type: File
      outputBinding:
        glob: $(inputs.sample_name).vcf
