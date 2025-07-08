cwlVersion: v1.2
class: Workflow
label: Cohort-matcher Workflow

inputs:
  bamsheet:
    label: BAM Sheet (sample id, bam flie, and reference genome)
    type: File
  targets_bed_file:
    label: Targets BED file
    type: File
  reference_fasta:
    label: Reference FASTA file
    type: File

steps:
    parse_samplesheet:
        in:
            bamsheet: bamsheet
        out: [sample_list]
        run: tools/parse_samplesheet.cwl

    genotype_samples:
        in:
          sample_name: parse_samplesheet/sample_list.sample_id
          bamfile: parse_samplesheet/sample_list.bamfile
          targets_bed_file: inputs.targets_bed_file
          reference_fasta: inputs.reference_fasta
        run: tools/freebayes.cwl
        out: [vcf_file]
        scatter: [sample_name, bamfile]
        scatterMethod: dotproduct

    compare_genotypes:
        
    merge_results:

    plot_snps:

    report_matches:

outputs:
