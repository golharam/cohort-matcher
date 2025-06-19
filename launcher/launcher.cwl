cwlVersion: v1.2
class: CommandLineTool
label: Cohort-Matcher Launcher

requirements:
  WorkReuse:
    enableReuse: false
  DockerRequirement:
    dockerPull: ngsbioinformatics/examplelauncher:latest

baseCommand: ["/opt/launcher/launcher.py"]

inputs:
  sample_sheet:
    doc: Input sample sheet with the columns (subject id, sample id, bam file path, reference)
    type: File
    inputBinding:
      prefix: --sample_sheet

outputs:
  log:
    type: File
    outputBinding:
      glob: '*.log'
