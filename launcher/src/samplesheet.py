import csv

def read_samplesheet(sample_sheet):
    ''' Read the samplesheet '''
    with open(sample_sheet, 'r') as file_handle:
        reader = csv.DictReader(
            file_handle, delimiter='\t',
            fieldnames=['subject_id', 'sample_id', 'bamfile', 'reference']
        )
        samples = [row for row in reader]
    return samples
