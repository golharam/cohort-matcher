import csv

def read_samplesheet(sample_sheet):
    '''
    Read the samplesheet consisting of a tab-separated file 
    with a header with the following columns:
    'subject_id', 'sample_id', 'bamfile', 'reference', 'assay'
    '''
    with open(sample_sheet, 'r') as file_handle:
        reader = csv.DictReader(
            file_handle, delimiter='\t'
        )
        samples = [row for row in reader]
    return samples
