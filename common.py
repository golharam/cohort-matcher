import boto3
import logging
import os
import uuid

logger = logging.getLogger(__name__)

def downloadFile(srcFile, destFile):
    ''' Download file from S3 '''
    s3 = boto3.resource('s3')
    bucket, key = find_bucket_key(srcFile)
    s3.meta.client.download_file(bucket, key, destFile)

def uploadFile(srcFile, destFile):
    ''' Upload a file to S3 '''
    s3 = boto3.resource('s3')
    bucket, key = find_bucket_key(destFile)
    s3.meta.client.upload_file(srcFile, bucket, key,
                               ExtraArgs={'ServerSideEncryption': 'AES256'})

def find_bucket_key(s3path):
    """
    This is a helper function that given an s3 path such that the path is of
    the form: bucket/key
    It will return the bucket and the key represented by the s3 path, eg
    if s3path == s3://bmsrd-ngs-data/P-234
    """
    if s3path.startswith('s3://'):
        s3path = s3path[5:]
    s3components = s3path.split('/')
    bucket = s3components[0]
    s3key = ""
    if len(s3components) > 1:
        s3key = '/'.join(s3components[1:])
    return bucket, s3key

def listFiles(s3_cache_dir, suffix=None):
    ''' Return a list of files in s3_cache_dir '''
    files = []
    bucketName, folder = find_bucket_key(s3_cache_dir)
    s3bucket = 's3://%s/' % bucketName
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucketName)
    for obj in bucket.objects.filter(Prefix=folder + '/'):
        if suffix:
            if obj.key.endswith(suffix):
                files.append(s3bucket + obj.key)
        else:
            files.append(s3bucket + obj.key)
    return files

def readSamples(sampleSheetFile):
    '''
    readSamples reads in a sampleSheetFile consisting of two columns:
    name and bamfile
    :param sampleSheetFile: tab-delimited file of samplename and s3 bamfile path
    :return: list of {name, bam} dictionaries
    '''
    if os.path.isfile(sampleSheetFile) is False:
        logger.error("%s does not exist", sampleSheetFile)
        return False
    logger.info("Reading %s", sampleSheetFile)
    sampleNames = []
    samples = []
    with open(sampleSheetFile, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if len(line) == 0 or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) != 2:
                logger.error("Expect 2 fields (sampleName, bamFile) but encountered %d",
                             len(fields))
                return False

            sampleNames.append(fields[0])
            sample = {"name": fields[0],
                      "bam": fields[1]}
            samples.append(sample)
    if len(sampleNames) != len(set(sampleNames)):
        logger.error("Duplicate sampleids found in %s", sampleSheetFile)
        return False
    logger.info("Read %d samples.", len(samples))
    return samples

def generate_working_dir(working_dir_base):
    """
    Creates a unique working directory to combat job multitenancy
    :param working_dir_base: base working directory
    :return: a unique subfolder in working_dir_base with a uuid
    """

    working_dir = os.path.join(working_dir_base, str(uuid.uuid4()))
    try:
        os.mkdir(working_dir)
    except Exception as e:
        return working_dir_base
    return working_dir

def delete_working_dir(working_dir):
    """
    Deletes working directory
    :param working_dir:  working directory
    """

    try:
        shutil.rmtree(working_dir)
    except Exception as e:
        print ('Can\'t delete %s' % working_dir)
