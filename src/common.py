'''
Common functions
'''
import collections
import logging
import os
import shutil
import tempfile
import uuid

import botocore
import boto3

def downloadFile(srcFile, destFile):
    ''' Download file from S3 '''
    s3 = boto3.resource('s3')
    bucket, key = find_bucket_key(srcFile)
    try:
        s3.meta.client.download_file(bucket, key, destFile)
    except botocore.exceptions.ClientError as e:
        logging.warn(e)

def exists(s3path):
    ''' Return true is s3path is an object, else false '''
    s3 = boto3.resource('s3')
    bucket, key = find_bucket_key(s3path)
    try:
        s3.Object(bucket, key).load()
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == "404":
            # The object does not exist.
            return False
        else:
            # Something else has gone wrong.
            raise
    # The object does exist.
    return True

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
    readSamples reads in a sampleSheetFile consisting of three columns:
    name, bamfile, reference
    :param sampleSheetFile: tab-delimited file
    :return: list of {name, bam, reference} dictionaries
    '''
    if os.path.isfile(sampleSheetFile) is False:
        logging.error("%s does not exist", sampleSheetFile)
        return False
    logging.info("Reading %s", sampleSheetFile)
    sampleNames = []
    samples = []
    with open(sampleSheetFile, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) != 3:
                logging.error("Expected 3 columns in samplesheet, but found %s", len(fields))
                return False

            sampleNames.append(fields[0])
            sample = {"name": fields[0],
                      "bam": fields[1],
                      "reference": fields[2]}
            samples.append(sample)
    duplicates = [item for item, count in collections.Counter(sampleNames).items() if count > 1]
    if duplicates:
        logging.error("Duplicate sampleids found in %s", sampleSheetFile)
        for dup in duplicates:
            logging.error(dup)
        return False
    logging.info("Read %d samples.", len(samples))
    return samples

def generate_working_dir(working_dir_base=None):
    """
    Creates a unique working directory to combat job multitenancy
    :param working_dir_base: base working directory
    :return: a unique subfolder in working_dir_base
    """
    if os.path.isdir('/scratch'):
        tmp_dir = tempfile.mkdtemp(dir="/scratch")
    else:
        tmp_dir = tempfile.mkdtemp()
    logging.debug("Using tmp dir: %s", tmp_dir)
    return tmp_dir

def delete_working_dir(working_dir):
    """
    Deletes working directory
    :param working_dir:  working directory
    """

    try:
        shutil.rmtree(working_dir)
    except OSError:
        print ('Can\'t delete %s', working_dir)
