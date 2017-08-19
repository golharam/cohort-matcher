from __future__ import print_function
import tarfile
import os
import shutil
import uuid

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

def uncompress(tarball, dest_dir):
    """
    Uncompresses a tarbar
    :param tarball: path to tarball to uncompress
    :param dest_dir: destination directory to uncompress to
    """
    tar = tarfile.open(tarball, "r")
    tar.extractall(dest_dir)
    tar.close()
