"""
A Python wrapper for standar functions used by the Bioinformatics Core RNA-Seq
group.
"""

import re
import os
import sys

def get_sample_name_from_fastq_file(file=None):
    """
    Obtain the sample name from a FASTQ file.  A typical FASTQ file will
    consist of the sample name, lane number, indication of whether the file is
    read one or read two for paired end reads and/or the associated file
    number (e.g. KOM22_RNA_T_S6_L001_R1_001.fastq.gz)

    USAGE:
        object.get_sample_name_from_fastq_file(file=None)

    INPUT:
        * file: FASTQ filename

    OUTPUT:
        Returns a string with the sample name extrated.
    """
    # to extract the sample name, we need to remove the following from the:
    #   * _L00[1-8]_R[12][0-9_]*\.fastq\.gz
    return re.sub("[L0-8_]*_R[12][0-9_]*\.fastq\.gz", "", file)

def get_paired_fastq_files(dir=None):
    """
    Parse the files in a directory and extract the read 1 and 2 paired files
    and store as a dictionary with read 1 representing the key and read 2
    representing the corresponding value.

    USAGE:
        object.get_paired_fastq_files(
            dir
            )

    INPUT:
        * dir: the directory containing the FASTQ files

    OUTPUT:
        Returns a dictionary with read 1 as the key and read 2 as the value.
    """
    if dir:
        fastq_files = os.listdir(dir)
    else:
        print("Invalid directory provided")
        sys.exit()

    paired_fastq_files = {}
    for _fastq_file in fastq_files:
        _samplename = get_sample_name_from_fastq_file(file=_fastq_file)
        if _samplename in paired_fastq_files.keys():
            if re.search("_R1\.", _fastq_file):
                _read2 = re.sub("_R1.", "_R2.", _fastq_file)
                paired_fastq_files[_sample_name][_fastq_file] = _read2
            elif re.search("_R1_", _fastq_file):
                _read2 = re.sub("_R1_", "_R2_", _fastq_file)
                paired_fastq_files[_samplename][_fastq_file] = _read2

    return paired_fastq_files


def merge_fastq_files(dir=None):
    """
    Merge multiple read 1 and read 2 FASTQ files into a single read 1 FASTQ
    file and a single read 2 FASTQ file.

    USAGE:
        object.merge_fastq_files(dir)

    INPUT:
        * dir: full path to the directory containing FASTQ files.

    OUTPUT:
        list of return values
    """
    