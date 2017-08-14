"""
A Python wrapper for the Bioinformatics Core RNA-Seq standard pipeline.
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
    return re.sub("_L00[1-8]_R[12][0-9_]*\.fastq\.gz", "", file)

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
        Returns a dictionary with read 1 as the key and read as the value.
    """
    if dir:
        fastq_files = os.listdir(dir)
    else:
        print("Invalid directory provided")
        sys.exit()

    for fastq_file in fastq_files:
        if re.search("_R1[._][0-9.]*fastq.gz", fastq_file):
            if (re.sub())
            paired_fastq{fastq_file:}