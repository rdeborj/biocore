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
        get_sample_name_from_fastq_file(file=None)

    INPUT:
        * file: FASTQ filename

    OUTPUT:
        Returns a string with the sample name extrated.
    """
    # to extract the sample name, we need to remove the following:
    #   * _L00[1-8]_R[12][0-9_]*\.fastq\.gz
    #   * _R[12]_[0-9]*\.fastq.gz
    #.  * _L00[1-8]_R[12]\.fastq.gz
    if re.search("_L[0-8]*_R[12][0-9_]*\.fastq\.gz", file):
        sample = re.sub("_L[0-8_]*_R[12][0-9_]*\.fastq\.gz", "", file)
    elif re.search("_R[12]_[0-9]*\.fastq\.gz", file):
        sample = re.sub("_R[12]_[0-9]*\.fastq\.gz", "", file)
    elif re.search("_L[0-8]*_R[12]\.fastq\.gz", file):
        sample = re.sub("_L[0-8]*_R[12]\.fastq\.gz", "", file)
    elif re.search("_R[12]\.fastq\.gz", file):
        sample = re.sub("_R[12]\.fastq\.gz", "", file)
    else:
        print("File " + file + " is not a paired-end FASTQ file")
        sys.exit()
    return sample


def get_paired_fastq_files(dir=None):
    """
    Parse the files in a directory and extract the read 1 and 2 paired files
    and store as a dictionary with read 1 representing the key and read 2
    representing the corresponding value.

    USAGE:
        get_paired_fastq_files(dir)

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
        if re.search("\.fastq.gz$", _fastq_file):
            _samplename = get_sample_name_from_fastq_file(file=_fastq_file)
            if _samplename in paired_fastq_files.keys():
                if re.search("_R1\.", _fastq_file):
                    _read2 = re.sub("_R1.", "_R2.", _fastq_file)
                    paired_fastq_files[_samplename][_fastq_file] = _read2
                elif re.search("_R1_", _fastq_file):
                    _read2 = re.sub("_R1_", "_R2_", _fastq_file)
                    paired_fastq_files[_samplename][_fastq_file] = _read2
            else:
                if re.search("_R1\.", _fastq_file):
                    _read2 = re.sub("_R1.", "_R2.", _fastq_file)
                    paired_fastq_files[_samplename] = {_fastq_file: _read2}
                elif re.search("_R1_", _fastq_file):
                    _read2 = re.sub("_R1_", "_R2_", _fastq_file)
                    paired_fastq_files[_samplename] = {_fastq_file: _read2}
        else:
            continue
    return paired_fastq_files


def merge_fastq_files(dir='.'):
    """
    Merge multiple read 1 and read 2 FASTQ files into a single read 1 FASTQ
    file and a single read 2 FASTQ file.

    USAGE:
        merge_fastq_files(dir)

    INPUT:
        * dir: full path to the directory containing FASTQ files

    OUTPUT:
        Creates a directory based on the sample name and constructs the command
        to concatenate files of the same sample for read1 and read2.
    """
    read1 = []
    read2 = []
    _fastq_files = get_paired_fastq_files(dir=dir)
    if len(_fastq_files) == 0:
        print("No FASTQ files found...")
        sys.exit()
    else:
        for sample in _fastq_files.keys():
            os.makedirs(sample, exist_ok=True)
            for _fastq_read1 in _fastq_files[sample].keys():
                read1.append(_fastq_read1)
                read2.append(_fastq_files[sample][_fastq_read1])
            read1_files = ' '.join(read1)
            read2_files = ' '.join(read2)
            read1_merged_file = ''.join([sample, '_R1.fastq.gz'])
            read2_merged_file = ''.join([sample, '_R2.fastq.gz'])
            cmd_read1 = ' '.join(["cat", read1_files, '>', '/'.join([sample, read1_merged_file])]) 
            cmd_read2 = ' '.join(["cat", read2_files, '>', '/'.join([sample, read2_merged_file])])
            return_val = {sample: {"cmd_read1":cmd_read1, "cmd_read2":cmd_read2}}
    return return_val
