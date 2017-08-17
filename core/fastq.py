"""
A set of functions for use with FASTQ files.
"""

import os
import re

def concat_fastq_files(dir):
    """
    Concatenate a set of FASTQ files for a single sample.

    USAGE:
        object.concat_fastq_files(
            dir
            )

    INPUT:
        * dir: a directory containing FASTQ files

    OUTPUT:
        Returns a dictionary of samples and corresponding FASTQ files
    """
    files = os.listdir(dir)
