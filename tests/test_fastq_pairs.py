#!/usr/bin/env python3


from core.general import get_paired_fastq_files

paired_fastq_files = get_paired_fastq_files(dir="./tests")
print(paired_fastq_files)
