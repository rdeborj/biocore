#!/bin/env python3

"""
This script runs STAR on gzipped FASTQ files, which are found by searching the source directory for files ending in "*R1.fastq.gz"
"""

import re
import sys
import glob
import os
import argparse

import qsub
import pmgctools
import process_bam


def init():
    parser = argparse.ArgumentParser(
        description='This script runs STAR on gzipped FASTQ files, which are found by searching the source '
                    'directory for files ending in "*R1.fastq.gz"')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-f', '--star-fusion', action='store_true', dest='fusion', default=True)
    parser.add_argument('-c', '--chimSegmentMin', required=False, default="20")
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsubdir', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit fusion jobs to.')
    parser.add_argument('-k', '--fusion-queue', default=None, required=False ,dest='fqueue', help='Launch STAR fusion jobs on a different queue')
    parser.add_argument('-i', '--star-index', dest='starindex', help='STAR reference genome')
    parser.add_argument('-g', '--gtf', default=None,  help='GTF file.')
    parser.add_argument('-O', '--one-folder', action='store_true', dest="one_folder", default=False,
                        help='No sample directory structure. All fastq files are in the source directory.')
    args = parser.parse_args()
    return args


def star (sample, outputdir, chimSegmentMin, read1, read2=None, index=None, **other_qsub_options):
    
    if read2 is None: read2 = ''
    tools = ['STAR']
    modules = pmgctools.check_vars(tools)

    if index is None:
        pmgctools.check_vars(['STARINDEX'])
        index = pmgctools.VARS['STARINDEX']

    cmd = 'STAR --runMode  alignReads '
    cmd += "--readFilesIn {}  {} ".format(read1, read2)
    cmd += "--outFileNamePrefix {} ".format(os.path.join(outputdir, sample))
    cmd += "--genomeDir {}  ".format(index)
    cmd += "--runThreadN 4 --chimSegmentMin {} --readFilesCommand zcat --twopassMode Basic --outSAMprimaryFlag AllBestScore ".format(chimSegmentMin)
    cmd += "--outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate "
    cmd += "--quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --genomeSAsparseD 2 \n"

    name = 'star_' + sample
    return qsub.qsub(name, cmd, modules=modules, other='cpu=4|mem=45gb|walltime=72:00:00',  **other_qsub_options)


def star_fusion (sample, outputdir, index, gtf,waitlist, **other_qsub_options):

    tools = ['STAR']
    modules = pmgctools.check_vars(tools)

    if index is None:
        pmgctools.check_vars(['STARINDEX'])
        index = pmgctools.VARS['STARINDEX']

    if gtf is None:
        pmgctools.check_vars(['GTF'])
        gtf = pmgctools.VARS['GTF']
    
    fus_outputdir=outputdir+"_fusion"
    # create output dir if not exists
    if not os.path.exists(fus_outputdir):
        os.makedirs(fus_outputdir, 0o770)

    cmd = "STAR-Fusion --chimeric_junction  {} ".format(os.path.join(outputdir, sample+'Chimeric.out.junction'))
    cmd += "--chimeric_out_sam {} ".format(os.path.join(outputdir, sample+'Chimeric.out.sam'))
    cmd += "--ref_GTF {} ".format(gtf)
    cmd += "--out_prefix {}  \n".format(os.path.join(fus_outputdir, sample))

    name = 'star_fusion_' + sample
    return qsub.qsub(name, cmd, modules=modules, other='cpu=4|mem=8gb|walltime=24:00:00',waitlist=waitlist, **other_qsub_options)

def star_sample(source, outputdir, fusion, gtf, chimSegmentMin, fqueue, samplename=None, index=None, **other_qsub_options):
    if samplename:
        sample_id = samplename
    else:
        matched = re.match('.*Sample_(.*)', source)
        if matched is None:
            print("Error: {} is not standard Sample fastq directory. (Sample_<sample_name>)".format(source),
                  file=sys.stderr)
            exit(1)
        sample_id = matched.group(1)

    waitlist = []

    # create output dir if not exists
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o770)

    fastqlist = glob.glob(os.path.join(source, sample_id + "*R1.fastq.gz")) 
    samplelist = []
    for read1 in fastqlist:
        read2 = read1.replace('R1.fastq', 'R2.fastq')
        if not os.path.exists(read2):
            read2 = None
        samplelist.append(sample)
        job_name = star(sample=sample, outputdir=outputdir, chimSegmentMin=chimSegmentMin, read1=read1, read2=read2, index=index, **other_qsub_options)
        #if dry is true on TORQUE machine ,user need to manually manage dependencies
        if not (not pmgctools.get_var('SGE_ROOT') and  other_qsub_options.get("dry")):
            waitlist.append(job_name)
            print (job_name)
        if fusion:
            if fqueue:
                tmp_queue=other_qsub_options.get("queue")
                other_qsub_options["queue"]=fqueue
            star_fusion(sample=sample, outputdir=outputdir, index=index, gtf=gtf,waitlist=waitlist , **other_qsub_options)
            if fqueue:
                other_qsub_options["queue"]=tmp_queue

    return job_name


if __name__ == '__main__':
    args = init()
	 
    pmgctools.check_vars(['STAR'])
    if not args.one_folder:
        for sample in glob.glob(os.path.join(args.source, 'Sample_*')):
            if os.path.isdir(sample):
                star_sample(source=sample, outputdir=args.output, fusion=args.fusion, gtf=args.gtf, chimSegmentMin=args.chimSegmentMin,
                            fqueue=args.fqueue, index=args.starindex, log=args.log,qsub=args.qsubdir, dry=args.dry, queue=args.queue)
                
    else:
        for fastq in glob.glob(os.path.join(args.source, "*R1.fastq.gz")) :
            sample = os.path.basename(fastq).split('_R1.fastq.gz')[0]
            star_sample(source=args.source, outputdir=args.output, fusion=args.fusion, gtf=args.gtf, chimSegmentMin=args.chimSegmentMin,
                        fqueue=args.fqueue, samplename=sample, index=args.starindex, log=args.log,qsub=args.qsubdir, dry=args.dry,queue=args.queue)
      
