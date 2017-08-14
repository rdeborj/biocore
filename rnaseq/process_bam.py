#!/bin/env python3

"""
  Finds all the BAM in the source directory and run the recommended GATK pre-processing steps: Picard
  MarkDuplicates, realignment around indels, and base quality score recalibration. Options
  and filters chosen based on GATK best practices.
"""

import argparse
import os
import glob

import pmgctools
import qsub



def init():
    parser = argparse.ArgumentParser(
        description='Finds all the BAM in the source directory and run the recommended GATK pre-processing steps')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-i', '--sample_list' , help='Sample list, requiered if STAR output is does not exist yet')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-d', '--dbsnp', help='current dbSNP reference in VCF format')
    parser.add_argument('-k', '--known_1000g', help=', path to known 1000 genome variant , --known_1000g None if non available')
    parser.add_argument('-K', '--known_mills', help=', path to known mills  genome variant , --known_mills None if non available')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')
    options = parser.parse_args()
    return options

def add_rg(source, outputdir, sample, waitlist=None, **other_qsub_options):

    tools = ['picard']
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    rgadd = os.path.join(outputdir, sample +'.readgroup.added.sorted.bam')

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    # picard command changes after 1.118
    if modules['picard'] > '1.118':
        cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $picard_dir/picard.jar AddOrReplaceReadGroups".format(tmpdir)
    else:
        cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $picard_dir/AddOrReplaceReadGroups.jar".format(tmpdir)
    cmd += " INPUT={} OUTPUT={} SORT_ORDER= coordinate RGID={} RGLB={} RGPL={} RGPU={} RGSM={} " \
 	    .format(os.path.join(source, sample + 'Aligned.sortedByCoord.out.bam'), rgadd, sample,"RNA-Seq", "Illumina", "Project", sample)
    cmd += "\nrm -rf {}\n".format(tmpdir)

    return qsub.qsub('add_RG_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=4|mem=20gb|walltime=72:00:00', **other_qsub_options)


def mark_duplicate(source, outputdir, sample, waitlist=None, **other_qsub_options):
    
    tools = ['picard']
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    dedup = os.path.join(outputdir, sample +'.readgroup.added.marked.sorted.bam')
    metrics = os.path.join(outputdir, sample + '.dedup')

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    # picard command changes after 1.118
    if modules['picard'] > '1.118':
        cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $picard_dir/picard.jar MarkDuplicates".format(tmpdir)
    else:
        cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $picard_dir/MarkDuplicates.jar".format(tmpdir)
    cmd += " INPUT={} OUTPUT={} METRICS_FILE={} CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=150000 "\
           " VALIDATION_STRINGENCY=SILENT".format(os.path.join(source, sample +'.readgroup.added.sorted.bam'), dedup, metrics)
    cmd += "\nrm -rf {}\n".format(tmpdir)

    return qsub.qsub('markdup_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=4|mem=20gb|walltime=72:00:00', **other_qsub_options)

def split_ncigar_reads(source, outputdir, sample, waitlist=None, **other_qsub_options):

    tools = ['gatk']
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())
    envs = pmgctools.check_vars(['REF'])
    ref = os.path.abspath(envs["REF"])

    split = os.path.join(outputdir, sample +'.readgroup.added.marked.sorted.split.bam')

    cmd = 'mkdir -p {}\n'.format(tmpdir)


    cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T SplitNCigarReads  -I {}" \
           " -R  {} -o {} -rf ReassignOneMappingQuality -rf UnmappedRead -RMQF 255 -RMQT 60 -U"\
           " ALLOW_N_CIGAR_READS".format(tmpdir,os.path.join(source, sample+'.readgroup.added.marked.sorted.bam'),ref,split)

    cmd += "\nrm -rf {}\n".format(tmpdir)
    
    return qsub.qsub('splitNCigarReads_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=4|mem=20gb|walltime=72:00:00', **other_qsub_options)

def indel_realignment(source, output, sample, known_1000g, known_mills , dbsnp=None, waitlist=None, **other_qsub_options):
    
    tools = ['gatk', 'samtools']
    envs = pmgctools.check_vars(['REF'])

    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    source = os.path.abspath(source)
    tmpdir = pmgctools.tmpdir()

    ref = os.path.abspath(envs["REF"])
    if not dbsnp:
        dbsnp = os.path.abspath(pmgctools.check_vars(["dbSNP"])["dbSNP"])
    if not known_1000g:
        pmgctools.check_vars(['KNOWN_1000G'])
        known_1000g = pmgctools.VARS['KNOWN_1000G'] 
    if not  known_mills:
        pmgctools.check_vars(['KNOWN_MILLS'])
        known_1000g = pmgctools.VARS['KNOWN_MILLS']   
    
    known = "--known {} --known {}".format(known_1000g, known_mills)
    if known_1000g is "None" and known_mills  :
        known = " --known {}".format( known_mills)
    if known_mills is "None" and known_1000g  :
        known = " --known {}".format( known_1000g)
    if known_1000g == "None" and known_1000g == "None" :
        known = "--known {}".format(dbsnp)


    input = os.path.join(source, sample+'.readgroup.added.marked.sorted.split.bam')
    intervals = os.path.join(outputdir, sample +'.intervals')
    realign = os.path.join(source, sample+'.readgroup.added.marked.sorted.split.realigned.bam')

    # Start to work in the directory
    if not os.path.exists(output):
        os.makedirs(output)
    cmd = 'mkdir {}\n'.format(os.path.join(output, tmpdir))
    
    # There is a bug when using the -nt option in the new versions of GATK , RealignerTargetCreator can run faster 
    # when this bug will be fixed
    cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator  -nt 8  -R {}" \
           " -I  {} -o {} {}".format(os.path.join(output, tmpdir), ref, input, intervals, known)
    cmd += "\njava -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner -I {} "\
           "  -o {} -targetIntervals {} -R {} ".format(os.path.join(output, tmpdir), input,realign, intervals, ref)

    # remove intermediate files
    cmd += '\nrm -rf {}'.format(os.path.join(output, tmpdir))
    cmd += '\nrm {}'.format(intervals)

    cmd += '\n'

    return qsub.qsub('realignment_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=8|mem=20gb|walltime=72:00:00', **other_qsub_options)


def bqsr(source, outputdir,  sample, known_1000g, known_mills, dbsnp=None , waitlist=None, **other_qsub_options):
    
    tools = ['gatk']
    envs = pmgctools.check_vars(['REF'])
    ref = os.path.abspath(envs["REF"])
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(source, pmgctools.tmpdir())

    if not known_1000g:
        pmgctools.check_vars(['KNOWN_1000G'])
        known_1000g = pmgctools.VARS['KNOWN_1000G']
    if not  known_mills:
        pmgctools.check_vars(['KNOWN_MILLS'])
        known_1000g = pmgctools.VARS['KNOWN_MILLS']

    known = "-knownSites {} -knownSites {}".format(known_1000g, known_mills)
    if known_1000g is "None" and known_mills  :
        known = " -knownSites {}".format( known_mills)
    if known_mills is "None" and known_1000g  :
        known = " -knownSites {}".format( known_1000g)
    if known_1000g is "None" and known_1000g is "None" :
        known = "-knownSites {}".format(dbsnp)

    input = os.path.join(source, sample+'.readgroup.added.marked.sorted.split.realigned.bam')
    recaldata = os.path.join(source, sample+'.recal_data.table')
    cmd = 'mkdir {}\n'.format(tmpdir)
    cmd += "java -Xmx4g -Djava.io.tmpdir={}  -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -I {} -o {}" \
          " -R {} {}  -U ALLOW_N_CIGAR_READS ".format(tmpdir, input, recaldata, ref,known)
    
    recal = os.path.join(source, sample + '.readgroup.added.marked.sorted.split.realigned.recal.bam')
    cmd += "\njava -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads -nct 8 -I {} -R {}" \
           " -BQSR {} -o {} ".format(tmpdir, input, ref, recaldata, recal)
    cmd += '\nrm {} \n'.format(recaldata)
    cmd += '\nrm -rf {}'.format(tmpdir)

    return qsub.qsub('bqsr_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=8|mem=12gb|walltime=72:00:00', **other_qsub_options)


def process_bam(source, outputdir, sample, known_1000g, known_mills, dbsnp=None, waitlist=None, **other_qsub_options):
    
    wait_add = []
    wait_add.append(add_rg(source, outputdir, sample, waitlist=wait_add, **other_qsub_options))
    source=outputdir

    wait_dup = []
    wait_dup.append(mark_duplicate(source, outputdir, sample, waitlist=wait_add, **other_qsub_options))
    
    wait_split = []
    wait_split.append(split_ncigar_reads(source, outputdir, sample, waitlist=wait_dup, **other_qsub_options))

    indel_wait = []
    indel_wait.append(indel_realignment(source, outputdir, sample, known_1000g, known_mills, dbsnp, waitlist=wait_split, **other_qsub_options))
    
    bqsr_wait = []
    bqsr_wait.append(bqsr(source, outputdir, sample, known_1000g, known_mills, dbsnp, waitlist=indel_wait, **other_qsub_options))

    return bqsr_wait + indel_wait


if __name__ == '__main__':
    args = init()
    if args.ini:
        pmgctools.read_vars(args.ini)
    pmgctools.check_vars(['gatk', 'samtools', 'picard', 'REF','GTF'])
    if not args.dbsnp:
            pmgctools.check_vars(['dbSNP'])
    if not args.known_1000g:
            pmgctools.check_vars(['KNOWN_1000G'])
    if not args.known_mills:
            pmgctools.check_vars(['KNOWN_MILLS'])

    source, outputdir, sample_list, log, qsubdir, dry, dbsnp, known_1000g , known_mills = \
    args.source,args.output, args.sample_list, args.log, args.qsub, args.dry, args.dbsnp, args.known_1000g, args.known_mills


    if sample_list is None:
        bams = glob.glob(os.path.join(source, '*Aligned.sortedByCoord.out.bam'))
        for bamfile in bams:
            sample = os.path.basename(bamfile).split("Aligned.sortedByCoord.out.bam")[0]
            process_bam(source=source, outputdir=outputdir, sample=sample, known_1000g=known_1000g, known_mills=known_mills, dbsnp=dbsnp,
                        log=log, qsub=qsubdir, dry=dry, queue=args.queue)

    else:
        with open(sample_list) as f:
            for line in f:
                sample = line.rstrip()
                process_bam(source=source, outputdir=outputdir, sample=sample, 
                             known_1000g=known_1000g, known_mills=known_mills, dbsnp=dbsnp, log=log, qsub=qsubdir, dry=dry, queue=args.queue)
