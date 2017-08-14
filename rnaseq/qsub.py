#!/bin/env python3

"""
This module is a python wrapper of qsub to submit jobs to the cluster.
"""

import argparse
import os
import sys
import time
import pmgctools
import subprocess
from jinja2 import Environment, FileSystemLoader
import logging

def init():
    """
    Initialize the arguments

    Input:
        NA

    Output:
        Returns the parsed arguments.
    """
    parser = argparse.ArgumentParser(description='Submit command to the cluster.')
    parser.add_argument('command', help='Command you want to run on the cluster.')
    parser.add_argument('-Q', '--queue', help='Cluster queue you want to submit to.')
    parser.add_argument('-l', '--log', default='process.log', help='Log file.')
    parser.add_argument('-q', '--qsub', default='qsub', help='Directory to store the qusb script.')
    parser.add_argument('-e', '--engine', default='SGE', choices=['SGE', 'Torque'], help='Grid engine (SGE or Torque).')
    parser.add_argument('-i', action='store_true', dest='no_overwrite', default=False,
                        help='Not overwrite the qsub script.')
    parser.add_argument('-n', '--name', required=True, help='Name of the job.')
    parser.add_argument('-m', '--module', action='append', default=[],
                        help='modules to load, use format "tool=version".')
    parser.add_argument('-w', '--wait', help='wait for jobs to finish before launch this job.(separate by ",")')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False, help='Dry run.')
    parser.add_argument('-o', '--other-options', dest='other', help='Other qsub options.')

    return parser.parse_args()


def waitlist_to_str(waitlist, spliter=','):
    """
    Create a comma separated string with job IDs to wait for.  This will be used
    with grid engine job dependencies.

    Input:
        * waitlist: a list of job IDs
        * spliter: delimiter to be used in the return string (default: ,)

    Output:
        Returns a string with delimited job IDs
    """
    if isinstance(waitlist, str):
        return waitlist
    else:
        wait = []
        for i in waitlist:
            if isinstance(i, list):
                wait.append(waitlist_to_str(i, spliter))
            else:
                wait.append(i)

        return ','.join(wait)

def create_qsub_template(name, engine, queue=None, waitlist=None, other=None, must_run=False):
    """
    Create a shell script to submit to a high performance compute cluster using
    Torque PBS (HPC for Health) or SGE (UHN local HPC).

    """
    waitlist = waitlist_to_str(waitlist) if waitlist else None

    if engine == 'SGE':
        qsub_content = '''\
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N {}
'''.format(name)

        if queue:
            qsub_content += '#$ -q {}\n'.format(queue)

        if waitlist:
            qsub_content += '#$ -hold_jid {}\n'.format(waitlist)


    else:
        qsub_content = '''\
#!/bin/bash
#PBS -S /bin/bash
#PBS -N {}
'''.format(name)

        if queue:
            qsub_content += '#PBS -q {}\n'.format(queue)

        if waitlist:
            if must_run:
                qsub_content += '#PBS -W depend=afterany:{}\n'.format(waitlist.replace(',', ':'))
            else:
                qsub_content += '#PBS -W depend=afterok:{}\n'.format(waitlist.replace(',', ':'))

        if other:
            for i in other.split('|'):
                k, v = i.split('=')
                k = k.strip()
                v = v.strip()
                if k == 'cpu':
                    qsub_content += '#PBS -l nodes=1:ppn={}\n'.format(v)
                if k == 'mem':
                    qsub_content += '#PBS -l vmem={}\n'.format(v)
                if k == 'walltime':
                    qsub_content += '#PBS -l walltime={}\n'.format(v)

        qsub_content += 'cd $PBS_O_WORKDIR\n'

    return qsub_content


def qsub(name, command, queue=None, log='process.log', meg='', qsub='qsub', no_overwrite=False, modules=None,
         waitlist=None, dry=False, other=None, must_run=False):
    """

    Input:
        * name:
            name of the job, also used as name of the script file
        * command:
            command to qsub
        * queue:
            queue to submit the job, None=default queue
        * log:
            markdown file to record the command
        * msg:
            extra message to write into log file
        * qsub:
            directory to store qsub script
        * engine:
            SGE or Torgue
        * no_overwrite:
            if True and there is already a script with the same name, error
            will raise
        * modules:
            a dictionary containing modules and corresponding version to be
            loaded in the submitted job shell
        * waitlist:
            a comma separated string of job IDs to wait for
        * dry:
            dry-run
        * other:
            other qsub options

    Return:
        Returns the jobname if the job is submitted to SGE or
        the extracted job ID if job is submitted to PBS.
    """

    if pmgctools.get_var('SGE_ROOT'):
        engine = 'SGE'
    else:
        engine = 'Torque'

    if not os.path.exists(qsub):
        os.makedirs(qsub)

    script = os.path.join(qsub, name + '.sh')
    if os.path.exists(script) and no_overwrite:
        sys.stderr.write('Error: qsub script {} already exists.\n'.format(script))

    cmds = command.split('\n')

    if log is not None and log != 'None':
        with open(log, 'a') as logfile:
            logfile.write('**' + name + '**  \n')
            logfile.write('Submitted at: {}  \n'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
            logfile.write(meg + '\n')
            if modules is not None:
                logfile.write("Module loaded:  \n")
                for key, value in modules.items():
                    logfile.write('- {}/{}\n'.format(key, value))

            logfile.write("\n```\n")
            for line in cmds:
                if not line.startswith('mkdir') and not line.startswith('rm'):
                    logfile.write(line + '\n')
            logfile.write('\n```\n\n')

    qsub_content = create_qsub_template(name=name, engine=engine, queue=queue, waitlist=waitlist, must_run=must_run, other=other)

    exclude_cmds = ['rm', 'mkdir', 'cd', 'echo', 'mv', 'if', 'export', 'let', 'exit', 'sleep']
    cmds_file = open('cmds_submitted', 'a')
    with open(script, 'w') as scriptfile:
        scriptfile.write(qsub_content)

        if modules is not None:
            for key, value in modules.items():
                scriptfile.write('module load {}/{}\n'.format(key, value))

        #scriptfile.write(command + '\n')
        if waitlist:
            scriptfile.write('sleep 120\n')
        for line in cmds:
            scriptfile.write(line + '\n')
            if len(line) > 0 and not line.split()[0] in exclude_cmds:
                scriptfile.write('exit_code=$?;final_exit_code=$((final_exit_code|exit_code))\n')
                scriptfile.write("echo \"TIME: $(date); JOB NAME: {}; COMMAND: '{}'; EXIT STATUS: $exit_code\" >> {}\n".format(name, line.replace("'", "\""), os.path.join(pmgctools.getwd(), 'exit_status.txt')))
                #scriptfile.write("echo EXIT STATUS: $?\n")
                print(line, file=cmds_file)
                pmgctools.cmds_more()

        scriptfile.write('exit $final_exit_code\n')

    cmds_file.close()

    if not dry:
        job_id = subprocess.check_output('qsub {}'.format(script).split())
        print(job_id.decode("utf-8").strip())
    else:
        job_id = b''

    return name if engine=='SGE' else job_id.decode("utf-8").strip()


if __name__ == '__main__':

    args = init()

    modules = dict()
    for i in args.module:
        for j in i.split(','):
            if '=' not in j:
                sys.stderr.write("Error: can not parse {}\n".format(j))
                continue
            else:
                tool, version = j.split('=')
                modules[tool] = version

    qsub(args.name, args.command, args.queue, args.log, args.qsub, args.engine, args.no_overwrite, modules, args.wait,
         args.dry, args.other)
