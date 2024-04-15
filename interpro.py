"""--------------------------------------------------------------------
Submit multiple sequences to interproscan and return the domain type
and number

BIOL 595 final project
Morgan Gyger  04/15/2024
-------------------------------------------------------------------"""

import sys
import time
import argparse
from sequence.fasta import Fasta
from api.interpro.interpro import Interpro
import textwrap as _textwrap

def setup_arguments():
    minlen_def = 80
    loglevel_def = 1
    batch_lim_def = 5
    cl = argparse.ArgumentParser(description='InterProScan of sequences')
    cl.add_argument('--logfile',type=argparse.FileType'w',default=sys.stderr,help='Output file for log information')
    cl.add_argument('-m','--minlen',type=int,default=minlen_def,help='Minimum length FASTA to run')
    cl.add_argument('--batch_limit',type=int,default=batch_lim_def,help='Number of sequences to submit per batch')
    cl.add_argument('--loglevel',type=int,default=loglevel_def,help='detail for REST queries report')
    cl.add_argument('fasta_in',type=argparse.FileType('r'))

    return cl.parse_args()

def reformat(job):
    str = ''
    parsed = job.parse_json()
    motifs = parsed['motifs']
    path = parsed['pathway'] #may not need?

    for m in motifs:
        str += '{}\t{}\t{}\n'.format(m['ipr_accession'],
                                     m['src_accession'],
                                     m['description'])
    for p in path:
        str += '{}\t{}\t{}\n'.format(p,path[p]['name'],
                                    path[p]['source'])

    return str

def save_fin(joblist, reformat=None,remove=True):
    delete_list=[]
    for job in joblist:
        if joblist[job] != 'finished':
            continue

        joblist[job]=job.result()

        if reformat:
            text = reformat(job)

        if fh:
            if text:
                fh.write('!{} - {}s\n'.format(job.jobname, job.jobid))
                fh.write('{}\n'.format(text))
            else:
                fh.write('!{} - {} no hits\n'.format(job.jobname, job.jobid))

        if remove:
            delete_list.append(job)

    for job in delete_list:
        del joblist[job]

    return text

"""--------------------------------------------------------------------------------------------------------
Main Code
need to be moved to one central main code?
--------------------------------------------------------------------------------------------------------"""
fasta = Fasta(fh=args.fasta_in)

joblist = {}

template = Interpro(loglevel=1)
template.log_fh = args.logfile
template.email = 'mgyger@purdue.edu'
template.application_select(['Pfam','Panther','SignalP'])
template.output_select = 'json'

sequence_limit = 5
batch_limit = 5
n_sequence = 0
nskip = 60
s = 0

while fasta.next():
    while s < nskip:
        s += 1
        fasta.next()

    if n_sequence >= sequence_limit: break
    n_sequence += 1

    ips = template.clone()
    ips.sequence = fasta.translate().seq.rstrip('*').format(linelen=60)
    ips.jobname = fasta.id
    ips.title = 'FASTA{}'.format(n_sequence)
    joblist[ips] = 'new'

    if ips.run():
        joblist[ips] = 'submitted'

    if n_sequence % batch_limit and n_sequence < sequence_limit:
        continue

    save_fin(joblist,reformat, sys.stdout,True)

args.logfile.close()
exit(0)