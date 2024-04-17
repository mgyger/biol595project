"""--------------------------------------------------------------------
Submit multiple sequences to interproscan and return the domain type
and number

BIOL 595 final project
Morgan Gyger  04/15/2024
-------------------------------------------------------------------"""
import argparse
import sys
import pickle
from os.path import basename, exists
from os import mkdir
from interpro.class_interpro import Interpro

def process_command_line():
    cl = argparse.ArgumentParser(
        description='Run interproscan on selected orthogroups',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-g', '--orthogroup',
                    help='Orthogroups file',
                    type=str,
                    default='<not specified>')

    cl.add_argument('-o', '--out',
                    help='output directory for results',
                    type=str,
                    default='./')

    cl.add_argument('-s', '--skip',
                    help='skip (default: %(default)s)',
                    action='store_true',
                    default=False)

    return cl.parse_args()

def read_fasta():
    ###### pull and separate the NCBI results

'''ALL BELOW NEEDS TO BE ADJUSTED FOR NCBI AND OUR FASTA FORMATTING'''
def interpro_setup(fasta, og, pickledir):
    ips = Interpro()
    ips.email = 'mgyger@purdue.edu'
    ips.application_select(['PfamA', 'SMART', 'PrositePatterns', 'CDD', 'NCBIfam', 'PIRSF', 'SuperFamily'])
    ips.output_select('json')
    ips.parameter_select({'goterms': True, 'pathways': False})
    ips.title = f'{og}_{fasta.id}'
    ####fasta.seq = fasta.seq.rstrip('*')
    ####ips.sequence = fasta.format()
    ips.outputfile = pickledir + ips.title.replace('|', '_') + '.pkl'

class InterproscanManager:
    n = 0

    def __init__(self, opt, batch_limit=30, poll_time=30, poll_max=50, pkl='pkl'):
        self.submitted = []
        self.finished = []
        self.failed = []
        self.save = []
        self.opt = opt
        self.batch_limit = batch_limit
        self.poll_time = poll_time
        self.poll_max = poll_max
        self.pkl = f'{pkl}/'
        if not exists(pkl):
            mkdir(pkl)
        self.log_fh = InterproscanManager.getlog('og_interpro.log')

    def submit(self, fasta, og):
        ip_submitted = self.submitted
        ip_failed = self.failed
        if len(ip_submitted) < self.batch_limit:
            ips = interpro_setup(fasta, og, self.pkl)

            if self.opt.skip and exists(ips.outputfile):
                # this sequence exists in the output, skip
                print(f'\tskipping {og}:{fasta.id}')
                return False

            print(f'\tJob {self.n} - {og}:{fasta.id}')
            self.n += 1  # total number of jobs submitted (class variable)

            tries = 1
            success = ips.submit(show_query=False)
            while not success and tries < 3:
                # try three times to submit with 5 seconds between tries
                time.sleep(5)
                success = ips.submit(show_query=False)
                tries += 1

            if success:
                # success
                self.log('SUBMIT', ips.title)
                ip_submitted.append(ips)
            else:
                # failure
                ip_failed.append(ips)
                sys.stderr.write(f'{ips.title} failed\n')
                self.log('FAIL-SUB', ips.title)

        return True

    def getresult(self):
        if not finished:
            # sys.stderr.write('no finished jobs to retrieve\n')
            return

        # parse finished jobs and extract desired information
        while finished:
            thisjob = finished.pop()
            self.log('RETRIEVE', thisjob.title)
            thisjob.result()
            self.save.append(thisjob)
            self.save_as_pickle(thisjob)
            # print(thisjob.content)

        # a cooldown period appears to be necessary
        # time.sleep(75)
        return len(self.save)

