import sys
import json
import requests
import time

class Interpro():
    available = {'applications': ['NCBIfam', 'SFLD', 'Phobius', 'SignalP', 'SignalP_EUK',
                                  'SignalP_GRAM_POSITIVE', 'SignalP_GRAM_NEGATIVE', 'SUPERFAMILY',
                                  'Panther', 'Gene3d', 'HAMAP', 'PrositeProfiles',
                                  'PrositePatterns', 'Coils', 'SMART', 'CDD', 'PRINTS', 'PfamA',
                                  'MobiDBLite', 'PIRSF', 'TMHMM', 'FunFam', 'PIRSR'],
                 'commands': ['run', 'status', 'result'],
                 'outputs': ['out', 'log', 'tsv', 'xml', 'gff', 'json',
                             'htmltarball', 'sequence', 'submission']
                 }
    def __init__(self):
        self.email = ''
        self.title = ''
        self.outputfile = ''
        self.applications = []
        self.sequence = ''
        self.output = 'json'
        self.parameters = {}

        self.url = u'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
        self.jobid = ''
        self.jobname = ''
        self.jobstatus = ''

        self.response = None
        self.content = ''

    def pick_app(self,selected,keep=False):
        if not keep:
            self.applications = []
        for app in selected:
            if app == 'Pfam':
                app = 'PfamA'
            if app in Interpro.available['applications']:
                self.applications.append(app)
            else:
                self.message = {'type':     'not available',
                                'text':     f'application={app}',
                                'loglevel':2}
        return len(self.applications)

    def pick_out(self, selected):
        if selected in Interpro.available['outputs']:
            self.output = selected
        else:
            self.message = {'type': 'not_available',
                            'text': f'output={selected}',
                            'loglevel': 2}

            return False

        return True

    def parse_json(self):
        pjson = json.loads(self.content)
        matches = pjson['results'][0]['matches']
        path_all= {}
        motifs = []
        for m in matches:
            entry = m['signature']['entry']
            signature = m['signature']['signatureLibraryRelease']
            source_accession = m['signature']['accession']
            source = '{} {}'.format(signature['library'],
                                    signature['version'])

            if not entry:
                continue

            motifs.append({'ipr_accession':entry['accession'],
                           'src_accession':source_accession,
                           'description':  entry['description'] or ''})

            if 'pathwayXRefs' in entry:
                pathstr = ''
                for path in entry['pathwayXRefs']:
                    pathstr += '{} ({})'.format(path['id'], path['name'])
                    if path['databaseName'] == 'Reactome':
                        field = path['id'].split('-')
                        id = 'Reactome:{}'.format(field[2])
                    elif path['databaseName'] == 'MetaCyc':
                        id = 'Metacyc:{}'.format(path['id'])
                    elif path['databaseName'] == 'KEGG':
                        id = 'KEGG:{}'.format(path['id'])
                    else:
                        print('unknown pathway {} | {} | {}'.format(path['databaseName'],
                                                                    path['id'], path['name']))

                    if id in path_all:
                        if source_accession not in path_all[id]['source']:
                            path_all[id]['source'].append(source_accession)
                    else:
                        path_all[id] = {'name': path['name'],
                                        'source': [source_accession]}

            return {'jobname': self.jobname, 'motifs': motifs, 'pathway': path_all}

   def submit(self, show_query=False):
        is_success = False

        # general fields for all queries
        param = {u'email': self.email, u'title':self.title, u'sequence':self.sequence,
                 u'output':self.output}

        if self.applications:
            # add selected applications
            param['appl'] = ','.join(self.applications)

        if self.parameters:
            for para in self.parameters:
                param[para] = self.parameters[para]

        command = self.url + 'run'
        try:
            self.response = requests.post(command, files=param, headers={'User-Agent':'ips-client'})
        except self.ConnectionError:
            sys.error.write( f'Connection Error: {command}')
            return False

        if show_query:
            # print out query if requested
            print(self.response.request.headers, '\n')
            print(self.response.request.body, '\n')

        if self.response_is_error('submitting job'):
            self.jobstatus = 'failed'
        else:
            # success
            self.jobid = self.response.text
            self.jobstatus = 'submitted'
            self.message = {'type':    'submitted',
                            'text':    f'job_name={self.title};job_id={self.jobid}',
                            'loglevel':1}

            is_success = True

        return is_success
    def result(self):
        # get the final result
        command = self.url + 'result/' + self.jobid + '/' + self.output
        self.response = requests.get(command)
        if not self.response_is_error('retrieving result'):
            # success
            self.content = self.response.text
            self.message = {'type':    'retrieved',
                            'text':    f'job_id={self.jobid};output_len={len(self.output)}',
                            'loglevel':1}
            return 'retrieved'

        return ''

    def response_is_error(self, task):

        if self.response.status_code == 200:
            # success
            is_error = False

        else:
            # error
            is_error = True
            self.jobstatus = 'error'
            self.message = {'type':    task,
                            'text':    f'job_id={self.jobid};status={self.response.status_code}',
                            'loglevel':1}

        return is_error
