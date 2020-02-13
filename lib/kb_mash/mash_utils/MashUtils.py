import time
import os
# import errno
import json
import subprocess
import csv
import requests
import time
import traceback
from requests.exceptions import RequestException


def log(message, prefix_newline=False):
    """
    Logging function, provides a hook to suppress or redirect log messages.
    """
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def logerr(err):
    '''
    Log an exception's traceback, split over several lines.
    '''
    for l in traceback.format_exc(err).split('\n'):  # I asssume there's a better way to do this
        log(l)


class SketchException(Exception):
    '''
    General exception thrown when trying to complete a genome comparison
    '''

mash_bin = '/kb/module/mash-Linux64-v2.0/mash'
# ahs_url  = 'https://homology.kbase.us/namespace/%s/search'

class MashUtils:

    def __init__(self, config, auth_token):
        self.scratch = os.path.abspath(config['scratch'])
        if config.get('srv-wiz-url'):
            self.sw_url = config.get('srv-wiz-url')
        else:
            self.sw_url = config['kbase-endpoint'] + '/service_wizard'
        if config.get('id-mapper-url'):
            self.id_mapper_url = config.get('id-mapper-url')
        else:
            self.id_mapper_url = config['kbase-endpoint'] + '/idmapper/api/v1'
        self.auth_token = auth_token
        endpoint = config['kbase-endpoint'].split('/services')[0]
        if 'appdev' in endpoint:
            endpoint = endpoint.replace("appdev", "narrative")
        elif '/kbase' in endpoint:
            endpoint = endpoint.replace("kbase", "narrative.kbase")
        self.endpoint = endpoint

    def mash_sketch(self, genome_file_path, paired_ends=False):
        """
        Generate a sketch file for a given fasta/fastq file path, saving the output to a tempfile.

        Documentation: http://mash.readthedocs.io/en/latest/tutorials.html
        """
        assert os.path.exists(genome_file_path), 'genome_file_path must exist'
        output_path = genome_file_path + '.msh'
        args = [mash_bin, 'sketch', genome_file_path, '-o', output_path]
        if paired_ends:
            # Sketch the reads using `-m 2` to improve results by ignoring single-copy k-mers, which
            # are more likely to be erroneous:
            args = args + ['-m', '2']
        self._run_command(' '.join(args))
        return output_path

    def get_sketch_service_url_with_service_wizard(self):
        '''Get Sketch Serivce url from the service wizard client.'''
        payload = {
            "method":"ServiceWizard.get_service_status",
            "id":'',
            "params":[{"module_name":"sketch_service","version":"beta"}],
            "version":"1.1"
        }

        sw_resp  = requests.post(url=self.sw_url, data=json.dumps(payload))
        sketch_resp = sw_resp.json()
        if sketch_resp.get('error'):
            raise RuntimeError("ServiceWizard Error: "+ str(sketch_resp['error']))
        sketch_url  = sketch_resp['result'][0]['url']

        return sketch_url

    def sketch_service_query(self, input_upas, n_max_results, search_db, cache):
        '''Query assembly homology service to leverage its caching and mash implementation

        params:
            input_upas - list of references to assembly or genome
            n_max_results - number of results to return
            search_db - string to specify search database
            cache - boolean of whether to use cache or not. True = use cache
        '''
        # get current sketch_service url from service wizard
        sketch_url = self.get_sketch_service_url_with_service_wizard()
        results = []
        for inp in input_upas:
            upa = inp['upa']
            input_name = inp['name']
            payload = {
                "method":"get_homologs",
                "params":{
                    'ws_ref':upa,
                    'n_max_results':n_max_results,
                    'search_db': search_db,
                    'bypass_caching': not cache
                }
            }
            print("Processing input: %s" % input_name)
            attempts = 1
            max_attempts = 5
            while True:
                # there's really no good way to test retries w/o a ton of work
                start = time.time()
                try:
                    resp = requests.post(
                        url=sketch_url,
                        data=json.dumps(payload),
                        headers={
                            'content-type': 'application/json-rpc',
                            'Authorization': self.auth_token
                            },
                        timeout=600)
                except RequestException as e:
                    if attempts > max_attempts:
                        # dunno if e.args[0] is right, no docs  
                        # https://requests.readthedocs.io/en/master/api/#exceptions
                        raise SketchException('Sketch service communications failed: ' +
                                              str(e.args[0])) from e
                    log('Sketch service comms failed after ' + str(time.time() - start) + 's')
                    logerr(e)
                else:
                    if not resp.ok:
                        # could try deserializing from json, YAGNI for now
                        if attempts > max_attempts:
                            raise SketchException(
                                'Recieved bad response from the sketch service: {}\n{}'.format(
                                    resp.status_code, resp.text))
                        log('Bad response from Sketch service after ' +
                            str(time.time() - start) + 's')
                        log('Error from sketch service: {}\n{}'.format(
                            resp.status_code, resp.text))
                    else:
                        break
                # possibly multiple runs in succession cause the assembly homology service to
                # time out per Sebastian LeBras. IIUC this was originally fixed by increasing a
                # timeout in nginx, but it's occurring again so try an exponential backoff
                # to let the service CTFO
                time.sleep(2 ** attempts)
                attempts += 1

            if len(input_upas) == 1:
                results = self.parse_results(resp.json())
            else:
                curr = self.parse_results(resp.json(), input_name=input_name, input_id=upa)
                results += curr

        return results

    def parse_results(self, results_data, input_name=None, input_id=None):
        '''
        params:
            results_data: dictionary response from sketch_service
        '''
        if results_data.get('error'):
            raise RuntimeError("Sketch Service Error: " + str(results_data['error']))
        if not results_data.get('result'):
            raise ValueError("No results or results empty in JSON response body")
        if not results_data['result'].get('distances'):
            raise ValueError("No Distances in results JSON response")

        results = []
        distances = results_data['result']['distances']
        # id_to_similarity, id_to_upa, id_to_sciname, id_to_strain = {}, {}, {}, {}
        for d in distances:
            curr = {}
            curr['Id'] = d['sourceid']
            sciname = ""
            if d.get('sciname'):
                sciname += d['sciname']
            if d.get('kbase_id'):
                curr['item_link'] = self.endpoint +  "/#dataview/" + d['kbase_id']
                curr['kbase_id'] = d['kbase_id']
            if d.get('strain'):
                sciname = sciname + " " + d['strain']
            curr['sciname'] = sciname
            curr['dist'] = float(d['dist'])
            if input_name != None:
                curr['input_name'] = input_name
            if input_id != None:
                curr['input_id'] = input_id
            results.append(curr)
        return results

    def id_mapping_query(self, ids):
        """
        """
        payload = {
            "ids":ids
        }
        id_mapper_url = self.id_mapper_url + '/mapping/RefSeq'

        resp = requests.get(url=id_mapper_url, data=json.dumps(payload),\
                    headers={'content-type':"application/json-rpc",'Authorization':self.auth_token})
        return self.parse_mapper_response(resp.json())

    def parse_mapper_response(self, resp):
        """
        """
        if resp.get('error'):
            raise RuntimeError("ID Mapper Error: "+ str(resp.get('error', "unknown error")))

        id_to_upa = {}
        for id_ in resp:
            mappings = resp[id_]["mappings"]
            # default upa is no upa
            id_to_upa[id_] = ""
            for mapping in mappings:
                if mapping['ns'] == "KBase":
                    id_to_upa[id_] = mapping['id']
        return id_to_upa

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """

        log('Start executing command:\n{}'.format(command))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\nOutput:\n{}'.format(exitCode, output))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}'.format(exitCode, output)
            raise ValueError(error_msg)
