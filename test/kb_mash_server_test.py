# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import time
import shutil
import json
from installed_clients.ReadsUtilsClient import ReadsUtils

from os import environ

from configparser import ConfigParser

from Workspace.WorkspaceClient import Workspace as workspaceService
from kb_mash.kb_mashImpl import kb_mash
from kb_mash.kb_mashServer import MethodContext
from kb_mash.authclient import KBaseAuth as _KBaseAuth
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil

class kb_mashTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_mash'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        # authServiceUrl = cls.cfg['auth-service-url']
        # auth_client = _KBaseAuth(authServiceUrl)
        # user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        # 'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_mash',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        # suffix = int(time.time() * 1000)
        # wsName = "test_kb_mash_" + str(suffix)
        # cls.ws_info = cls.wsClient.create_workspace({'workspace': wsName})
        cls.serviceImpl = kb_mash(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        # cls.au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        # cls.setup_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def get_genome_ref(self, ws_name, tf='ecoliMG1655.fa'):
        if hasattr(self.__class__, 'genomeInfo'):
            return self.__class__.genomeInfo
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        target = os.path.join(self.scratch, tf)
        self.genome_path = target
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        shutil.copy(os.path.join(curr_dir, 'data', tf), target)
        self.__class__.genomeInfo = au.save_assembly_from_fasta({
            'file': {'path': target},
            'workspace_name': ws_name,
            'assembly_name': tf.split('.fa')[0]
        })
        return self.__class__.genomeInfo

    # def get_genome_set_ref(self, ws_name, genome_refs):
    #     items = []
    #     for ref in genome_refs:
    #         item = {
    #             "ws_genome_id": ref,
    #             "label":""}
    #         items.append(item)

    #     genome_set = {
    #         "description": "",
    #         "items":items
    #     }
    #     dfu  = DataFileUtil(self.__class__.callback_url)
    #     info = dfu.save_objects({
    #         "objects": [{
    #             "data":genome_set,
    #             "type":"KBaseSets.GenomeSet"
    #         }],
    #         "id": self.__class__.wsId
    #     })[0]
    #     gs_ref = '/'.join([str(info[6]),str(info[0]),str(info[4])])
    #     return gs_ref

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_mash_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})# noqa
        self.__class__.wsName = wsName
        # print('-'*80)
        # print(ret)
        # print('-'*80)
        self.__class__.wsId = ret[0]
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_mash_search_genome_set(self):
        ws_name = self.getWsName()
        # genome_refs = ['38012/42/1','38012/38/1','38012/33/1','38012/36/2']

        params = {'input_upa': "38012/45/1", 'workspace_name': ws_name,
                 'search_db': "NCBI_Refseq", 'n_max_results': 10}
        self.getImpl().run_mash_dist_search(self.getContext(), params)

    def test_another_genomeset(self):
        ws_name = self.getWsName()
        genome_set_ref = "38485/7/3"
        params = {"input_upa": genome_set_ref, "workspace_name": ws_name,
                  "search_db": "NCBI_Refseq", "n_max_results": 100}
        ret = self.getImpl().run_mash_dist_search(self.getContext(), params)

    def test_mash_search_jgi(self):
        ws_name = self.getWsName()
        params = {'input_upa': self.get_genome_ref(ws_name, tf='3300011599_1.fa'), 'workspace_name': ws_name,
                  'search_db':'JGI_MAGS', 'n_max_results':10}
        self.getImpl().run_mash_dist_search(self.getContext(), params)

    def test_mash_search(self):
        ws_name = self.getWsName()
        params = {'input_upa': self.get_genome_ref(ws_name), 'workspace_name': ws_name,
                  'search_db':'NCBI_Refseq', 'n_max_results':10}
        self.getImpl().run_mash_dist_search(self.getContext(), params)

    def test_mash_search_with_caching(self):
        ws_name = self.getWsName()
        params = {
            'input_upa': self.get_genome_ref(ws_name),
            'workspace_name': ws_name,
            'search_db':'NCBI_Refseq',
            'n_max_results':10,
            'cache': 1
        }
        ret = self.getImpl().run_mash_dist_search(self.getContext(), params)

    def test_mash_sketch_valid_local(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        params = {'input_path': os.path.join(dir_path, 'data', 'ecoliMG1655.fa')}
        self.getImpl().run_mash_sketch(self.getContext(), params)
        output_path = os.path.join(dir_path, 'data', 'ecoliMG1655.fa.msh')
        with open(output_path, 'rb') as output_file:
            num_lines = sum(1 for line in output_file)
        self.assertTrue(os.path.exists(output_path))
        self.assertEqual(num_lines, 103)

    def test_mash_sketch_valid_assembly_ref(self):
        ws_name = self.getWsName()
        assembly_ref = self.get_genome_ref(ws_name)
        params = {'assembly_ref': assembly_ref}
        self.getImpl().run_mash_sketch(self.getContext(), params)
        output_path = os.path.join(self.scratch, 'ecoliMG1655.fa.msh')
        with open(output_path, 'rb') as output_file:
            num_lines = sum(1 for line in output_file)
        self.assertTrue(os.path.exists(output_path))
        self.assertEqual(num_lines, 103)

    def test_mash_sketch_valid_reads_ref(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        reads_file_name = 'reads-example.fastq'
        reads_test_path = os.path.join(dir_path, 'data', reads_file_name)
        reads_scratch_path = os.path.join(self.scratch, reads_file_name)
        shutil.copy(reads_test_path, reads_scratch_path)
        reads_utils = ReadsUtils(self.callback_url)
        upload_result = reads_utils.upload_reads({
            'wsname': self.getWsName(),
            'interleaved': 'true',
            'fwd_file': reads_scratch_path,
            'name': 'example-reads',
            'sequencing_tech': 'illumina'
        })
        reads_ref = upload_result['obj_ref']
        params = {'reads_ref': reads_ref, 'paired_ends': True}
        result = self.getImpl().run_mash_sketch(self.getContext(), params)
        output_path = result[0]['sketch_path']
        with open(output_path, 'rb') as output_file:
            num_lines = sum(1 for line in output_file)
        self.assertTrue(os.path.exists(output_path))
        self.assertEqual(num_lines, 25)
