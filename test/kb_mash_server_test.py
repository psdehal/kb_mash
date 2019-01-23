# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import time
import shutil

from ReadsUtils.ReadsUtilsClient import ReadsUtils

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from biokbase.workspace.client import Workspace as workspaceService
from kb_mash.kb_mashImpl import kb_mash
from kb_mash.kb_mashServer import MethodContext
from kb_mash.authclient import KBaseAuth as _KBaseAuth
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil


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
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
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

    def get_genome_ref(self, ws_name):
        if hasattr(self.__class__, 'genomeInfo'):
            return self.__class__.genomeInfo
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        tf = 'ecoliMG1655.fa'
        target = os.path.join(self.scratch, tf)
        self.genome_path = target
        shutil.copy('data/' + tf, target)
        self.__class__.genomeInfo = au.save_assembly_from_fasta({
            'file': {'path': target},
            'workspace_name': ws_name,
            'assembly_name': 'ecoliMG1655'
        })
        return self.__class__.genomeInfo

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_mash_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_mash_search(self):
        ws_name = self.getWsName()
        params = {'input_assembly_upa': self.get_genome_ref(ws_name), 'workspace_name': ws_name,
                  'search_db':'NCBI_Refseq', 'n_max_results':10}
        self.getImpl().run_mash_dist_search(self.getContext(), params)

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
