# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from kb_mash.kb_object_utils.KBObjectUtils import KBObjectUtils
from kb_mash.mash_utils.MashUtils import MashUtils
#END_HEADER


class kb_mash:
    '''
    Module Name:
    kb_mash

    Module Description:
    A KBase module: kb_mash
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/kbaseapps/kb_mash.git"
    GIT_COMMIT_HASH = "ed944a0524a71cc88dc95fea77a6a1d9d5c0eeac"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = os.path.abspath(config['scratch'])
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.auth_token = os.environ['KB_AUTH_TOKEN']
        #self.SEARCH_DBS = {'Ecoli': '/kb/module/test/data/ecolidb.msh',
        #                   'KBaseRefseq': '/data/kb_refseq_ci.msh'}
        #END_CONSTRUCTOR
        pass


    def run_mash_dist_search(self, ctx, params):
        """
        :param params: instance of type "MashParams" -> structure: parameter
           "input_upa" of String, parameter "workspace_name" of
           String, parameter "search_db" of String, parameter "n_max_results"
           of Long
        :returns: instance of type "MashResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: results
        #BEGIN run_mash_dist_search

        if params.get('search_db'):
            search_db = params.get('search_db')
        else:
            raise ValueError("search_db not present in params")
        if params.get('n_max_results'):
            n_max_results = params.get('n_max_results', 10)
        else:
            raise ValueError("n_max_results not present as an argument in params")
        if params.get('input_upa'):
            upa = params.get('input_upa')
        else:
            raise ValueError("Assembly or Genome workspace reference must be specified")

        cache = bool(params.get('cache', 1))

        os.chdir(self.scratch)
        kb_obj_helper = KBObjectUtils(self.config)

        input_upas = kb_obj_helper.input_upa_parse(upa)

        # [file_list] = kb_obj_helper.stage_assembly_files([params['input_upa']])
        # print(file_list)
        mash_utils = MashUtils(self.config, self.auth_token)

        query_results = mash_utils.sketch_service_query(input_upas, n_max_results, search_db, cache)

        report = kb_obj_helper.create_search_report(
            params['workspace_name'],
            query_results,
            False if len(input_upas) == 1 else True
        )

        results = {'report_name': report['name'], 'report_ref': report['ref']}
        #END run_mash_dist_search

        # At some point might do deeper type checking...
        if not isinstance(results, dict):
            raise ValueError('Method run_mash_dist_search return value ' +
                             'results is not type dict as required.')
        # return the results
        return [results]

    def run_mash_sketch(self, ctx, params):
        """
        Generate a sketch file from a fasta/fastq file
        :param params: instance of type "MashSketchParams" (* * Pass in **one
           of** input_path, assembly_ref, or reads_ref *   input_path -
           string - local file path to an input fasta/fastq *   assembly_ref
           - string - workspace reference to an Assembly type *   reads_ref -
           string - workspace reference to a Reads type * Optionally, pass in
           a boolean indicating whether you are using paired-end reads. *
           paired_ends - boolean - whether you are passing in paired ends) ->
           structure: parameter "input_path" of String, parameter
           "assembly_ref" of String, parameter "reads_ref" of String,
           parameter "paired_ends" of type "boolean" (params:
           input_upa: workspace reference to an assembly object
           workspace_name: name of current workspace search_db: database to
           search n_max_results: number of results to return, integer between
           1 and 100)
        :returns: instance of type "MashSketchResults" (* * Returns the local
           scratch file path of the generated sketch file. * Will have the
           extension '.msh') -> structure: parameter "sketch_path" of String
        """
        # ctx is the context object
        # return variables are: results
        #BEGIN run_mash_sketch
        if 'reads_ref' in params:
            reads_utils = ReadsUtils(self.callbackURL)
            result = reads_utils.download_reads({
                'read_libraries': [params['reads_ref']],
                'interleaved': 'true'
            })
            input_path = result['files'][params['reads_ref']]['files']['fwd']
        elif 'assembly_ref' in params:
            assembly_util = AssemblyUtil(self.callbackURL)
            result = assembly_util.get_assembly_as_fasta({'ref': params['assembly_ref']})
            input_path = result['path']
        elif 'input_path' in params:
            input_path = params['input_path']
        else:
            raise ValueError(
                'Invalid params; must provide one of `reads_ref`, `assembly_ref`, or `input_path`.'
            )
        mash_utils = MashUtils(self.config, self.auth_token)
        output_file_path = mash_utils.mash_sketch(input_path, paired_ends=params.get('paired_ends'))
        results = {'sketch_path': output_file_path}
        #END run_mash_sketch

        # At some point might do deeper type checking...
        if not isinstance(results, dict):
            raise ValueError('Method run_mash_sketch return value ' +
                             'results is not type dict as required.')
        # return the results
        return [results]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
