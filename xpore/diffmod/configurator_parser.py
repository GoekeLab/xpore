import configparser
import os
from collections import defaultdict
from ..utils import misc


def read_config_file(config_file):
    config = configparser.ConfigParser()
    config.optionxform = str  # Setting to str would make option names case sensitive:
    config.read(config_file)

    return config


def config_section_map(config, section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


class Configurator(object):
    def __init__(self, config_filepath):
        self.filepath = os.path.abspath(config_filepath)
        self.filename = self.filepath.split('/')[-1]
        self.config = read_config_file(self.filepath)

    def get_info(self):
        """ Return basic information of the datasets specified in the configuration file. """
        
        info = dict()
        info['run_names'] = config_section_map(self.config, 'Info')['run_names'].split(',')
        info['condition_names'] = config_section_map(self.config, 'Info')['condition_names'].split(',')

        info['cond2run_dict'] = defaultdict(list)
        for condition_name, run_name in zip(info['condition_names'], info['run_names']):
            info['cond2run_dict'][condition_name]  += [run_name]

        return info

    def get_paths(self):
        """ Return all file paths specified in the configuration file. """

        paths = dict()

        paths['data_dirs'] = config_section_map(self.config, 'Path')['data_dirs'].split()
        paths['out_dir'] = os.path.join(config_section_map(self.config, 'Path')['out_dir'], self.filename)
        paths.update(misc.makedirs(paths['out_dir'],sub_dirs=['models']))
        paths['model_filepath'] = os.path.join(paths['out_dir'], 'models', '%s.model')
        # paths['result_filepath'] = os.path.join(paths['out_dir'], 'out', '%s.table')
        # paths['gt_mapping_filepath'] = config_section_map(self.config, 'Path')['gt_mapping_filepath']
        # paths['gt_mapping_dir'] = config_section_map(self.config, 'Path')['gt_mapping_dir']
        # paths['bamtx_dir'] = config_section_map(self.config, 'Path')['bamtx_dir']
        # paths['bamgenome_dir'] = config_section_map(self.config, 'Path')['bamgenome_dir']
        paths['model_kmer'] = config_section_map(self.config, 'Path')['model_kmer_filepath']
        # paths['modification_filepath'] = config_section_map(self.config, 'Path')['modification_filepath']
        # paths['isoformTSSTES_filepath'] = config_section_map(self.config, 'Path')['isoformTSSTES_filepath']

        return paths

    def get_method(self):
        method = dict()
        method['name'] = config_section_map(self.config, 'Method')['name']

        if 'gmm' in method['name']:
            method['max_iters'] = int(config_section_map(self.config, 'Method')['max_iters'])
            method['stopping_criteria'] = float(config_section_map(self.config, 'Method')['stopping_criteria'])
            method['pooling'] = config_section_map(self.config,'Method')['pooling'] == 'True'
            # method['compute_elbo'] = config_section_map(self.config,'Method')['compute_elbo'] == 'True'
            # method['config']['verbose'] = config_section_map(self.config,'Method')['verbose'] == 'True'
        return method

    def get_priors(self):
        prior_params = {}
        # mu_tau
        prior_params['location'] = config_section_map(self.config, 'Prior')['location'].split(',')
        prior_params['lambda'] = list(map(float, config_section_map(self.config, 'Prior')['lambda'].split(',')))
        prior_params['alpha'] = config_section_map(self.config, 'Prior')['alpha'].split(',')
        prior_params['beta'] = config_section_map(self.config, 'Prior')['beta'].split(',')
        prior_params['beta_scale'] = list(map(float, config_section_map(self.config, 'Prior')['beta_scale'].split(',')))

        # w
        prior_params['concentration'] = config_section_map(self.config, 'Prior')['concentraion'].split(',')

        return prior_params

    def get_criteria(self):
        criteria = dict()
        criteria['read_count_min'] = int(config_section_map(self.config, 'Criteria')['read_count_min'])
        criteria['read_count_max'] = int(config_section_map(self.config, 'Criteria')['read_count_max'])
        return criteria
