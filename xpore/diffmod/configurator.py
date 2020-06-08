import yaml
import os
from collections import defaultdict

from ..utils import misc

def get_condition_run_name(condition_name,run_name):
    return '-'.join([condition_name,run_name])

class Configurator(object):
    def __init__(self, config_filepath):
        self.filepath = os.path.abspath(config_filepath)
        self.filename = self.filepath.split('/')[-1]
        self.yaml = yaml.safe_load(open(self.filepath, 'r'))
        
    def get_paths(self):
        paths = {}
        
        if 'prior' in self.yaml:
            paths['model_kmer'] = os.path.abspath(self.yaml['prior'])
        else:
            paths['model_kmer'] = os.path.join(os.path.dirname(__file__),'model_kmer.csv')

        paths['out_dir'] = os.path.join(os.path.abspath(self.yaml['out']))
        paths.update(misc.makedirs(paths['out_dir'],sub_dirs=['models']))
        paths['model_filepath'] = os.path.join(paths['out_dir'], 'models', '%s.model')        
        return paths
        
    def get_data_info(self):
        data = defaultdict(dict)
        for condition_name, run_names in self.yaml['data'].items():
            for run_name, dirpath in run_names.items():
                data[condition_name][get_condition_run_name(condition_name,run_name)] = dirpath
        return data
    
    def get_criteria(self):
        criteria = {}
        if 'criteria' in self.yaml.keys():
            criteria = self.yaml['criteria']
        else:
            criteria['readcount_min'] = 15
            criteria['readcount_max'] = 1000

        return criteria
        
    def get_method(self):
        if 'method' in self.yaml.keys():
            method = self.yaml['method']
        else:
            method = {}

        method.setdefault('name', 'gmm')
        method.setdefault('max_iters', 500)
        method.setdefault('stopping_criteria', 0.0001)
        method.setdefault('compute_elbo', True)
        method.setdefault('verbose', False)
        method.setdefault('update', ['z','y','w','mu_tau'])
        method.setdefault('pooling', False)
        method.setdefault('prefiltering',False)
        return method
    
    def get_priors(self):
        prior_params = defaultdict(dict)
        if 'priors' not in self.yaml.keys():
            # mu_tau
            prior_params['mu_tau']['location'] = ['model_kmer_mean','model_kmer_mean']
            prior_params['mu_tau']['lambda'] = [1,1]
            prior_params['mu_tau']['alpha'] = [0.5,0.5]
            prior_params['mu_tau']['beta'] = ['model_kmer_tau','model_kmer_tau']
            prior_params['mu_tau']['beta_scale'] = [0.5,0.5]

            # w
            prior_params['w']['concentration'] = [0.001,0.001]        
        else:
            pass #todo

        return prior_params
        
