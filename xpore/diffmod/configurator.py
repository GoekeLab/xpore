import yaml
import os
from collections import defaultdict

from ..utils import misc

class Configurator(object):
    def __init__(self, config_filepath):
        self.filepath = os.path.abspath(config_filepath)
        self.filename = self.filepath.split('/')[-1]
        self.yaml = yaml.safe_load(open(self.filepath, 'r'))
        
    def get_paths(self):
        paths = {}
        
        paths['model_kmer'] = self.yaml['paths']['model_kmer']

        paths['out_dir'] = os.path.join(self.yaml['paths']['out_dir'], self.filename)
        paths.update(misc.makedirs(paths['out_dir'],sub_dirs=['models']))
        paths['model_filepath'] = os.path.join(paths['out_dir'], 'models', '%s.model')        
        return paths
        
    def get_data_info(self):
        return self.yaml['data']
    
    def get_criteria(self):
        return self.yaml['criteria']
        
    def get_method(self):
        method = {}
        if 'priors' not in self.yaml.keys():
            method['name'] = 'gmm'
            method['max_iters'] = 500
            method['stopping_criteria'] = 0.0001
            method['compute_elbo'] = True
            method['verbose'] = False
            method['update'] = 'z','y','w','mu_tau'
            method['pooling'] = False
        else:
            pass #todo

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
        