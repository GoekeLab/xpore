import argparse
import numpy as np
import pandas
import os
import multiprocessing 
import json
from collections import defaultdict
import csv

from . import helper
from ..diffmod.configurator import Configurator
from ..diffmod.gmm import GMM
from ..diffmod import io
from ..diffmod.statstest import StatsTest

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    
    # Required arguments
    required.add_argument('--config', dest='config', help='yaml configuraion filepath.',required=True)

    # Optional arguments
    optional.add_argument('--n_processes', dest='n_processes', help='number of processes to run.',type=int,default=1)
    optional.add_argument('--save_models', dest='save_models', help='save the models.',default=False,action='store_true') # todo
    optional.add_argument('--resume', dest='resume', help='resume from the previous run.',default=False,action='store_true') 
    
    optional.add_argument('--ids', dest='ids', help='gene ids or transcript ids.',default=[],nargs='*')

    parser._action_groups.append(optional)
    return parser.parse_args()
        
def execute(idx, data_dict, data_info, method, criteria, model_kmer, prior_params, out_paths, save_models,locks):
    """
    Run the model on each posiiton across the given idx.
    """
    data = io.load_data(idx,data_dict,min_count=criteria['readcount_min'],max_count=criteria['readcount_max'],pooling=method['pooling']) 
    models = dict()
    for key,data_at_pos in data.items(): # For each position
        idx, pos, kmer = key
        kmer_signal = {'mean':model_kmer.loc[kmer,'model_mean'],'std':model_kmer.loc[kmer,'model_stdv']}
        kmer_signal['tau'] = 1./(kmer_signal['std']**2)
        y_mean = data_at_pos['y'].mean()
        y_tau = 1./(data_at_pos['y'].std()**2)

        K = 2
        if method['pooling']:
            n_groups = len(data_at_pos['condition_names'])
        else:
            n_groups = len(data_at_pos['run_names'])

        ### Set up priors.
        priors = {'mu_tau':defaultdict(list),'w':dict()}

        for k in range(K):
            priors['mu_tau']['location'] += [kmer_signal['mean']]
            priors['mu_tau']['lambda'] += [prior_params['mu_tau']['lambda'][k]]
            priors['mu_tau']['alpha'] += [kmer_signal['tau']]
            priors['mu_tau']['beta'] += [prior_params['mu_tau']['beta_scale'][k]*1./kmer_signal['tau']]
        
        for k,v in priors['mu_tau'].items():
            priors['mu_tau'][k] = np.array(v)
        
        priors['w']['concentration'] = np.ones([n_groups,K])*1. #GK
        priors['w']['concentration'][:,0] = float(prior_params['w']['concentration'][0])
        priors['w']['concentration'][:,1] = float(prior_params['w']['concentration'][1])
        ###

        ### Fit a model.
        if method['prefiltering']:
            pval = StatsTest(data_at_pos).fit(method=method['prefiltering']['method'])
            if np.isnan(pval) | (pval < method['prefiltering']['threshold']):
                prefiltering = {method['prefiltering']['method']:pval}
                models[key] = GMM(method,data_at_pos,priors=priors,kmer_signal=kmer_signal).fit(), prefiltering
        else:
            models[key] = GMM(method,data_at_pos,priors=priors,kmer_signal=kmer_signal).fit(), None

        
    if save_models & (len(models)>0): #todo: 
        print(out_paths['model_filepath'],idx)
        io.save_models(models,out_paths['model_filepath'])
    if len(models)>0:
        # Generating the result table.
        table = io.generate_result_table(models,data_info)
        with locks['table'], open(out_paths['table'],'a') as f:
            csv.writer(f,delimiter=',').writerows(table)
        # # Logging
        # log_str = '%s: Saving the result table ... Done.' %(idx)
        # with locks['log'], open(out_paths['log'],'a') as f:
        #     f.write(log_str + '\n')
        
    # Logging
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(idx + '\n')
                        
def main():
    args = get_args()
    
    n_processes = args.n_processes       
    config_filepath = args.config
    save_models = args.save_models
    resume = args.resume
    ids = args.ids

    config = Configurator(config_filepath) 
    paths = config.get_paths()
    data_info = config.get_data_info()
    method = config.get_method()
    criteria = config.get_criteria()
    prior_params = config.get_priors()

    print('Using the signal of unmodified RNA from',paths['model_kmer'])
    model_kmer = pandas.read_csv(paths['model_kmer']).set_index('model_kmer')
    ###

    ###
    # Get gene ids for modelling
    # todo
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['model','table','log']:
        out_paths[out_filetype] = os.path.join(paths['out_dir'],'diffmod.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
        
    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Writing the starting of the files.
    ids_done = []
    if resume and os.path.exists(out_paths['log']):
        ids_done = [line.rstrip('\n') for line in open(out_paths['log'],'r')]  
    else:
        with open(out_paths['table'],'w') as f:
            csv.writer(f,delimiter=',').writerow(io.get_result_table_header(data_info,method))
        with open(out_paths['log'],'w') as f:
            f.write(helper.decor_message('diffmod'))


    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=execute,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    ### Load tasks in to task_queue. ###
    f_index,f_data = {},{}
    for condition_name, run_names in data_info.items():
        for run_name, dirpath in run_names.items():
            # Read index files
            df_index = pandas.read_csv(os.path.join(dirpath,'data.index'),sep=',') 
            f_index[run_name] = dict(zip(df_index['idx'],zip(df_index['start'],df_index['end'])))

            # Read readcount files
            # df_readcount[run_name] = pandas.read_csv(os.path.join(info['dirpath'],'readcount.csv')).groupby('gene_id')['n_reads'].sum() # todo: data.readcount

            # Open data files
            f_data[run_name] = open(os.path.join(dirpath,'data.json'),'r') 
    
    # Load tasks into task_queue.
    # gene_ids = helper.get_gene_ids(config.filepath)
#    gene_ids = ['ENSG00000168496','ENSG00000204388','ENSG00000123989','ENSG00000170144'] #test data; todo
    # gene_ids = ['ENSG00000159111']    
    
    if len(ids) == 0:
        ids = helper.get_ids(f_index,data_info)


    print(len(ids),'ids to be testing ...')
    
    for idx in ids:
        if resume and (idx in ids_done):
            continue
        
        data_dict = dict()
        for condition_name, run_names in data_info.items():
            for run_name, dirpath in run_names.items():
                try:
                    pos_start,pos_end = f_index[run_name][idx]
                except KeyError:
                    data_dict[(condition_name,run_name)] = None
                else:
                    # print(idx,run_name,pos_start,pos_end,df_readcount[run_name].loc[idx])
                    f_data[run_name].seek(pos_start,0)
                    json_str = f_data[run_name].read(pos_end-pos_start)
                    # print(json_str[:50])
                    # json_str = '{%s}' %json_str # used for old dataprep
                    data_dict[(condition_name,run_name)] = json.loads(json_str) # A data dict for each gene.
                
        # tmp
        out_paths['model_filepath'] = os.path.join(paths['models'],'%s.hdf5' %idx)
        #
        # if data_dict[run_name][idx] is not None: # todo: remove this line. Fix in dataprep
        task_queue.put((idx, data_dict, data_info, method, criteria, model_kmer, prior_params, out_paths,save_models)) # Blocked if necessary until a free slot is available.

        
    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

    # Close data files
    for f in f_data.values():
        f.close()   
        
    with open(out_paths['log'],'a+') as f:
        f.write(helper.decor_message('successfully finished'))
        

if __name__ == '__main__':
    """
    Usage:
        xpore-diffmod --config CONFIG [--n_processes N_PROCESSES] \
                     [--save_models] [--resume] \
                     [--ids [IDS [IDS ...]]]
    """
    main()
