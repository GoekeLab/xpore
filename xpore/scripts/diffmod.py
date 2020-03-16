import argparse
import numpy
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

def get_args():
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('--config', dest='config', help='Specify.',required=True)

    # Optional arguments
    parser.add_argument('--n_processes', dest='n_processes', help='Specify the number of processes.',type=int,default=1)
    parser.add_argument('--save_table', dest='save_table', help='Save the result table.',default=False,action='store_true') 
    parser.add_argument('--save_models', dest='save_models', help='Save the result table.',default=False,action='store_true') # todo
    parser.add_argument('--resume', dest='resume', help='Resume.',default=False,action='store_true') #todo


    return parser.parse_args()
        
def execute(idx, data_dict, info, method, criteria, model_kmer, prior_params, out_paths, save_models, save_table,locks):
    """
    """
    data = io.load_data(idx,data_dict,info['condition_names'],info['run_names'],min_count=criteria['read_count_min'],max_count=criteria['read_count_max']) #,pooling=False
    models = dict()
    for key,data_at_pos in data.items():
        idx, pos, kmer = key
        model_kmer_mean = model_kmer.loc[kmer,'model_mean']
        model_kmer_tau = 1./(model_kmer.loc[kmer,'model_stdv']**2)
        y_mean = data_at_pos['y'].mean()
        y_tau = 1./(data_at_pos['y'].std()**2)

        K = 2
        group_names = sorted(info['run_names'])

        ### Set up priors.
        priors = {'mu_tau':defaultdict(list),'w':dict()}

        for k in range(K):
            priors['mu_tau']['location'] += [model_kmer_mean]
            priors['mu_tau']['lambda'] += [prior_params['lambda'][k]]
            priors['mu_tau']['alpha'] += [model_kmer_tau]
            priors['mu_tau']['beta'] += [prior_params['beta_scale'][k]*1./model_kmer_tau]
        
        for k,v in priors['mu_tau'].items():
            priors['mu_tau'][k] = numpy.array(v)
        
        priors['w']['concentration'] = numpy.ones([len(data_at_pos['run_names']),K])*1. #GK
        priors['w']['concentration'][:,0] = float(prior_params['concentration'][0])
        priors['w']['concentration'][:,1] = float(prior_params['concentration'][1])
        ###

        ### Fit models.
        models[key] = GMM(method,data_at_pos,priors=priors).fit()
        
    if save_models: #todo: 
        io.save_models(models,out_paths['model_filepath'])
    if save_table:
        # Generating the result table.
        table = io.generate_result_table(models,info['cond2run_dict'])
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
    save_table = args.save_table
    resume = args.resume

    config = Configurator(config_filepath) 
    paths = config.get_paths()
    info = config.get_info()
    method = config.get_method()
    criteria = config.get_criteria()
    prior_params = config.get_priors()

    model_kmer = pandas.read_csv(paths['model_kmer']).set_index('model_kmer')
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
    if save_table:
        if resume:
            gene_ids_done = [line for line in open(out_paths['log'],'r')]                
        else:
            with open(out_paths['table'],'w') as f:
                csv.writer(f,delimiter=',').writerow(io.get_result_table_header(info['cond2run_dict']))
            with open(out_paths['log'],'w') as f:
                f.write(helper.decor_message('diffmod'))


    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=execute,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    ### Load tasks in to task_queue. ###
    # Read index files
    f_index = dict()
    for run_name in info['run_names']:
        df_index = pandas.read_csv(os.path.join(paths['data_dir'],run_name,'dataprep','data.index'),sep=',') # todo
        f_index[run_name] = dict(zip(df_index['gene_id'],zip(df_index['start'],df_index['end'])))
    # Open data files
    f_data = dict()
    for run_name in info['run_names']:
        f_data[run_name] = open(os.path.join(paths['data_dir'],run_name,'dataprep','data.json'),'r') # todo
        
    # Load tasks into task_queue.
    # gene_ids = helper.get_gene_ids(config.filepath)
#    gene_ids = ['ENSG00000168496','ENSG00000204388','ENSG00000123989','ENSG00000170144'] #test data; todo
    # gene_ids = ['ENSG00000159111']
    gene_ids = helper.get_gene_ids(f_index,info)
    print(len(gene_ids),'genes to be testing ...')
    
    for idx in gene_ids:
        if resume and (idx in gene_ids_done):
            continue
            
        data_dict = dict()
        for run_name in info['run_names']:
            try:
                pos_start,pos_end = f_index[run_name][idx]
            except KeyError:
                data_dict[run_name] = None
            else:
                f_data[run_name].seek(pos_start,0)
                json_str = f_data[run_name].read(pos_end-pos_start)
                json_str = '{%s}' %json_str
                data_dict[run_name] = json.loads(json_str) # A data dict for each gene.
                
        # tmp
        out_paths['model_filepath'] = os.path.join(paths['out_dir'],'models','%s.hdf5' %idx)
        #
        # if data_dict[run_name][idx] is not None: # todo: remove this line. Fix in dataprep
        task_queue.put((idx, data_dict, info, method, criteria, model_kmer, prior_params, out_paths,save_models,save_table)) # Blocked if necessary until a free slot is available.

        
    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

    # Close data files
    for run_name in info['run_names']:
        f_data[run_name].close()   
        
    if save_table:
        with open(out_paths['log'],'a+') as f:
            f.write(helper.decor_message('successfully finished'))
        

if __name__ == '__main__':
    """
    Usage:
        xpore-diffmod --config --n_processes --save_table
    """
    main()
