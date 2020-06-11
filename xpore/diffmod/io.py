import numpy as np
import h5py
import os
from collections import OrderedDict, defaultdict
import itertools
import scipy.stats

from ..utils import stats
from .gmm import GMM


def get_dummies(x):
    X = []
    labels = sorted(list(set(x)))
    for label in labels:
        X += [x == label]
    # labels = [label.encode('UTF-8') for label in labels]
    return np.array(X).T, labels


def load_data(idx, data_dict, min_count, max_count, pooling=False): 
    """
    Parameters
    ----------
    data_dict: dict of ...
        Data for a gene 
    """
    
    # Create all (pos,kmer) pairs from all runs.
    position_kmer_pairs = []
    for (condition_name,run_name), d_dict in data_dict.items(): # data_dict[run_name][idx][position][kmer]
        pairs = []
        if d_dict is not None:
            for pos in d_dict[idx].keys():
                for kmer in d_dict[idx][pos].keys():
                    pairs += [(pos,kmer)]                
            position_kmer_pairs += [pairs]

    position_kmer_pairs = set(position_kmer_pairs[0]).intersection(*position_kmer_pairs)

    data = OrderedDict()
    for pos,kmer in position_kmer_pairs:  
        y, read_ids, condition_labels, run_labels = [], [], [], []
        n_reads = defaultdict(list)
        
        for (condition_name,run_name), d_dict in data_dict.items():
            if d_dict is not None:
                norm_means = d_dict[idx][pos][kmer]#['norm_means']
                n_reads_per_run = len(norm_means)
                # In case of pooling==False, if not enough reads, don't include them. 
                if (not pooling) and ((n_reads_per_run < min_count) or (n_reads_per_run > max_count)):
                    continue
                #
                n_reads[condition_name] += [n_reads_per_run]
                y += norm_means
                # read_ids += list(data_dict[run_name][idx][pos][kmer]['read_ids'][:])
                condition_labels += [condition_name]*n_reads_per_run
                run_labels += [run_name]*n_reads_per_run

        y = np.array(y)
        # read_ids = np.array(read_ids)
        condition_labels = np.array(condition_labels)
        run_labels = np.array(run_labels)
        
        # Filter those sites that don't have enough reads.
        if len(y) == 0:  # no reads at all.
            continue
        conditions_incl = []
        unique_condition_names = {condition_name for condition_name,_ in data_dict.keys()}
        if pooling: # At the modelling step all the reads from the same condition will be combined.
            for condition_name in unique_condition_names:
                if (sum(n_reads[condition_name]) >= min_count) and (sum(n_reads[condition_name]) <= max_count):
                    conditions_incl += [condition_name]
        else:
            for condition_name in unique_condition_names:
                if (np.array(n_reads[condition_name]) >= min_count).any() and (np.array(n_reads[condition_name]) <= max_count).any():
                    conditions_incl += [condition_name]
                    
        if len(conditions_incl) < 2:
            continue

        # Get dummies
        x, condition_names_dummies = get_dummies(condition_labels)
        r, run_names_dummies = get_dummies(run_labels)

        key = (idx, pos, kmer)

        data[key] = {'y': y, 'x': x, 'r': r, 'condition_names': condition_names_dummies, 'run_names': run_names_dummies, 'y_condition_names': condition_labels, 'y_run_names': run_labels} # , 'read_ids': read_ids

    return data


def save_result_table(table, out_filepath):
    out_file = h5py.File(out_filepath, 'w')
    out_file['result'] = table  # Structured np array.
    out_file.close()


def save_models(models, model_filepath):  # per gene/transcript
    """
    Save model parameters.

    Parameters
    ----------
    models
        Learned models.
    model_filepath: str
        Path to save the models.
    """
    model_file = h5py.File(model_filepath, 'w')
    for model_key, model in models.items():  
        idx, position, kmer = model_key

        position = str(position)
        if idx not in model_file:
            model_file.create_group(idx)
        model_file[idx].create_group(position)
        model_file[idx][position].attrs['kmer'] = kmer.encode('UTF-8')
        model_file[idx][position].create_group('info')
        for key, value in model.info.items():
            model_file[idx][position]['info'][key] = value

        model_file[idx][position].create_group('nodes')  # ['x','y','z','w','mu_tau'] => store only their params
        for node_name in model.nodes:
            model_file[idx][position]['nodes'].create_group(node_name)
            if model.nodes[node_name].data is not None: # To make it optional.
                value = model.nodes[node_name].data
                if node_name in ['y_run_names','y_condition_names']:
                    value = [val.encode('UTF-8') for val in value]
                model_file[idx][position]['nodes'][node_name]['data'] = value
            if model.nodes[node_name].params is None:
                continue
            for param_name, value in model.nodes[node_name].params.items():
                if param_name == 'group_names':
                    value = [val.encode('UTF-8') for val in value]
                model_file[idx][position]['nodes'][node_name][param_name] = value

    model_file.close()


def load_models(model_filepath):  # per gene/transcript #Todo: refine.
    """
    Construct a model and load model parameters.

    Parameters
    ----------
    model_filepath: str
        Path where the model is stored.

    Return
    ------
    models
        Models for each genomic position.
    """

    model_file = h5py.File(model_filepath, 'r')
    models = {}
    data = defaultdict(dict)
    for idx in model_file:
        for position in model_file[idx]:
            inits = {'info': None, 'nodes': {'x': {}, 'y': {}, 'w': {}, 'mu_tau': {}, 'z': {}}}
            kmer = model_file[idx][position].attrs['kmer']
            key = (idx, position, kmer)
            # for k in model_file[idx][position]['info']:
            #     inits['info'] = model_file[idx][position]['info'][k]
            for node_name, params in model_file[idx][position]['nodes'].items():
                for param_name, value in params.items():
                    if param_name == 'data':
                        data[key][node_name] = value[:]
                    else:
                        inits['nodes'][node_name][param_name] = value[:]
                # for param_name, value in priors.items():
                #     inits['nodes'][node_name][param_name] = value[:]

            models[key] = GMM(data,inits=inits)

    model_file.close()

    return models,data  # {(idx,position,kmer): GMM obj}

def get_ordered_condition_run_names(data_info):
    unique_condition_names = list(data_info.keys())
    unique_condition_names_sorted = sorted(unique_condition_names)
    unique_run_names_sorted = sum([list(run_names.keys()) for _, run_names in data_info.items()],[])

    return unique_condition_names_sorted, unique_run_names_sorted

def get_result_table_header(data_info,method):
    condition_names,run_names = get_ordered_condition_run_names(data_info)
    ### stats header
    stats_pairwise = []
    for cond1, cond2 in itertools.combinations(condition_names, 2):
        pair = '_vs_'.join((cond1, cond2)) # cond1 - cond2
        stats_pairwise += ['diff_mod_rate_%s' % pair, 'pval_%s' % pair, 'z_score_%s' % pair]
    stats_one_vs_all = []
    for condition_name in condition_names:
        pair = '%s_vs_all' %condition_name # condition_name - the rest
        stats_one_vs_all += ['diff_mod_rate_%s' % pair, 'pval_%s' %pair, 'z_score_%s' %pair]

    ### position header
    header = ['id', 'position', 'kmer']
    
    # header += ['p_overlap']
    # header += ['x_x1', 'y_x1', 'x_x2', 'y_x2']
    
    if method['pooling']:
        names = condition_names
    else:
        names = run_names
    
    ### stats header
    header += stats_pairwise
    if len(condition_names) > 2:
        header += stats_one_vs_all
    ###
    ### model header
    for name in names:
        header += ['mod_rate_%s' % name]
    for name in names:
        header += ['coverage_%s' % name]
    header += ['mu_unmod', 'mu_mod', 'sigma2_unmod', 'sigma2_mod', 'conf_mu_unmod', 'conf_mu_mod','mod_assignment']
    ###
    
    ### prefiltering header
    if method['prefiltering']:
        header += [method['prefiltering']['method']]
    return header

def generate_result_table(models, data_info):  # per idx (gene/transcript)
    """
    Generate a table containing learned model parameters and statistic tests.

    Parameters
    ----------
    models
        Learned models for individual genomic positions of a gene.
    group_labels
        Labels of samples.
    data_inf
        Dict

    Returns
    -------
    table
        List of tuples.
    """

    ###
    condition_names,run_names = get_ordered_condition_run_names(data_info) # information from the config file used for modelling.
    ###

    ###
    table = []
    for key, (model,prefiltering) in models.items():
        idx, position, kmer = key
        mu = model.nodes['mu_tau'].expected()  # K
        sigma2 = 1./model.nodes['mu_tau'].expected(var='gamma')  # K
        var_mu = model.nodes['mu_tau'].variance(var='normal')  # K
        # mu = model.nodes['y'].params['mean']
        # sigma2 = model.nodes['y'].params['variance']
        w = model.nodes['w'].expected()  # GK
        N = model.nodes['y'].params['N'].round()  # GK
        N0 = N[:, 0].squeeze()
        N1 = N[:, 1].squeeze()
        w0 = w[:, 0].squeeze()
        coverage = np.sum(model.nodes['y'].params['N'], axis=-1)  # GK => G # n_reads per group

        p_overlap, list_cdf_at_intersections = stats.calc_prob_overlapping(mu, sigma2)

        model_group_names = model.nodes['x'].params['group_names'] #condition_names if pooling, run_names otherwise.
        
        ### Cluster assignment ###
        conf_mu = [calculate_confidence_cluster_assignment(mu[0],model.kmer_signal),calculate_confidence_cluster_assignment(mu[1],model.kmer_signal)]
    
        cluster_idx = {}
        if conf_mu[0] > conf_mu[1]:
            cluster_idx['unmod'] = 0
            cluster_idx['mod'] = 1
        else:
            cluster_idx['unmod'] = 1
            cluster_idx['mod'] = 0

        mu_assigned = [mu[cluster_idx['unmod']],mu[cluster_idx['mod']]]
        sigma2_assigned = [sigma2[cluster_idx['unmod']],sigma2[cluster_idx['mod']]] 
        conf_mu = [conf_mu[cluster_idx['unmod']],conf_mu[cluster_idx['mod']]]
        w_mod = w[:,cluster_idx['mod']]
        mod_assignment = [['higher','lower'][(mu[0]<mu[1])^cluster_idx['mod']]]
            
        
        ### calculate stats_pairwise
        stats_pairwise = []
        for cond1, cond2 in itertools.combinations(condition_names, 2):
            if model.method['pooling']:
                cond1, cond2 = [cond1], [cond2]
            else:
                cond1, cond2 = list(data_info[cond1].keys()), list(data_info[cond2].keys())
            if any(r in model_group_names for r in cond1) and any(r in model_group_names for r in cond2):
                w_cond1 = w[np.isin(model_group_names, cond1), cluster_idx['mod']].flatten()
                w_cond2 = w[np.isin(model_group_names, cond2), cluster_idx['mod']].flatten()
                n_cond1 = coverage[np.isin(model_group_names, cond1)]
                n_cond2 = coverage[np.isin(model_group_names, cond2)]

                z_score, p_ws = stats.z_test(w_cond1, w_cond2, n_cond1, n_cond2) # two=tailed
                w_mod_mean_diff = np.mean(w_cond1)-np.mean(w_cond2)

                stats_pairwise += [w_mod_mean_diff, p_ws, z_score]
            else:
                stats_pairwise += [None, None, None]

        if len(condition_names) > 2:
            ### calculate stats_one_vs_all
            stats_one_vs_all = []
            for cond in condition_names:
                if model.method['pooling']:
                    cond = [cond]
                else:
                    cond = list(data_info[cond].keys())
                if any(r in model_group_names for r in cond):
                    w_cond1 = w[np.isin(model_group_names, cond), cluster_idx['mod']].flatten()
                    w_cond2 = w[~np.isin(model_group_names, cond), cluster_idx['mod']].flatten()
                    n_cond1 = coverage[np.isin(model_group_names, cond)]
                    n_cond2 = coverage[~np.isin(model_group_names, cond)]

                    z_score, p_ws = stats.z_test(w_cond1, w_cond2, n_cond1, n_cond2)
                    w_mod_mean_diff = np.mean(w_cond1)-np.mean(w_cond2)

                    stats_one_vs_all += [w_mod_mean_diff, p_ws, z_score]
                else:
                    stats_one_vs_all += [None, None, None]

        ###
        w_mod_ordered, coverage_ordered = [], [] # ordered by conditon_names or run_names based on headers.        
        if model.method['pooling']:
            names = condition_names
        else:
            names = run_names
        for name in names:
            if name in model_group_names:
                w_mod_ordered += list(w_mod[np.isin(model_group_names, name)])
                coverage_ordered += list(coverage[np.isin(model_group_names, name)])
            else:
                w_mod_ordered += [None]
                coverage_ordered += [None]

        ###
        ### prepare values to write
        row = [idx, position, kmer]
        row += stats_pairwise
        if len(condition_names) > 2:
            row += stats_one_vs_all

        # row += [p_overlap]
        # row += list_cdf_at_intersections
        row += list(w_mod_ordered)
        row += list(coverage_ordered)
        row += mu_assigned + sigma2_assigned + conf_mu + mod_assignment


        if prefiltering is not None:
            row += [prefiltering[model.method['prefiltering']['method']]]
        ### Filtering those positions with a nearly single distribution.
        cdf_threshold = 0.1
        x_x1, y_x1, x_x2, y_x2 = list_cdf_at_intersections
        is_not_inside = ((y_x1 < cdf_threshold) & (x_x1 < cdf_threshold)) | ((y_x2 < cdf_threshold) & (x_x2 < cdf_threshold)) | (( (1-y_x1) < cdf_threshold) & ((1-x_x1) < cdf_threshold)) | (( (1-y_x2) < cdf_threshold) & ((1-x_x2) < cdf_threshold))
        if (p_overlap <= 0.5) and (is_not_inside):
            table += [tuple(row)]

    return table


def calculate_confidence_cluster_assignment(mu,kmer_signal):
    cdf = scipy.stats.norm.cdf(kmer_signal['mean'] - abs(kmer_signal['mean']-mu), loc=kmer_signal['mean'], scale=kmer_signal['std'])
    return cdf*2
