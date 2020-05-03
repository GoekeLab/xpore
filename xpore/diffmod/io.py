import numpy
import h5py
import os
from tqdm import tqdm  # progress bar.
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
    return numpy.array(X).T, labels


def load_data(idx, data_dict, condition_names, run_names, min_count=30, max_count=3000, pooling=False): 
    """
    Parameters
    ----------
    data_dict: dict of ...
        Data for a gene 
    """
    
    # Create all (pos,kmer) pairs from all runs.
    position_kmer_pairs = []
    for run_name in run_names: # data_dict[run_name][idx][position][kmer]
        pairs = []
        if data_dict[run_name] is not None:
            for pos in data_dict[run_name][idx].keys():
                for kmer in data_dict[run_name][idx][pos].keys():
                    pairs += [(pos,kmer)]                
            position_kmer_pairs += [pairs]

    position_kmer_pairs = set(position_kmer_pairs[0]).intersection(*position_kmer_pairs)

    data = OrderedDict()
    for pos,kmer in position_kmer_pairs:  
        y, read_ids, condition_labels, run_labels = [], [], [], []
        n_reads = defaultdict(list)
        for condition_name, run_name in zip(condition_names, run_names):
            if data_dict[run_name] is not None:
                norm_means = data_dict[run_name][idx][pos][kmer]['norm_means']
                n_reads_per_run = len(norm_means)
                # In case of pooling==False, if not enough reads, don't include them. 
                if (not pooling) and ((n_reads_per_run < min_count) or (n_reads_per_run > max_count)):
                    continue
                #
                n_reads[condition_name] += [n_reads_per_run]
                y += norm_means
                read_ids += list(data_dict[run_name][idx][pos][kmer]['read_ids'][:])
                condition_labels += [condition_name]*n_reads_per_run
                run_labels += [run_name]*n_reads_per_run

        y = numpy.array(y)
        read_ids = numpy.array(read_ids)
        condition_labels = numpy.array(condition_labels)
        run_labels = numpy.array(run_labels)
        
        # Filter those sites that don't have enough reads.
        if len(y) == 0:  # no reads at all.
            continue
        conditions_incl = []
        if pooling: # At the modelling step all the reads from the same condition will be combined.
            for condition_name in set(condition_names):
                if (sum(n_reads[condition_name]) >= min_count) and (sum(n_reads[condition_name]) <= max_count):
                    conditions_incl += [condition_name]
        else:
            for condition_name in set(condition_names):
                if (numpy.array(n_reads[condition_name]) >= min_count).any() and (numpy.array(n_reads[condition_name]) <= max_count).any():
                    conditions_incl += [condition_name]
                    
        if len(conditions_incl) < 2:
            continue

        # Get dummies
        x, condition_names_dummies = get_dummies(condition_labels)
        r, run_names_dummies = get_dummies(run_labels)

        key = (idx, pos, kmer)

        data[key] = {'y': y, 'x': x, 'r': r, 'condition_names': condition_names_dummies, 'run_names': run_names_dummies, 'read_ids': read_ids, 'y_condition_names': condition_labels, 'y_run_names': run_labels}

    return data


def save_result_table(table, out_filepath):
    out_file = h5py.File(out_filepath, 'w')
    out_file['result'] = table  # Structured numpy array.
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
    for model_key, model in models.items():  # tqdm(models.items()):
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
            if model.nodes[node_name].data is not None: # To be optional.
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
        for position in tqdm(model_file[idx]):
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

def get_result_table_header(cond2run_dict,pooling=False):
    condition_names,run_names = get_ordered_condition_run_names(cond2run_dict)
    ### stats headers
    stats_pairwise = []
    for cond1, cond2 in itertools.combinations(condition_names, 2):
        pair = '_vs_'.join((cond1, cond2))
        stats_pairwise += ['p_ws_%s' % pair, 'ws_mean_diff_%s' % pair, 'abs_z_score_%s' % pair]
    stats_one_vs_all = []
    for condition_name in condition_names:
        stats_one_vs_all += ['p_ws_%s_vs_all' % condition_name, 'ws_mean_diff_%s_vs_all' % condition_name, 'abs_z_score_%s_vs_all' % condition_name]

    header = ['idx', 'position', 'kmer']
    header += ['p_overlap']
    header += ['x_x1', 'y_x1', 'x_x2', 'y_x2']
    
    if pooling:
        names = condition_names
    else:
        names = run_names
    # for name in names:
    #     header += ['w_min_%s' % name]
    for name in names:
        header += ['coverage_%s' % name]
    header += ['mu_unmod', 'mu_mod', 'sigma2_unmod', 'sigma2_mod', 'conf_mu_unmod', 'conf_mu_mod']
    for name in names:
        header += ['w_mod_%s' % name]
    ###
    header += stats_pairwise
    if len(condition_names) > 2:
        header += stats_one_vs_all
    ###

    return header

def get_ordered_condition_run_names(cond2run_dict):
    condition_names = sorted(list(set(cond2run_dict.keys())))
    run_names = sorted(list(set(sum(list(cond2run_dict.values()), []))))
    return condition_names,run_names

# def generate_result_table(models, cond2run_dict):  # per idx (gene/transcript)


def generate_result_table(models, cond2run_dict):  # per idx (gene/transcript)
    """
    Generate a table containing learned model parameters and statistic tests.

    Parameters
    ----------
    models
        Learned models for individual genomic positions of a gene.
    group_labels
        Labels of samples.
    cond2run_dict
        Dict mapping condition_names to list of run_names

    Returns
    -------
    table
        List of tuples.
    """

    ###
    condition_names,run_names = get_ordered_condition_run_names(cond2run_dict) # information from the config file used for modelling.
    ###

    ###
    table = []
    for key, model in models.items():
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
        coverage = numpy.sum(model.nodes['y'].params['N'], axis=-1)  # GK => G # n_reads per group

        p_overlap, list_cdf_at_intersections = stats.calc_prob_overlapping(mu, sigma2)

        model_group_names = model.nodes['x'].params['group_names'] #condition_names if pooling, run_names otherwise.

        ### calculate stats_pairwise
        stats_pairwise = []
        for cond1, cond2 in itertools.combinations(condition_names, 2):
            if model.method['pooling']:
                cond1, cond2 = [cond1], [cond2]
            else:
                cond1, cond2 = cond2run_dict[cond1], cond2run_dict[cond2]
            if any(r in model_group_names for r in cond1) and any(r in model_group_names for r in cond2):
                w_cond1 = w[numpy.isin(model_group_names, cond1), 0].flatten()
                w_cond2 = w[numpy.isin(model_group_names, cond2), 0].flatten()
                n_cond1 = coverage[numpy.isin(model_group_names, cond1)]
                n_cond2 = coverage[numpy.isin(model_group_names, cond2)]

                z_score, p_ws = stats.z_test(w_cond1, w_cond2, n_cond1, n_cond2)
                ws_mean_diff = abs(numpy.mean(w_cond1)-numpy.mean(w_cond2))
                abs_z_score = abs(z_score)

                stats_pairwise += [p_ws, ws_mean_diff, abs_z_score]
            else:
                stats_pairwise += [None, None, None]

        if len(condition_names) > 2:
            ### calculate stats_one_vs_all
            stats_one_vs_all = []
            for cond in condition_names:
                if model.method['pooling']:
                    cond = [cond]
                else:
                    cond = cond2run_dict[cond]
                if any(r in model_group_names for r in cond):
                    w_cond1 = w[numpy.isin(model_group_names, cond), 0].flatten()
                    w_cond2 = w[~numpy.isin(model_group_names, cond), 0].flatten()
                    n_cond1 = coverage[numpy.isin(model_group_names, cond)]
                    n_cond2 = coverage[~numpy.isin(model_group_names, cond)]

                    z_score, p_ws = stats.z_test(w_cond1, w_cond2, n_cond1, n_cond2)
                    ws_mean_diff = abs(numpy.mean(w_cond1)-numpy.mean(w_cond2))
                    abs_z_score = abs(z_score)

                    stats_one_vs_all += [p_ws, ws_mean_diff, abs_z_score]
                else:
                    stats_one_vs_all += [None, None, None]

        ### lower, higher clusters
        w_min = w0
        if mu[1] < mu[0]:
            mu = mu[::-1]
            sigma2 = sigma2[::-1]
            w_min = 1-w0
        w_max = 1-w_min
        ###
        w_min_ordered, w_max_ordered, coverage_ordered = [], [], [] # ordered by conditon_names / run_names based on headers.        
        if model.method['pooling']:
            names = condition_names
        else:
            names = run_names
        for name in names:
            if name in model_group_names:
                w_min_ordered += list(w_min[numpy.isin(model_group_names, name)])
                w_max_ordered += list(w_max[numpy.isin(model_group_names, name)])
                coverage_ordered += list(coverage[numpy.isin(model_group_names, name)])
            else:
                w_min_ordered += [None]
                w_max_ordered += [None]
                coverage_ordered += [None]
        
        ### Cluster assignment ###
        # Calculate confidence of mu_min, mu_max
        # Given that mu[0] < mu[1],
        conf_mu_min,conf_mu_max = calculate_confidence_cluster_assignment(mu[0],model.kmer_signal),calculate_confidence_cluster_assignment(mu[1],model.kmer_signal)
        if conf_mu_min > conf_mu_max:
            mu_assigned = [mu[0],mu[1]] 
            sigma2_assigned = [sigma2[0],sigma2[1]] 
            w_mod_ordered = w_max_ordered
            conf_mu = [conf_mu_min, conf_mu_max]
        else:
            mu_assigned = [mu[1],mu[0]] 
            sigma2_assigned = [sigma2[1],sigma2[0]] 
            w_mod_ordered = w_min_ordered
            conf_mu = [conf_mu_max, conf_mu_min]
        #
        
        ###
        ### prepare values to write
        row = [idx, position, kmer] + [p_overlap]
        row += list_cdf_at_intersections
        row += list(coverage_ordered)
        row += mu_assigned + sigma2_assigned + conf_mu
        row += list(w_mod_ordered)

        row += stats_pairwise
        if len(condition_names) > 2:
            row += stats_one_vs_all

        table += [tuple(row)]

    return table

def calculate_confidence_cluster_assignment(mu,kmer_signal):
    cdf = scipy.stats.norm.cdf(kmer_signal['mean'] - abs(kmer_signal['mean']-mu), loc=kmer_signal['mean'], scale=kmer_signal['std'])
    return cdf*2

