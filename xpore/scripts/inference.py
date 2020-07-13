import argparse
import numpy
import pandas
import os
import multiprocessing 
import h5py
import csv
import json
import subprocess
import pysam #0-based leftmost coordinate
from itertools import product
from pyensembl import EnsemblRelease
from pyensembl import Genome
from tqdm import tqdm
from operator import itemgetter
from collections import defaultdict

from . import helper
from ..utils import misc

def get_args():
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('--mapping', dest='mapping', default=None, help='gene-transcript mapping directory.')
    parser.add_argument('--out_dir', dest='out_dir', help='Output directory.',required=True)


    # Optional
    # parser.add_argument('--features', dest='features', help='Signal features to extract.',type=list,default=['norm_mean'])
    parser.add_argument('--n_processes', dest='n_processes', help='Number of processes.',type=int, default=1)


    
    return parser.parse_args()

def combine(read_name,eventalign_per_read,out_paths,locks):
    eventalign_result = pandas.DataFrame.from_records(eventalign_per_read)

    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    eventalign_result = eventalign_result[cond_successfully_eventaligned]

    keys = ['read_index','contig','position','reference_kmer'] # for groupby
    eventalign_result['length'] = pandas.to_numeric(eventalign_result['end_idx'])-pandas.to_numeric(eventalign_result['start_idx'])
    eventalign_result['sum_norm_mean'] = pandas.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']
    eventalign_result['sum_norm_std'] = pandas.to_numeric(eventalign_result['event_stdv']) * eventalign_result['length']
    eventalign_result['sum_dwell_time'] = pandas.to_numeric(eventalign_result['event_length']) * eventalign_result['length']

    eventalign_result = eventalign_result.groupby(keys)  
    sum_norm_mean = eventalign_result['sum_norm_mean'].sum() 
    sum_norm_std = eventalign_result["sum_norm_std"].sum()
    sum_dwell_time = eventalign_result["sum_dwell_time"].sum()

    start_idx = eventalign_result['start_idx'].min().astype('i8')
    end_idx = eventalign_result['end_idx'].max().astype('i8')
    total_length = eventalign_result['length'].sum()

    eventalign_result = pandas.concat([start_idx,end_idx],axis=1)
    eventalign_result['norm_mean'] = sum_norm_mean / total_length
    eventalign_result["norm_std"] = sum_norm_std / total_length
    eventalign_result["dwell_time"] = sum_dwell_time / total_length
    eventalign_result.reset_index(inplace=True)

    eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]
    eventalign_result['transcriptomic_position'] = pandas.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
    # eventalign_result = misc.str_encode(eventalign_result)
    eventalign_result['read_id'] = [read_name]*len(eventalign_result)

    features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','norm_std','dwell_time','start_idx','end_idx']
    # features_dtype = numpy.dtype([('read_id', 'object'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'),
    #                               ('end_idx', '<i8')])
    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean'] #original features that Ploy's using.

    df_events_per_read = eventalign_result[features]
    # print(df_events_per_read.head())
    
    # write to file.
    df_events_per_read = df_events_per_read.set_index(['transcript_id','read_id'])

    with locks['hdf5'], h5py.File(out_paths['hdf5'],'a') as hf:
        for tx_id,read_id in df_events_per_read.index.unique():
            df2write = df_events_per_read.loc[[tx_id,read_id],:].reset_index() 
            events = numpy.rec.fromrecords(misc.str_encode(df2write[features]),names=features) #,dtype=features_dtype
            
            hf_tx = hf.require_group('%s/%s' %(tx_id,read_id))
            if 'events' in hf_tx:
                continue
            else:
                hf_tx['events'] = events
    
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write('%s\n' %(read_name))    


def prepare_for_inference(tx,gt_dir,read_task,all_kmers,out_dir,locks):
    reads = numpy.concatenate(read_task)
    try:
        gt_map = pandas.read_csv(os.path.join(gt_dir,tx,"gt_mapping.csv.gz")).set_index("tx_pos")
    except FileNotFoundError:
        with locks['log'], open(os.path.join(out_dir, "prepare_for_inference.log"),'a') as f:
            f.write('Error at %s\n' %(tx))
        return

    gene = gt_map.iloc[0]["g_id"]
    reads = reads[numpy.isin(reads["reference_kmer"], all_kmers)] # Filter for motifs
    reads = reads[numpy.argsort(reads["transcriptomic_position"])] # sort by reference kmer
    positions, indices = numpy.unique(reads["transcriptomic_position"],return_index=True) # retrieve group indexing

    for i in range(len(positions)):
        pos = positions[i]
        gpos = gt_map.loc[pos]["g_pos"]
        kmer = gt_map.loc[pos]["kmer"]
        if len(positions) > 1:
            start_idx = indices[i]
            end_idx = indices[i + 1] if i < len(positions) - 1 else None
            read = reads[start_idx:end_idx]
            
            # Converting to numpy array
            
            X = read[["norm_mean", "norm_std", "dwell_time"]]\
                        .astype([('norm_mean', '<f8'), ('norm_std', '<f8'), ('dwell_time', '<f8')]).view('<f8')
            read_ids = read["read_id"].view('<S36')
            start_event_indices = read["start_idx"].view('<i8')
            end_event_indices = read["end_idx"].view('<i8')
        else:
            read = reads[0]

            # Converting to numpy array when there is only one entry
            
            X = numpy.array(read[["norm_mean", "norm_std", "dwell_time"]].tolist())
            read_ids = numpy.array(read["read_id"].tolist())
            start_event_indices = numpy.array(read["start_idx"].tolist())
            end_event_indices = numpy.array(read["end_idx"].tolist())
        
        # Reshaping columns
        X = X.reshape(-1, 3)
        read_ids = read_ids.reshape(-1, 1)
        start_event_indices = start_event_indices.reshape(-1, 1)
        end_event_indices = end_event_indices.reshape(-1, 1)

        # Saving output in hdf5 file format

        n_reads = len(X)
        fname = os.path.join(out_dir, '{}_{}_{}_{}_{}_{}.hdf5'.format(gene, tx, gpos, pos, kmer, n_reads))
        with h5py.File(fname, 'w') as f:
            assert(n_reads == len(read_ids))
            f['X'] = X
            f['read_ids'] = read_ids
            f['start_idx'] = start_event_indices
            f['end_idx'] = end_event_indices
        f.close()

    with locks['log'], open(os.path.join(out_dir, "prepare_for_inference.log"),'a') as f:
        f.write('%s\n' %(tx))    

        
def prepare_for_inference_tx(tx, read_task, all_kmers, out_dir, locks):    
    reads = numpy.concatenate(read_task)
    reads = reads[numpy.isin(reads["reference_kmer"], all_kmers)] # Filter for motifs
    reads = reads[numpy.argsort(reads["transcriptomic_position"])] # sort by reference kmer
    positions, indices = numpy.unique(reads["transcriptomic_position"],return_index=True) # retrieve group indexing

    for i in range(len(positions)):
        pos = positions[i]
        kmer = gt_map.loc[pos]["kmer"]
        if len(positions) > 1:
            start_idx = indices[i]
            end_idx = indices[i + 1] if i < len(positions) - 1 else None
            read = reads[start_idx:end_idx]
            
            # Converting to numpy array
            
            X = read[["norm_mean", "norm_std", "dwell_time"]]\
                        .astype([('norm_mean', '<f8'), ('norm_std', '<f8'), ('dwell_time', '<f8')]).view('<f8')
            read_ids = read["read_id"].view('<S36')
            start_event_indices = read["start_idx"].view('<i8')
            end_event_indices = read["end_idx"].view('<i8')
        else:
            read = reads[0]

            # Converting to numpy array when there is only one entry
            
            X = numpy.array(read[["norm_mean", "norm_std", "dwell_time"]].tolist())
            read_ids = numpy.array(read["read_id"].tolist())
            start_event_indices = numpy.array(read["start_idx"].tolist())
            end_event_indices = numpy.array(read["end_idx"].tolist())
        
        # Reshaping columns
        X = X.reshape(-1, 3)
        read_ids = read_ids.reshape(-1, 1)
        start_event_indices = start_event_indices.reshape(-1, 1)
        end_event_indices = end_event_indices.reshape(-1, 1)

        # Saving output in hdf5 file format

        n_reads = len(X)
        fname = os.path.join(out_dir, '{}_{}_{}_{}.hdf5'.format(tx, pos, kmer, n_reads))
        with h5py.File(fname, 'w') as f:
            assert(n_reads == len(read_ids))
            f['X'] = X
            f['read_ids'] = read_ids
            f['start_idx'] = start_event_indices
            f['end_idx'] = end_event_indices
        f.close()

    with locks['log'], open(os.path.join(out_dir, "prepare_for_inference.log"),'a') as f:
        f.write('%s\n' %(tx))   
        
        
def parallel_prepare_for_inference(eventalign_filepath,gt_dir,eventalign_prep_dir,n_processes):
    # Create output path and locks.
    out_dir = os.path.join(eventalign_prep_dir, "inference")
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    locks = {'log': multiprocessing.Lock()}
    log_path = os.path.join(out_dir, "prepare_for_inference.log")
    # Create empty files for logs.
    open(log_path,'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)
    
    task_func = prepare_for_inference if gt_dir is not None else prepare_for_inference_tx
    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=task_func,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()
        
    ## Load tasks into task_queue. A task is a read information from a specific site.

    # Only include reads that conform to DRACH motifs

    all_kmers = numpy.array(["".join(x) for x in product(['A', 'G', 'T'], ['G', 'A'], ['A'], ['C'], ['A', 'C', 'T'])], dtype='S5')
    with h5py.File(eventalign_filepath, 'r') as f:
        for tx in f:
            print("Preprocessing transcript {}".format(tx))
            read_task = []
            for read in f[tx]:
                read_task.append(f[tx][read]['events'][:])
            if gt_dir is not None:
                task_queue.put((tx,gt_dir,read_task,all_kmers,out_dir))
            else:
                task_queue.put((tx, read_task, all_kmers, out_dir))

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(log_path,'a+') as f:
        f.write(helper.decor_message('successfully finished'))

        
def main():
    args = get_args()
    #
    n_processes = args.n_processes        
    out_dir = args.out_dir
    gt_mapping_dir = args.mapping

    
    parallel_prepare_for_inference(os.path.join(out_dir, 'eventalign.hdf5'),gt_mapping_dir,out_dir,n_processes)

if __name__ == '__main__':
    main()


