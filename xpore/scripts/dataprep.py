import argparse
import numpy as np
import pandas as pd
import os
import multiprocessing 
import h5py
import csv
import json
from pyensembl import EnsemblRelease
from operator import itemgetter
from collections import defaultdict

from . import helper
from ..utils import misc

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required arguments
    required.add_argument('--eventalign', dest='eventalign', help='eventalign filepath, the output from nanopolish.',required=True)
    required.add_argument('--summary', dest='summary', help='eventalign summary filepath, the output from nanopolish.',required=True)
    required.add_argument('--out_dir', dest='out_dir', help='output directory.',required=True)
    
    required.add_argument('--ensembl', dest='ensembl', help='ensembl version for gene-transcript mapping.',type=int, default=91)
    required.add_argument('--species', dest='species', help='species for ensembl gene-transcript mapping.', default='homo_sapiens')

    # Optional
    # parser.add_argument('--features', dest='features', help='Signal features to extract.',type=list,default=['norm_mean'])
    optional.add_argument('--genome', dest='genome', help='to run on Genomic coordinates. Without this argument, the program will run on transcriptomic coordinates',default=False,action='store_true') 
    optional.add_argument('--n_processes', dest='n_processes', help='number of processes to run.',type=int, default=1)
    optional.add_argument('--readcount_max', dest='readcount_max', help='maximum read counts per gene.',type=int, default=1000)
    optional.add_argument('--resume', dest='resume', help='resume from the previous run.',default=False,action='store_true') #todo

    parser._action_groups.append(optional)
    return parser.parse_args()

def combine(read_name,eventalign_per_read,out_paths,locks):
    eventalign_result = pd.DataFrame.from_records(eventalign_per_read)

    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    eventalign_result = eventalign_result[cond_successfully_eventaligned]

    keys = ['read_index','contig','position','reference_kmer'] # for groupby
    eventalign_result['length'] = pd.to_numeric(eventalign_result['end_idx'])-pd.to_numeric(eventalign_result['start_idx'])
    eventalign_result['sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']

    eventalign_result = eventalign_result.groupby(keys)  
    sum_norm_mean = eventalign_result['sum_norm_mean'].sum() 
    start_idx = eventalign_result['start_idx'].min()
    end_idx = eventalign_result['end_idx'].max()
    total_length = eventalign_result['length'].sum()

    eventalign_result = pd.concat([start_idx,end_idx],axis=1)
    eventalign_result['norm_mean'] = sum_norm_mean/total_length

    eventalign_result.reset_index(inplace=True)

    # eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]
    eventalign_result['transcript_id'] = eventalign_result['contig']
    eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
    # eventalign_result = misc.str_encode(eventalign_result)
    eventalign_result['read_id'] = [read_name]*len(eventalign_result)

    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx']
    # features_dtype = np.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'), ('end_idx', '<i8')])
    features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean']

    df_events_per_read = eventalign_result[features]
    # print(df_events_per_read.head())
    
    # write to file.
    df_events_per_read = df_events_per_read.set_index(['transcript_id','read_id'])

    with locks['hdf5'], h5py.File(out_paths['hdf5'],'a') as hf:
        for tx_id,read_id in df_events_per_read.index.unique():
            df2write = df_events_per_read.loc[[tx_id,read_id],:].reset_index() 
            events = np.rec.fromrecords(misc.str_encode(df2write[features]),names=features) #,dtype=features_dtype

            hf_tx = hf.require_group('%s/%s' %(tx_id,read_id))
            if 'events' in hf_tx:
                continue
            else:
                hf_tx['events'] = events
    
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write('%s\n' %(read_name))    


def parallel_combine(eventalign_filepath,summary_filepath,out_dir,n_processes,resume):
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['hdf5','log']:
        out_paths[out_filetype] = os.path.join(out_dir,'eventalign.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
        
        
    read_names_done = []
    if resume and os.path.exists(out_paths['log']):
        read_names_done = [line.rstrip('\n') for line in open(out_paths['log'],'r')]
    else:
        # Create empty files.
        open(out_paths['hdf5'],'w').close()
        open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=combine,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()
        
    ## Load tasks into task_queue. A task is eventalign information of one read.            
    with helper.EventalignFile(eventalign_filepath) as eventalign_file, open(summary_filepath,'r') as summary_file:
        
        reader_summary = csv.DictReader(summary_file, delimiter="\t")
        reader_eventalign = csv.DictReader(eventalign_file, delimiter="\t")

        row_summary = next(reader_summary)
        read_name = row_summary['read_name']
        read_index = row_summary['read_index']
        eventalign_per_read = []
        for row_eventalign in reader_eventalign:
            if (row_eventalign['read_index'] == read_index):
                eventalign_per_read += [row_eventalign]
            else: 
                # Load a read info to the task queue.
                if read_name not in read_names_done:
                    task_queue.put((read_name,eventalign_per_read,out_paths))
                # Next read.
                try:
                    row_summary = next(reader_summary)
                except StopIteration: # no more read.
                    break
                else:
                    read_index = row_summary['read_index']
                    read_name = row_summary['read_name']
                    assert row_eventalign['read_index'] == read_index 
                    eventalign_per_read = [row_eventalign]
    
    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(out_paths['log'],'a+') as f:
        f.write(helper.decor_message('successfully finished'))

    
def t2g(gene_id,ensembl):
    tx_ids = []
    t2g_dict = {}
    for tx in ensembl.gene_by_id(gene_id).transcripts:
        tx_seq = ensembl.transcript_sequence(tx.id)
        if tx_seq is None:
            continue
        for interval in tx.exon_intervals:
            for g_pos in range(interval[0],interval[1]+1): # Exclude the rims of exons.
                tx_pos = tx.spliced_offset(g_pos)
                if (interval[0] <= g_pos < interval[0]+2) or (interval[1]-2 < g_pos <= interval[1]): # Todo: To improve the mapping
                    kmer = 'XXXXX'
                else:
                    kmer = tx_seq[tx_pos-2:tx_pos+3]
                t2g_dict[(tx.id,tx_pos)] = (tx.contig,gene_id,g_pos,kmer) # tx.contig is chromosome.
        tx_ids += [tx.id]
                
    return tx_ids, t2g_dict
            
def parallel_preprocess_gene(ensembl,out_dir,n_processes,readcount_max,resume):
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['json','index','log','readcount']:
        out_paths[out_filetype] = os.path.join(out_dir,'data.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
                
    # Writing the starting of the files.
    gene_ids_done = []
    if resume and os.path.exists(out_paths['index']):
        df_index = pd.read_csv(out_paths['index'],sep=',')
        gene_ids_done = list(df_index['idx'].unique())
    else:
        # with open(out_paths['json'],'w') as f:
        #     f.write('{\n')
        #     f.write('"genes":{')
        open(out_paths['json'],'w').close()
        with open(out_paths['index'],'w') as f:
            f.write('idx,start,end\n') # header
        with open(out_paths['readcount'],'w') as f:
            f.write('idx,n_reads\n') # header
        open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=preprocess_gene,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Get all gene ids.
    gene_ids = set()
    tx_ensembl = dict()
    with h5py.File(os.path.join(out_dir,'eventalign.hdf5'),'r') as f:
        for tx_id in f.keys():
            tx_id,tx_version = tx_id.split('.') # Based on Ensembl
            tx_ensembl[tx_id] = tx_version
            try:
                g_id = ensembl.transcript_by_id(tx_id).gene_id 
            except ValueError:
                continue
            else:
                gene_ids = gene_ids.union([g_id])
    #

    # Load tasks into task_queue.
    gene_ids_processed = []
    with h5py.File(os.path.join(out_dir,'eventalign.hdf5'),'r') as f:
        for gene_id in gene_ids:
            if resume and (gene_id in gene_ids_done):
                continue
                
            # mapping a gene <-> transcripts
            tx_ids, t2g_mapping = t2g(gene_id,ensembl)
            #
            read_ids = []
            data_dict = dict()
            n_reads = 0
            for tx_id in tx_ids:
                
                if tx_id not in tx_ensembl:
                    continue
                tx_id += '.' + tx_ensembl[tx_id]
        
                if tx_id not in f: # no eventalign for tx_id
                    continue
                    
                # n_reads += len(f[tx_id])                
                for read_id in f[tx_id].keys():
                    if (n_reads < readcount_max) and (read_id not in read_ids):
                        data_dict[read_id] = f[tx_id][read_id]['events'][:]
                        read_ids += [read_id]
                        n_reads += 1
                    elif n_reads >= readcount_max:
                        break
                    
            if len(read_ids) > 0:
                task_queue.put((gene_id,data_dict,t2g_mapping,out_paths)) # Blocked if necessary until a free slot is available. 
                gene_ids_processed += [gene_id]

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

    # Write the ending of the json file.
    # with open(out_paths['json'],'a+') as f:
    #     f.seek(0,2)  # end of file
    #     f.truncate(f.tell()-1) 
    #     f.write('\n}\n}\n')
    ###
    
    with open(out_paths['log'],'a+') as f:
        f.write('Total %d genes.\n' %len(gene_ids_processed))
        f.write(helper.decor_message('successfully finished'))

def preprocess_gene(gene_id,data_dict,t2g_mapping,out_paths,locks):  
    """
    Convert transcriptomic to genomic coordinates for a gene.
    
    Parameters
    ----------
        gene_id: str
            Gene ID.
        data_dict: {read_id:events_array}
            Events for each read.
        t2g_mapping: {(,):()}
            A dict to map transcriptomic coordinates (transcript id and transcriptomic position) to genomic (gene id and genomic position).
        features: [str] # todo
            A list of features to collect from the reads that are aligned to each genomic coordinate in the output.
    Returns
    -------
    dict
        A dict of all specified features collected for each genomic coordinate.
    """
    
    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx'] # columns in the eventalign file per read.

    events = []
    condition_labels = []
    run_labels = []
    read_ids = []
    genomic_coordinates = []
    
    # Concatenate
    if len(data_dict) == 0:
        return

    for read_id,events_per_read in data_dict.items(): 
        # print(read_id)

        # ===== transcript to gene coordinates ===== # TODO: to use gtf.
        
        tx_ids = [tx_id.decode('UTF-8').split('.')[0] for tx_id in events_per_read['transcript_id']] 
        tx_positions = events_per_read['transcriptomic_position']
        
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- np structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = np.array(genomic_coordinate,dtype=np.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        # ===== 
        
        # Based on Ensembl, remove transcript version.
        events_per_read['transcript_id'] = [tx_id.decode('UTF-8').split('.')[0] for tx_id in events_per_read['transcript_id']] 
        events_per_read = np.array(events_per_read,dtype=np.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')]))
        #
        
        events += [events_per_read]
        genomic_coordinates += [genomic_coordinate]
        n_events_per_read = len(events_per_read)
        
    events = np.concatenate(events)
    genomic_coordinates = np.concatenate(genomic_coordinates)
   
    # Sort and split 
    idx_sorted = np.lexsort((events['reference_kmer'],genomic_coordinates['genomic_position'],genomic_coordinates['gene_id']))
    key_tuples, index = np.unique(list(zip(genomic_coordinates['gene_id'][idx_sorted],genomic_coordinates['genomic_position'][idx_sorted],events['reference_kmer'][idx_sorted])),return_index = True,axis=0) #'chr',
    y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
    read_id_arrays = np.split(events['read_id'][idx_sorted], index[1:])
    g_kmer_arrays = np.split(genomic_coordinates['g_kmer'][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
    for key_tuple,y_array,read_id_array,g_kmer_array in zip(key_tuples,y_arrays,read_id_arrays,g_kmer_arrays):
        gene_id,position,kmer = key_tuple
        if (len(set(g_kmer_array)) == 1) and ('XXXXX' in set(g_kmer_array)) or (len(y_array) == 0):
            continue
                        
        data[position] = {kmer: list(np.around(y_array,decimals=2))} #,'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}
        
    # write to file.
    log_str = '%s: Data preparation ... Done.' %(gene_id)

    with locks['json'], open(out_paths['json'],'a') as f:

#         f.write('\n')
#         pos_start = f.tell()
#         f.write('"%s":' %gene_id)
#         json.dump(data, f)
#         pos_end = f.tell()
#         f.write(',')
        pos_start = f.tell()
        f.write('{')
        f.write('"%s":' %gene_id)
        json.dump(data, f)
        f.write('}\n')
        pos_end = f.tell()

    with locks['index'], open(out_paths['index'],'a') as f:
        f.write('%s,%d,%d\n' %(gene_id,pos_start,pos_end))
        
    with locks['readcount'], open(out_paths['readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
        n_reads = len(data_dict)
        f.write('%s,%d\n' %(gene_id,n_reads))
        
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(log_str + '\n')


def parallel_preprocess_tx(out_dir,n_processes,readcount_max,resume):
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['json','index','log','readcount']:
        out_paths[out_filetype] = os.path.join(out_dir,'data.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
                
    # Writing the starting of the files.
    tx_ids_done = []
    if resume and os.path.exists(out_paths['index']):
        df_index = pd.read_csv(out_paths['index'],sep=',')
        tx_ids_done = list(df_index['transcript_id'].unique())
    else:
        open(out_paths['json'],'w').close()
        with open(out_paths['index'],'w') as f:
            f.write('idx,start,end\n') # header
        with open(out_paths['readcount'],'w') as f:
            f.write('idx,n_reads\n') # header
        open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=preprocess_tx,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Load tasks into task_queue.
    tx_ids_processed = []
    with h5py.File(os.path.join(out_dir,'eventalign.hdf5'),'r') as f:
        for tx_id in f.keys():
            if resume and (tx_id in tx_ids_done):
                continue
                
            read_ids = []
            n_reads = 0
            data_dict = dict()
            for read_id in f[tx_id].keys():
                if n_reads < readcount_max:
                    data_dict[read_id] = f[tx_id][read_id]['events'][:]
                    read_ids += [read_id]
                    n_reads += 1
                else:
                    break
            task_queue.put((tx_id,data_dict,out_paths)) # Blocked if necessary until a free slot is available. 
            tx_ids_processed += [tx_id]

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(out_paths['log'],'a+') as f:
        f.write('Total %d transcripts.\n' %len(tx_ids_processed))
        f.write(helper.decor_message('successfully finished'))

def preprocess_tx(tx_id,data_dict,out_paths,locks):  # todo
    """
    Convert transcriptomic to genomic coordinates for a gene.
    
    Parameters
    ----------
        tx_id: str
            Transcript ID.
        data_dict: {read_id:events_array}
            Events for each read.
        features: [str] # todo
            A list of features to collect from the reads that are aligned to each genomic coordinate in the output.
    Returns
    -------
    dict
        A dict of all specified features collected for each genomic coordinate.
    """
    
    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx'] # columns in the eventalign file per read.

    events = []
    condition_labels = []
    run_labels = []
    read_ids = []
    transcriptomic_coordinates = []
    
    # Concatenate
    if len(data_dict) == 0:
        return

    for read_id,events_per_read in data_dict.items(): 
        # print(read_id)
        events += [events_per_read]
        
    events = np.concatenate(events)
   
    # Sort and split 
    idx_sorted = np.lexsort((events['reference_kmer'],events['transcriptomic_position'],events['transcript_id']))
    key_tuples, index = np.unique(list(zip(events['transcript_id'][idx_sorted],events['transcriptomic_position'][idx_sorted],events['reference_kmer'][idx_sorted])),return_index = True,axis=0) #'chr',
    y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
    read_id_arrays = np.split(events['read_id'][idx_sorted], index[1:])
    reference_kmer_arrays = np.split(events['reference_kmer'][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
    for key_tuple,y_array,read_id_array,reference_kmer_array in zip(key_tuples,y_arrays,read_id_arrays,reference_kmer_arrays):
        idx,position,kmer = key_tuple
        position = int(position.decode('UTF-8'))
        kmer = kmer.decode('UTF-8')
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(y_array) == 0):
            continue
                        
        data[position] = {kmer: list(np.around(y_array,decimals=2))} #,'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}
        
    # write to file.
    log_str = '%s: Data preparation ... Done.' %(tx_id)
    with locks['json'], open(out_paths['json'],'a') as f:
        pos_start = f.tell()
        f.write('{')
        f.write('"%s":' %tx_id)
        json.dump(data, f)
        f.write('}\n')
        pos_end = f.tell()
        
    with locks['index'], open(out_paths['index'],'a') as f:
        f.write('%s,%d,%d\n' %(tx_id,pos_start,pos_end))
        
    with locks['readcount'], open(out_paths['readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
        n_reads = len(data_dict)
        f.write('%s,%d\n' %(tx_id,n_reads))
        
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(log_str + '\n')
        

def main():
    args = get_args()
    #
    n_processes = args.n_processes        
    eventalign_filepath = args.eventalign
    summary_filepath = args.summary
    out_dir = args.out_dir
    ensembl_version = args.ensembl
    ensembl_species = args.species
    readcount_max = args.readcount_max
    resume = args.resume
    genome = args.genome


    misc.makedirs(out_dir) #todo: check every level.
    
    # (1) For each read, combine multiple events aligned to the same positions, the results from nanopolish eventalign, into a single event per position.
    eventalign_log_filepath = os.path.join(out_dir,'eventalign.log')
    if not helper.is_successful(eventalign_log_filepath):
        parallel_combine(eventalign_filepath,summary_filepath,out_dir,n_processes,resume)
    
    # (2) Create a .json file, where the info of all reads are stored per position, for modelling.
    if genome:
        ensembl = EnsemblRelease(ensembl_version,ensembl_species) # Default: human reference genome GRCh38 release 91 used in the ont mapping.    
        parallel_preprocess_gene(ensembl,out_dir,n_processes,readcount_max,resume)

    else:
        parallel_preprocess_tx(out_dir,n_processes,readcount_max,resume)


if __name__ == '__main__':
    """
    Usage:
        --eventalign EVENTALIGN --summary SUMMARY --out_dir OUT_DIR \
                      [--ensembl ENSEMBL] [--genome] \
                      [--n_processes N_PROCESSES] \
                      [--readcount_max READCOUNT_MAX] [--resume] \
    """
    main()


