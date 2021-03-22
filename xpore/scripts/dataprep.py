import argparse
import numpy as np
import pandas as pd
import os,re
import multiprocessing 
import h5py
import csv
import ujson
import pickle
from operator import itemgetter
from collections import defaultdict
from itertools import groupby
from io import StringIO

from . import helper
from .constants import M6A_KMERS, NUM_NEIGHBORING_FEATURES
from ..utils import misc

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required arguments
    required.add_argument('--eventalign', dest='eventalign', help='eventalign filepath, the output from nanopolish.',required=True)
    required.add_argument('--summary', dest='summary', help='eventalign summary filepath, the output from nanopolish.',required=True)
    required.add_argument('--out_dir', dest='out_dir', help='output directory.',required=True)
    optional.add_argument('--gtf_path_or_url', dest='gtf_path_or_url', help='gtf file path or url.',type=str)
    optional.add_argument('--transcript_fasta_paths_or_urls', dest='transcript_fasta_paths_or_urls', help='transcript fasta paths or urls.',type=str)
    optional.add_argument('--program', dest='program', help='program name (xpore and/or m6anet).',type=str, default='xpore,m6anet')

    # Optional
    optional.add_argument('--skip_eventalign_indexing', dest='skip_eventalign_indexing', help='skip indexing the eventalign nanopolish output.',default=False,action='store_true')

    # parser.add_argument('--features', dest='features', help='Signal features to extract.',type=list,default=['norm_mean'])
    optional.add_argument('--genome', dest='genome', help='to run on Genomic coordinates. Without this argument, the program will run on transcriptomic coordinates',default=False,action='store_true') 
    optional.add_argument('--n_processes', dest='n_processes', help='number of processes to run.',type=int, default=1)
    optional.add_argument('--chunk_size', dest='chunk_size', help='number of lines from nanopolish eventalign.txt for processing.',type=int, default=1000000)
    optional.add_argument('--readcount_min', dest='readcount_min', help='minimum read counts per gene.',type=int, default=1)
    optional.add_argument('--readcount_max', dest='readcount_max', help='maximum read counts per gene.',type=int, default=1000)
    optional.add_argument('--resume', dest='resume', help='with this argument, the program will resume from the previous run.',default=False,action='store_true') #todo
    ##for m6anet
    optional.add_argument('--n_neighbors', dest='n_neighbors', help='number of neighboring features to extract.',type=int, default=NUM_NEIGHBORING_FEATURES)
    parser._action_groups.append(optional)
    return parser.parse_args()

####### m6anet processing specifics
def partition_into_continuous_positions(arr, window_size=1):
    arr = arr[np.argsort(arr["transcriptomic_position"])]
    float_features = ['dwell_time', 'norm_std', 'norm_mean']
    float_dtypes = [('norm_mean', '<f8'), ('norm_std', '<f8'), ('dwell_time', '<f8')]
    
    float_arr = arr[float_features].astype(float_dtypes).view('<f8').reshape(-1, 3)
    kmer_arr = arr["reference_kmer"].reshape(-1, 1)
    tx_pos_arr = arr["transcriptomic_position"]
    tx_id_arr = arr["transcript_id"]

    partitions = [list(map(itemgetter(0), g)) for k, g in groupby(enumerate(tx_pos_arr), 
                                                                  lambda x: x[0] - x[1])]
    return [(float_arr[partition],
             kmer_arr[partition], tx_id_arr[partition], tx_pos_arr[partition]) 
            for partition in partitions if len(partition) > 2 * window_size + 1]

def filter_by_kmer(partition, kmers, window_size):
    feature_arr, kmer_arr, tx_id_arr, tx_pos_arr = partition
    kmers_5 = kmer_arr[:, (2 * window_size + 1) // 2]
    mask = np.isin(kmers_5, kmers)
    filtered_feature_arr = feature_arr[mask, :]
    filtered_kmer_arr = kmer_arr[mask, :]
    filtered_tx_pos_arr = tx_pos_arr[mask]
    filtered_tx_id_arr = tx_id_arr[mask]

    if len(filtered_kmer_arr) == 0:
        return []
    else:
        return filtered_feature_arr, filtered_kmer_arr, filtered_tx_id_arr, filtered_tx_pos_arr

def filter_partitions(partitions, window_size, kmers):
    windowed_partition = [create_features(partition, window_size) for partition in partitions]
    filtered_by_kmers = [filter_by_kmer(partition, kmers, window_size) for partition in windowed_partition]
    final_partitions = [x for x in filtered_by_kmers if len(x) > 0]
    return final_partitions

def roll(to_roll, window_size=1):
    nex = np.concatenate([np.roll(to_roll, i, axis=0) for i in range(-1, - window_size - 1, -1)],
                          axis=1)
    prev = np.concatenate([np.roll(to_roll, i, axis=0) for i in range(window_size, 0, -1)], axis=1)
    return np.concatenate((prev, to_roll, nex), axis=1)[window_size: -window_size, :]

def create_features(partition, window_size=1):
    float_arr, kmer_arr, tx_id_arr, tx_pos_arr = partition
    return roll(float_arr, window_size), roll(kmer_arr, window_size), \
        tx_id_arr[window_size: -window_size], tx_pos_arr[window_size: -window_size]

def filter_events(events, window_size, kmers):
    events = partition_into_continuous_positions(events)
    events = filter_partitions(events, window_size, kmers)
    return events

def combine_sequence(kmers):
    kmer = kmers[0]
    for _kmer in kmers[1:]:
        kmer += _kmer[-1]
    return kmer
########################################################################

def index(eventalign_result,pos_start,out_paths,locks):
    eventalign_result = eventalign_result.set_index(['contig','read_index'])
    pos_end=pos_start
    with locks['index'], open(out_paths['index'],'a') as f_index:
        for index in list(dict.fromkeys(eventalign_result.index)):
            transcript_id,read_index = index
            pos_end += eventalign_result.loc[index]['line_length'].sum()
            
            try: # sometimes read_index is nan
                f_index.write('%s,%d,%d,%d\n' %(transcript_id,read_index,pos_start,pos_end))
            except:
                pass
            pos_start = pos_end

def parallel_index(eventalign_filepath,summary_filepath,chunk_size,out_dir,n_processes,resume):
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['index']:
        out_paths[out_filetype] = os.path.join(out_dir,'eventalign.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
        
        
#     read_names_done = []
#     if resume and os.path.exists(out_paths['log']):
#         read_names_done = [line.rstrip('\n') for line in open(out_paths['log'],'r')]
#     else:
        # Create empty files.
    with open(out_paths['index'],'w') as f:
        f.write('transcript_id,read_index,pos_start,pos_end\n') # header


    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=index,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()
        
    ## Load tasks into task_queue. A task is eventalign information of one read.
    eventalign_file = open(eventalign_filepath,'r')
    pos_start = len(eventalign_file.readline()) #remove header
    chunk_split = None
    index_features = ['contig','read_index','line_length']
    for chunk in pd.read_csv(eventalign_filepath, chunksize=chunk_size,sep='\t'):
        chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
        chunk_concat = pd.concat([chunk_split,chunk_complete])
        chunk_concat_size = len(chunk_concat.index)
        ## read the file at where it left off because the file is opened once ##
        lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
        chunk_concat['line_length'] = np.array(lines)
        task_queue.put((chunk_concat[index_features],pos_start,out_paths))
        pos_start += sum(lines)
        chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']]
    ## the loop above leaves off w/o adding the last read_index to eventalign.index
    chunk_split_size = len(chunk_split.index)
    lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
    chunk_split['line_length'] = np.array(lines)
    task_queue.put((chunk_split[index_features],pos_start,out_paths))

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

def readFasta(transcript_fasta_paths_or_urls):
    fasta=open(transcript_fasta_paths_or_urls,"r")
    entries=""
    for ln in fasta:
        entries+=ln
    entries=entries.split(">")
    dict={}
    for entry in entries:
        entry=entry.split("\n")
        id=entry[0].split(".")[0]
        seq="".join(entry[1:])
        dict[id]=seq
    with open(transcript_fasta_paths_or_urls+'.pickle', 'wb') as fasta_pickle:
        pickle.dump(dict, fasta_pickle)
    return dict

def readGTF(gtf_path_or_url):
    gtf=open(gtf_path_or_url,"r")
    for i in range(5):
        gtf.readline()
    dict={}
    for ln in gtf:
        ln=ln.split("\t")
        if ln[2] == "transcript" or ln[2] == "exon":
            chr,type,start,end=ln[0],ln[2],int(ln[3]),int(ln[4])
            tx_id=ln[-1].split('; transcript_id "')[1].split('";')[0]
            g_id=ln[-1].split('gene_id "')[1].split('";')[0]
            if tx_id not in dict:
                dict[tx_id]={'chr':chr,'g_id':g_id,'strand':ln[6]}
                if type not in dict[tx_id]:
                    if type == "transcript":
                        dict[tx_id][type]=(start,end)
            else:
                if type == "exon":
                    if type not in dict[tx_id]:
                        dict[tx_id][type]=[(start,end)]
                    else:
                        dict[tx_id][type].append((start,end))
    #convert genomic positions to tx positions
    for id in dict:
        tx_pos,tx_start=[],0
        for pair in dict[id]["exon"]:
            tx_end=pair[1]-pair[0]+tx_start
            tx_pos.append((tx_start,tx_end))
            tx_start=tx_end+1
        dict[id]['tx_exon']=tx_pos
    with open(gtf_path_or_url+'.pickle', 'wb') as gtf_pickle:
        pickle.dump(dict, gtf_pickle)
    return dict

def t2g(gene_id,fasta_dict,gtf_dict,g2t_mapping,df_eventalign_index,readcount_min):
    tx_ids = []
    t2g_dict = {}
    transcripts = [tx for tx in gtf_dict if tx in g2t_mapping[gene_id]]
    n_reads = sum([len(df_eventalign_index.loc[tx]) for tx in transcripts])
    if n_reads >= readcount_min:
        for tx in transcripts:
            tx_seq = fasta_dict[tx]
            tx_contig = gtf_dict[tx]['chr']
            if tx_seq is None:
                continue
            for exon_num in range(len(gtf_dict[tx]['exon'])):
                g_interval=gtf_dict[tx]['exon'][exon_num]
                tx_interval=gtf_dict[tx]['tx_exon'][exon_num]
                for g_pos in range(g_interval[0],g_interval[1]+1): # Exclude the rims of exons.
                    dis_from_start = g_pos - g_interval[0]
                    if gtf_dict[tx]['strand'] == "+":
                        tx_pos = tx_interval[0] + dis_from_start
                    elif gtf_dict[tx]['strand'] == "-":
                        tx_pos = tx_interval[1] - dis_from_start
                    if (g_interval[0] <= g_pos < g_interval[0]+2) or (g_interval[1]-2 < g_pos <= g_interval[1]): # Todo: To improve the mapping
                        kmer = 'XXXXX'
                    else:
                        kmer = tx_seq[tx_pos-2:tx_pos+3]
                    t2g_dict[(tx,tx_pos)] = (tx_contig,gene_id,g_pos,kmer) # tx.contig is chromosome.
            tx_ids += [tx]

                
    return n_reads, tx_ids, t2g_dict

def combine_xpore(events_str):
    f_string = StringIO(events_str)
    eventalign_result = pd.read_csv(f_string,delimiter='\t',names=['contig','position','reference_kmer','read_index',
                         'strand','event_index','event_level_mean','event_stdv','event_length','model_kmer',
                         'model_mean', 'model_stdv', 'standardized_level', 'start_idx', 'end_idx'])
    f_string.close()
    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    if cond_successfully_eventaligned.sum() != 0:

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
        eventalign_result['norm_mean'] = (sum_norm_mean/total_length).round(1)

        eventalign_result.reset_index(inplace=True)


        eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]    #### CHANGE MADE ####
        #eventalign_result['transcript_id'] = eventalign_result['contig']

        eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
        # eventalign_result = misc.str_encode(eventalign_result)
#         eventalign_result['read_id'] = [read_name]*len(eventalign_result)

        # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx']
        # features_dtype = np.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'), ('end_idx', '<i8')])
        
#         features = ['transcript_id','transcriptomic_position','reference_kmer','norm_mean']

#         df_events = eventalign_result[['read_index']+features]
#         # print(df_events.head())

        features = ['transcript_id','transcriptomic_position','reference_kmer','norm_mean']
#        np_events = eventalign_result[features].reset_index().values.ravel().view(dtype=[('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')])
        df_events = eventalign_result[features]
        np_events = np.rec.fromrecords(df_events, names=[*df_events])
        return np_events

def combine_m6anet(events_str):
    f_string = StringIO(events_str)
    eventalign_result = pd.read_csv(f_string,delimiter='\t',names=['contig','position','reference_kmer','read_index','strand','event_index','event_level_mean','event_stdv','event_length','model_kmer','model_mean','model_stdv','standardized_level','start_idx','end_idx'])
    f_string.close()
    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    if cond_successfully_eventaligned.sum() != 0:

        eventalign_result = eventalign_result[cond_successfully_eventaligned]

        keys = ['read_index','contig','position','reference_kmer'] # for groupby
        eventalign_result['length'] = pd.to_numeric(eventalign_result['end_idx'])-pd.to_numeric(eventalign_result['start_idx'])
        eventalign_result['sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']
        eventalign_result['sum_norm_std'] = pd.to_numeric(eventalign_result['event_stdv']) * eventalign_result['length']
        eventalign_result['sum_dwell_time'] = pd.to_numeric(eventalign_result['event_length']) * eventalign_result['length']
            
        eventalign_result = eventalign_result.groupby(keys)  
        sum_norm_mean = eventalign_result['sum_norm_mean'].sum() 
        sum_norm_std = eventalign_result["sum_norm_std"].sum()
        sum_dwell_time = eventalign_result["sum_dwell_time"].sum()

        start_idx = eventalign_result['start_idx'].min()
        end_idx = eventalign_result['end_idx'].max()
        total_length = eventalign_result['length'].sum()

        eventalign_result = pd.concat([start_idx,end_idx],axis=1)
        eventalign_result['norm_mean'] = (sum_norm_mean/total_length).round(1)
        eventalign_result["norm_std"] = sum_norm_std / total_length
        eventalign_result["dwell_time"] = sum_dwell_time / total_length
        eventalign_result.reset_index(inplace=True)


        eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]    #### CHANGE MADE ####
        #eventalign_result['transcript_id'] = eventalign_result['contig']

        eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
        # eventalign_result = misc.str_encode(eventalign_result)
#         eventalign_result['read_id'] = [read_name]*len(eventalign_result)

        # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx']
        # features_dtype = np.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'), ('end_idx', '<i8')])
        
#         features = ['transcript_id','transcriptomic_position','reference_kmer','norm_mean']

#         df_events = eventalign_result[['read_index']+features]
#         # print(df_events.head())

        features = ['transcript_id','read_index','transcriptomic_position','reference_kmer','norm_mean','norm_std','dwell_time']
#        np_events = eventalign_result[features].reset_index().values.ravel().view(dtype=[('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')])
        df_events = eventalign_result[features]
        np_events = np.rec.fromrecords(df_events, names=[*df_events])
        return np_events
    else:
        return np.array([])

def writeOutputs_xpore(id,asserted,data,data_dict,locks,out_paths):
    log_str = '%s: %s' %(id,asserted)
    with locks['xpore.json'], open(out_paths['xpore.json'],'a') as f:

        pos_start = f.tell()
        f.write('{')
        f.write('"%s":' %id)
        ujson.dump(data, f)
        f.write('}\n')
        pos_end = f.tell()

    with locks['xpore.index'], open(out_paths['xpore.index'],'a') as f:
        f.write('%s,%d,%d\n' %(id,pos_start,pos_end))
        
    with locks['xpore.readcount'], open(out_paths['xpore.readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
        n_reads = len(data_dict)
        f.write('%s,%d\n' %(id,n_reads))
        
    with locks['xpore.log'], open(out_paths['xpore.log'],'a') as f:
        f.write(log_str + '\n')

def writeOutputs_m6anet(id,data,locks,out_paths):
    log_str = '%s: Data preparation ... Done.' %(id)
    with locks['m6anet.json'], open(out_paths['m6anet.json'],'a') as f, \
            locks['m6anet.index'], open(out_paths['m6anet.index'],'a') as g, \
            locks['m6anet.readcount'], open(out_paths['m6anet.readcount'],'a') as h:
        
        for pos, dat in data.items():
            pos_start = f.tell()
            f.write('{')
            f.write('"%s":{"%d":' %(id,pos))
            ujson.dump(dat, f)
            f.write('}}\n')
            pos_end = f.tell()
        
            # with locks['m6anet.index'], open(out_paths['m6anet.index'],'a') as f:
            g.write('%s,%d,%d,%d\n' %(id,pos,pos_start,pos_end))
        
            # with locks['m6anet.readcount'], open(out_paths['m6anet.readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
            n_reads = 0
            for kmer, features in dat.items():
                n_reads += len(features)
            h.write('%s,%d,%d\n' %(id,pos,n_reads))
        
    with locks['m6anet.log'], open(out_paths['m6anet.log'],'a') as f:
        f.write(log_str + '\n')

def parallel_preprocess_gene(program,eventalign_filepath,fasta_dict,gtf_dict,fun_dict,out_dir,n_processes,readcount_min,readcount_max,n_neighbors,resume):
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    gene_ids_done = []
    for p in program:
        for base_out_filetype in ['json','index','log','readcount']:
            out_filetype= p+"."+base_out_filetype ##out_filetype specific to either xpore or m6anet
            out_paths[out_filetype] = os.path.join(out_dir,'data.%s' %out_filetype)
            locks[out_filetype] = multiprocessing.Lock()
        json_tag,index_tag,readcount_tag,log_tag=p+'.json',p+'.index',p+'.readcount',p+'.log'
        if resume and os.path.exists(out_paths[index_tag]):
            df_index = pd.read_csv(out_paths[index_tag],sep=',')
            gene_ids_done = list(df_index['transcript_id'].unique())
        else:
            open(out_paths[json_tag],'w').close()
            with open(out_paths[index_tag],'w') as f:
                f.write('idx,start,end\n') # header
            with open(out_paths[readcount_tag],'w') as f:
                f.write('idx,n_reads\n') # header
            open(out_paths[log_tag],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=fun_dict['preprocess_gene'],locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Get all gene ids and create a dict of eventalign.combine index.
#     gene_ids = set()

    df_eventalign_index = pd.read_csv(os.path.join(out_dir,'eventalign.index'))
    df_eventalign_index['transcript_id'] = [tx_id.split('.')[0] for tx_id in  df_eventalign_index['transcript_id']]
    df_eventalign_index.set_index('transcript_id',inplace=True)
    g2t_mapping = defaultdict(list)

    for tx_id in set(df_eventalign_index.index):
        try:
##           g_id = ensembl.transcript_by_id(tx_id).gene_id 
            g_id = gtf_dict[tx_id]['g_id'] 
        except ValueError:
            continue
        else:
#             gene_ids = gene_ids.union([g_id])
            g2t_mapping[g_id] += [tx_id]

#     f = open(os.path.join(out_dir,'eventalign.index'))
#     for ln in f:
#         tx_id,read_index,pos_start,pos_end = ln.split(',')
#         tx_id,tx_version = tx_id.split('.') # Based on Ensembl
#         eventalign_index[tx_id] += [(int(read_index),int(pos_start),int(pos_end))]
#         tx_ensembl[tx_id] = tx_version
#         try:
#             g_id = ensembl.transcript_by_id(tx_id).gene_id 
#         except ValueError:
#             continue
#         else:
#             gene_ids = gene_ids.union([g_id])
            

    # Load tasks into task_queue.    
    gene_ids_processed = []

    with open(eventalign_filepath,'r') as eventalign_result:

        for gene_id in g2t_mapping:
                        
            if resume and (gene_id in gene_ids_done):
                continue
            # mapping a gene <-> transcripts

            n_reads, tx_ids, t2g_mapping = t2g(gene_id,fasta_dict,gtf_dict,g2t_mapping,df_eventalign_index,readcount_min)
            #
            if n_reads >= readcount_min: 
                data_dict = dict()
                readcount = 0
                for tx_id in tx_ids:
                    for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                        read_index,pos_start,pos_end = row['read_index'],row['pos_start'],row['pos_end']
                        eventalign_result.seek(pos_start,0)
                        events_str = eventalign_result.read(pos_end-pos_start)
                        data = fun_dict['combine'](events_str)
                        #data = np.genfromtxt(f_string,delimiter=',',dtype=np.dtype([('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')]))
                        if (data is not None) and (data.size > 1):
                            data_dict[read_index] = data
                        readcount += 1 
                        if readcount > readcount_max:
                            break

                    if readcount > readcount_max:
                        break
                if len(data_dict)>=readcount_min:
#                     print(gene_id,len(data_dict)) #len(data_dict) is the number of reads to be processed.
                    task_queue.put((gene_id,data_dict,t2g_mapping,n_neighbors,out_paths)) # Blocked if necessary until a free slot is available. 
                    gene_ids_processed += [gene_id]


    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    for p in program:
        log_tag = p+'.log'
        with open(out_paths[log_tag],'a+') as f:
            f.write('Total %d genes.\n' %len(gene_ids_processed))
            f.write(helper.decor_message('successfully finished'))

def preprocess_gene_xpore(gene_id,data_dict,t2g_mapping,n_neighbors,out_paths,locks):  
    """
    Convert transcriptomic to genomic coordinates for a gene.
    
    Parameters
    ----------
        gene_id: str
            Gene ID.
        data_dict: {tx_id:events_array}
            Events for each read.
        t2g_mapping: {(,):()}
            A dict to map transcriptomic coordinates (transcript id and transcriptomic position) to genomic (gene id and genomic position).
        db_type: 
            Type of gene-tx mapping either EnsemblRelease or (customised) Genome 
        features: [str] # todo
            A list of features to collect from the reads that are aligned to each genomic coordinate in the output.
    Returns
    -------
    dict
        A dict of all specified features collected for each genomic coordinate.
    """
    
    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx'] # columns in the eventalign file per read.

    events = []
#    condition_labels = []
#    run_labels = []
#    read_ids = []
    genomic_coordinates = []
    
    # Concatenate
#     if len(data_dict) == 0:
#         return


    for read_index,events_per_read in data_dict.items():
#         if len(events_per_read) > 0:
        # ===== transcript to gene coordinates ===== # TODO: to use gtf.
#        tx_ids = [tx_id.decode('UTF-8').split('.')[0] for tx_id in events_per_read['transcript_id']]
        tx_ids = [tx_id for tx_id in events_per_read['transcript_id']] 
        tx_positions = events_per_read['transcriptomic_position']
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- np structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = np.array(genomic_coordinate,dtype=np.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        # ===== 

        # Based on Ensembl, remove transcript version.

        events_per_read['transcript_id'] = tx_ids
        events_per_read = np.array(events_per_read,dtype=np.dtype([('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')]))

        #

        events += [events_per_read]
        genomic_coordinates += [genomic_coordinate]
#        n_events_per_read = len(events_per_read)
#         else:
#             print(read_index,len(events_per_read))

    events = np.concatenate(events)
    genomic_coordinates = np.concatenate(genomic_coordinates)
   
    # Sort and split # 
#     idx_sorted = np.lexsort((events['reference_kmer'],genomic_coordinates['genomic_position'],genomic_coordinates['gene_id']))
#     key_tuples, index = np.unique(list(zip(genomic_coordinates['gene_id'][idx_sorted],genomic_coordinates['genomic_position'][idx_sorted],events['reference_kmer'][idx_sorted])),return_index = True,axis=0) #'chr',
#     y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
# #     read_id_arrays = np.split(events['read_id'][idx_sorted], index[1:])
#     g_kmer_arrays = np.split(genomic_coordinates['g_kmer'][idx_sorted], index[1:])

    idx_sorted = np.argsort(genomic_coordinates['genomic_position'])
    unique_positions, index = np.unique(genomic_coordinates['genomic_position'][idx_sorted],return_index = True)
    y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
    #     read_id_arrays = np.split(events['read_id'][idx_sorted], index[1:])
    g_kmer_arrays = np.split(genomic_coordinates['g_kmer'][idx_sorted], index[1:])
    g_positions_arrays = np.split(genomic_coordinates['genomic_position'][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
#     data = dict(zip(key_tuples, y_arrays))

    asserted = True
#     for key_tuple,y_array,g_kmer_array in zip(key_tuples,y_arrays,g_kmer_arrays):
    for position,y_array,g_kmer_array,g_positions_array in zip(unique_positions,y_arrays,g_kmer_arrays,g_positions_arrays):
#         gene_id,position,kmer = key_tuple            
        if (len(set(g_kmer_array)) == 1) and ('XXXXX' in set(g_kmer_array)) or (len(y_array) == 0):
            continue
            
        if 'XXXXX' in set(g_kmer_array):
            y_array = y_array[g_kmer_array != 'XXXXX']  
            assert len(y_array) == len(g_kmer_array) - (g_kmer_array=='XXXXX').sum()
            g_kmer_array = g_kmer_array[g_kmer_array != 'XXXXX']  
            
        try:
            assert len(set(g_kmer_array)) == 1
            assert {position} == set(g_positions_array)
        except:
            asserted = False
            break
        kmer = set(g_kmer_array).pop()

        data[position] = {kmer: list(y_array)} #,'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}
        
    # write to file.
    writeOutputs_xpore(gene_id,asserted,data,data_dict,locks,out_paths)

def preprocess_gene_m6anet(gene_id,data_dict,t2g_mapping,n_neighbors,out_paths,locks):  # todo
    """
    Convert transcriptomic to genomic coordinates for a gene.
    
    Parameters
    ----------
        gene_id: str
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

    # Concatenate
    if len(data_dict) == 0:
        return
    
    features_arrays = []
    reference_kmer_arrays = []
    g_positions_arrays = []

    for read_id,events_per_read in data_dict.items(): 
        # print(read_id)
        tx_ids = [tx_id for tx_id in events_per_read['transcript_id']]
        tx_positions = events_per_read['transcriptomic_position']
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- np structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = np.array(genomic_coordinate,dtype=np.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        gene_ids = [gene_id for gene_id in genomic_coordinate['gene_id']]
        genomic_positions = [genomic_position for genomic_position in genomic_coordinate['genomic_position']]
        ## note the following changes transcript_id -> gene_id and transcriptomic_position -> genomic_position ##
        events_per_read['transcript_id'] = gene_ids  ##convert the transcript_id field as gene_id
        events_per_read['transcriptomic_position'] = genomic_positions ##convert the transcriptomic_position as genomic_position
        events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS)
        for event_per_read in events_per_read:
            features_arrays.append(event_per_read[0])
            reference_kmer_arrays.append([combine_sequence(kmer) for kmer in event_per_read[1]])
            g_positions_arrays.append(event_per_read[3])

    if len(features_arrays) == 0:
        return
    else:
        features_arrays = np.concatenate(features_arrays)
        reference_kmer_arrays = np.concatenate(reference_kmer_arrays)
        g_positions_arrays = np.concatenate(g_positions_arrays)
        assert(len(features_arrays) == len(reference_kmer_arrays) == len(g_positions_arrays))
    # Sort and split

    idx_sorted = np.argsort(g_positions_arrays)
    positions, index = np.unique(g_positions_arrays[idx_sorted], return_index = True,axis=0) #'chr',
    features_arrays = np.split(features_arrays[idx_sorted], index[1:])
    reference_kmer_arrays = np.split(reference_kmer_arrays[idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)


    # for each position, make it ready for json dump
    for position, features_array, reference_kmer_array in zip(positions, features_arrays, reference_kmer_arrays):
        kmer = set(reference_kmer_array)
     #   assert(len(kmer) == 1) ##AssertionError rose
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(features_array) == 0):
            continue

        data[int(position)] = {kmer.pop(): features_array.tolist()}

    # write to file.
    writeOutputs_m6anet(gene_id,data,locks,out_paths)

def preprocess_gene_xpore_m6anet(gene_id,data_dict,t2g_mapping,n_neighbors,out_paths,locks):  # todo
    """
    Convert transcriptomic to genomic coordinates for a gene.
    
    Parameters
    ----------
        gene_id: str
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

    # Concatenate
    if len(data_dict) == 0:
        return
    ##for xpore
    xpore_events = []
    xpore_genomic_coordinates = []
    ##for m6anet
    m6anet_features_arrays = []
    m6anet_reference_kmer_arrays = [] 
    m6anet_g_positions_arrays = []

    for read_id,events_per_read in data_dict.items(): 
        tx_ids = [tx_id for tx_id in events_per_read['transcript_id']]
        tx_positions = events_per_read['transcriptomic_position']
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- np structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = np.array(genomic_coordinate,dtype=np.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        ##for xpore
        events_per_read['transcript_id'] = tx_ids
#        xpore_events_per_read = np.array(events_per_read,dtype=np.dtype([('transcript_id', 'S15'), ('read_id', '<i8'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('stdv', '<f8'), ('dwelltime', '<f8')]))
        xpore_events += [events_per_read]
        xpore_genomic_coordinates += [genomic_coordinate]
        ##for m6anet
        gene_ids = [gene_id for gene_id in genomic_coordinate['gene_id']]
        genomic_positions = [genomic_position for genomic_position in genomic_coordinate['genomic_position']]
        ## note the following changes transcript_id -> gene_id and transcriptomic_position -> genomic_position ##
        events_per_read['transcript_id'] = gene_ids  ##convert the transcript_id field as gene_id
        events_per_read['transcriptomic_position'] = genomic_positions ##convert the transcriptomic_position as genomic_position
        events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS)
        for event_per_read in events_per_read:
            m6anet_features_arrays.append(event_per_read[0])
            m6anet_reference_kmer_arrays.append([combine_sequence(kmer) for kmer in event_per_read[1]])
            m6anet_g_positions_arrays.append(event_per_read[3])

    ##for xpore
    xpore_events = np.concatenate(xpore_events)
    xpore_genomic_coordinates = np.concatenate(xpore_genomic_coordinates)
    ##for m6anet
    if len(m6anet_features_arrays) == 0:
        return
    else:
        m6anet_features_arrays = np.concatenate(m6anet_features_arrays)
        m6anet_reference_kmer_arrays = np.concatenate(m6anet_reference_kmer_arrays)
        m6anet_g_positions_arrays = np.concatenate(m6anet_g_positions_arrays)
        assert(len(m6anet_features_arrays) == len(m6anet_reference_kmer_arrays) == len(m6anet_g_positions_arrays))

    ##for xpore # Sort and split 
    xpore_idx_sorted = np.argsort(xpore_genomic_coordinates['genomic_position'])
    xpore_unique_positions, xpore_index = np.unique(xpore_genomic_coordinates['genomic_position'][xpore_idx_sorted],return_index = True)
    xpore_y_arrays = np.split(xpore_events['norm_mean'][xpore_idx_sorted], xpore_index[1:])
    xpore_g_kmer_arrays = np.split(xpore_genomic_coordinates['g_kmer'][xpore_idx_sorted], xpore_index[1:])
    xpore_g_positions_arrays = np.split(xpore_genomic_coordinates['genomic_position'][xpore_idx_sorted], xpore_index[1:])

    ##for m6anet # Sort and split
    m6anet_idx_sorted = np.argsort(m6anet_g_positions_arrays)
    m6anet_positions, m6anet_index = np.unique(m6anet_g_positions_arrays[m6anet_idx_sorted], return_index = True,axis=0) #'chr',
    m6anet_features_arrays = np.split(m6anet_features_arrays[m6anet_idx_sorted], m6anet_index[1:])
    m6anet_reference_kmer_arrays = np.split(m6anet_reference_kmer_arrays[m6anet_idx_sorted], m6anet_index[1:])

    ##for xpore # Prepare
    xpore_data = defaultdict(dict)
    # for each position, make it ready for json dump
    asserted = True
    for position,y_array,g_kmer_array,g_positions_array in zip(xpore_unique_positions,xpore_y_arrays,xpore_g_kmer_arrays,xpore_g_positions_arrays):
        if (len(set(g_kmer_array)) == 1) and ('XXXXX' in set(g_kmer_array)) or (len(y_array) == 0):
            continue
            
        if 'XXXXX' in set(g_kmer_array):
            y_array = y_array[g_kmer_array != 'XXXXX']  
            assert len(y_array) == len(g_kmer_array) - (g_kmer_array=='XXXXX').sum()
            g_kmer_array = g_kmer_array[g_kmer_array != 'XXXXX']  
            
        try:
            assert len(set(g_kmer_array)) == 1
            assert {position} == set(g_positions_array)
        except:
            asserted = False
            break
        kmer = set(g_kmer_array).pop()

        xpore_data[position] = {kmer: list(y_array)} #,'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}

    ##for m6anet # Prepare
    m6anet_data = defaultdict(dict)
    # for each position, make it ready for json dump
    for position, features_array, reference_kmer_array in zip(m6anet_positions, m6anet_features_arrays, m6anet_reference_kmer_arrays):
        kmer = set(reference_kmer_array)
     #   assert(len(kmer) == 1) ##AssertionError rose
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(features_array) == 0):
            continue

        m6anet_data[int(position)] = {kmer.pop(): features_array.tolist()}

    ##for xpore # write to file.
    writeOutputs_xpore(gene_id,asserted,xpore_data,data_dict,locks,out_paths)

    ##for m6anet # write to file.
    writeOutputs_m6anet(gene_id,m6anet_data,locks,out_paths)

def parallel_preprocess_tx(program,eventalign_filepath,fun_dict,out_dir,n_processes,readcount_min,readcount_max,n_neighbors,resume):
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    tx_ids_done = []
    for p in program:
        for base_out_filetype in ['json','index','log','readcount']:
            out_filetype= p+"."+base_out_filetype ##out_filetype specific to either xpore or m6anet
            out_paths[out_filetype] = os.path.join(out_dir,'data.%s' %out_filetype)
            locks[out_filetype] = multiprocessing.Lock()
        json_tag,index_tag,readcount_tag,log_tag=p+'.json',p+'.index',p+'.readcount',p+'.log'
        if resume and os.path.exists(out_paths[index_tag]):
            df_index = pd.read_csv(out_paths[index_tag],sep=',')
            tx_ids_done = list(df_index['transcript_id'].unique())
        else:
            open(out_paths[json_tag],'w').close()
            with open(out_paths[index_tag],'w') as f:
                f.write('idx,start,end\n') # header
            with open(out_paths[readcount_tag],'w') as f:
                f.write('idx,n_reads\n') # header
            open(out_paths[log_tag],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=fun_dict['preprocess_tx'],locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Load tasks into task_queue.
    tx_ids_processed = []
    df_eventalign_index = pd.read_csv(os.path.join(out_dir,'eventalign.index'))
    df_eventalign_index['transcript_id'] = [tx_id.split('.')[0] for tx_id in  df_eventalign_index['transcript_id']]
    tx_ids = df_eventalign_index['transcript_id'].values.tolist()
    tx_ids = list(dict.fromkeys(tx_ids))
    df_eventalign_index.set_index('transcript_id',inplace=True)
    with open(eventalign_filepath,'r') as eventalign_result:
        for tx_id in tx_ids:
            data_dict = dict()
            readcount = 0
            for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                read_index,pos_start,pos_end = row['read_index'],row['pos_start'],row['pos_end']
                eventalign_result.seek(pos_start,0)
                events_str = eventalign_result.read(pos_end-pos_start)
                data = fun_dict['combine'](events_str)
                if data.size > 1:
                    data_dict[read_index] = data
                readcount += 1 
                if readcount > readcount_max:
                    break
            if readcount>=readcount_min:
                task_queue.put((tx_id,data_dict,n_neighbors,out_paths))
                tx_ids_processed += [tx_id]

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    for p in program:
        log_tag = p+'.log'
        with open(out_paths[log_tag],'a+') as f:
            f.write('Total %d transcripts.\n' %len(tx_ids_processed))
            f.write(helper.decor_message('successfully finished'))

def preprocess_tx_xpore(tx_id,data_dict,n_neighbors,out_paths,locks): 
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
#    condition_labels = []
#    run_labels = []
#    read_ids = []
#    transcriptomic_coordinates = []
    
    # Concatenate
    if len(data_dict) == 0:
        return

    for read_id,events_per_read in data_dict.items(): 
        # print(read_id)
        events += [events_per_read]
        
    events = np.concatenate(events)
   
    # Sort and split 
    idx_sorted = np.argsort(events['transcriptomic_position'])
    unique_positions, index = np.unique(events['transcriptomic_position'][idx_sorted],return_index = True)
    y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
#    read_id_arrays = np.split(events['read_id'][idx_sorted], index[1:])
    reference_kmer_arrays = np.split(events['reference_kmer'][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
    asserted = True
#     for key_tuple,y_array,reference_kmer_array in zip(key_tuples,y_arrays,reference_kmer_arrays):
    for position,y_array,reference_kmer_array in zip(unique_positions,y_arrays,reference_kmer_arrays):
        
        position = int(position)
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(y_array) == 0):
            continue
            
        if 'XXXXX' in set(reference_kmer_array):
            y_array = y_array[reference_kmer_array != 'XXXXX']  
            assert len(y_array) == len(reference_kmer_array) - (reference_kmer_array=='XXXXX').sum()
            reference_kmer_array = reference_kmer_array[reference_kmer_array != 'XXXXX']  
            
        try:
            assert len(set(reference_kmer_array)) == 1
        except:
            asserted = False
            break
        kmer = set(reference_kmer_array).pop()

        data[position] = {kmer: list(np.around(y_array,decimals=2))}
        
    # write to file.
    writeOutputs_xpore(tx_id,asserted,data,data_dict,locks,out_paths)

def preprocess_tx_m6anet(tx_id,data_dict,n_neighbors,out_paths,locks):  # todo
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

    # Concatenate
    if len(data_dict) == 0:
        return
    
    features_arrays = []
    reference_kmer_arrays = []
    transcriptomic_positions_arrays = []

    for read_id,events_per_read in data_dict.items(): 
        # print(read_id)
        events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS)
        for event_per_read in events_per_read:
            features_arrays.append(event_per_read[0])
            reference_kmer_arrays.append([combine_sequence(kmer) for kmer in event_per_read[1]])
            transcriptomic_positions_arrays.append(event_per_read[3])

    if len(features_arrays) == 0:
        return
    else:
        features_arrays = np.concatenate(features_arrays)
        reference_kmer_arrays = np.concatenate(reference_kmer_arrays)
        transcriptomic_positions_arrays = np.concatenate(transcriptomic_positions_arrays)
        assert(len(features_arrays) == len(reference_kmer_arrays) == len(transcriptomic_positions_arrays))
    # Sort and split

    idx_sorted = np.argsort(transcriptomic_positions_arrays)
    positions, index = np.unique(transcriptomic_positions_arrays[idx_sorted], return_index = True,axis=0) #'chr',
    features_arrays = np.split(features_arrays[idx_sorted], index[1:])
    reference_kmer_arrays = np.split(reference_kmer_arrays[idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)


    # for each position, make it ready for json dump
    for position, features_array, reference_kmer_array in zip(positions, features_arrays, reference_kmer_arrays):
        kmer = set(reference_kmer_array)
        assert(len(kmer) == 1)
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(features_array) == 0):
            continue

        data[int(position)] = {kmer.pop(): features_array.tolist()}

    # write to file.
    writeOutputs_m6anet(tx_id,data,locks,out_paths)     

def preprocess_tx_xpore_m6anet(tx_id,data_dict,n_neighbors,out_paths,locks):  # todo
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

    # Concatenate
    if len(data_dict) == 0:
        return
    ##for xpore
    xpore_events = []
    ##for m6anet
    m6anet_features_arrays = []
    m6anet_reference_kmer_arrays = [] 
    m6anet_transcriptomic_positions_arrays = []

    for read_id,events_per_read in data_dict.items(): 
        ##for xpore
        xpore_events += [events_per_read]
        ##for m6anet
        events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS)
        for event_per_read in events_per_read:
            m6anet_features_arrays.append(event_per_read[0])
            m6anet_reference_kmer_arrays.append([combine_sequence(kmer) for kmer in event_per_read[1]])
            m6anet_transcriptomic_positions_arrays.append(event_per_read[3])
    ##for xpore
    xpore_events = np.concatenate(xpore_events)
    ##for m6anet
    if len(m6anet_features_arrays) == 0:
        return
    else:
        m6anet_features_arrays = np.concatenate(m6anet_features_arrays)
        m6anet_reference_kmer_arrays = np.concatenate(m6anet_reference_kmer_arrays)
        m6anet_transcriptomic_positions_arrays = np.concatenate(m6anet_transcriptomic_positions_arrays)
        assert(len(m6anet_features_arrays) == len(m6anet_reference_kmer_arrays) == len(m6anet_transcriptomic_positions_arrays))

    ##for xpore # Sort and split 
    xpore_idx_sorted = np.argsort(xpore_events['transcriptomic_position'])
    xpore_unique_positions, xpore_index = np.unique(xpore_events['transcriptomic_position'][xpore_idx_sorted],return_index = True)
    xpore_y_arrays = np.split(xpore_events['norm_mean'][xpore_idx_sorted], xpore_index[1:])
    xpore_reference_kmer_arrays = np.split(xpore_events['reference_kmer'][xpore_idx_sorted], xpore_index[1:])

    ##for m6anet # Sort and split
    m6anet_idx_sorted = np.argsort(m6anet_transcriptomic_positions_arrays)
    m6anet_positions, m6anet_index = np.unique(m6anet_transcriptomic_positions_arrays[m6anet_idx_sorted], return_index = True,axis=0) #'chr',
    m6anet_features_arrays = np.split(m6anet_features_arrays[m6anet_idx_sorted], m6anet_index[1:])
    m6anet_reference_kmer_arrays = np.split(m6anet_reference_kmer_arrays[m6anet_idx_sorted], m6anet_index[1:])

    ##for xpore # Prepare
    xpore_data = defaultdict(dict)
    # for each position, make it ready for json dump
    asserted = True
    for position,y_array,reference_kmer_array in zip(xpore_unique_positions,xpore_y_arrays,xpore_reference_kmer_arrays):
        position = int(position)
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(y_array) == 0):
            continue
        if 'XXXXX' in set(reference_kmer_array):
            y_array = y_array[reference_kmer_array != 'XXXXX']  
            assert len(y_array) == len(reference_kmer_array) - (reference_kmer_array=='XXXXX').sum()
            reference_kmer_array = reference_kmer_array[reference_kmer_array != 'XXXXX']  
        try:
            assert len(set(reference_kmer_array)) == 1
        except:
            asserted = False
            break
        kmer = set(reference_kmer_array).pop()
        xpore_data[position] = {kmer: list(np.around(y_array,decimals=2))}

    ##for m6anet # Prepare
    m6anet_data = defaultdict(dict)
    # for each position, make it ready for json dump
    for position, features_array, reference_kmer_array in zip(m6anet_positions, m6anet_features_arrays, m6anet_reference_kmer_arrays):
        kmer = set(reference_kmer_array)
        assert(len(kmer) == 1)
        if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (len(features_array) == 0):
            continue
        m6anet_data[int(position)] = {kmer.pop(): features_array.tolist()}

    ##for xpore # write to file.
    writeOutputs_xpore(tx_id,asserted,xpore_data,data_dict,locks,out_paths)

    ##for m6anet # write to file.
    writeOutputs_m6anet(tx_id,m6anet_data,locks,out_paths)

# def index_nanopolish(eventalign_filepath,summary_filepath,out_dir,n_processes):
#     with helper.EventalignFile(eventalign_filepath) as eventalign_file, open(summary_filepath,'r') as summary_file:
        
#         reader_summary = csv.DictReader(summary_file, delimiter="\t")
#         reader_eventalign = csv.DictReader(eventalign_file, delimiter="\t")

#         row_summary = next(reader_summary)
#         read_name = row_summary['read_name']
#         read_index = row_summary['read_index']
#         eventalign_per_read = []
#         for row_eventalign in reader_eventalign:
#             if (row_eventalign['read_index'] == read_index):
#                 eventalign_per_read += [row_eventalign]
#             else: 
#                 # Load a read info to the task queue.
#                 if read_name not in read_names_done:
#                     task_queue.put((read_name,eventalign_per_read,out_paths))
#                 # Next read.
#                 try:
#                     row_summary = next(reader_summary)
#                 except StopIteration: # no more read.
#                     break
#                 else:
#                     read_index = row_summary['read_index']
#                     read_name = row_summary['read_name']
#                     assert row_eventalign['read_index'] == read_index 
#                     eventalign_per_read = [row_eventalign]



def main():
    args = get_args()
    #
    n_processes = args.n_processes        
    eventalign_filepath = args.eventalign
    summary_filepath = args.summary
    chunk_size = args.chunk_size
    out_dir = args.out_dir
    readcount_min = args.readcount_min    
    readcount_max = args.readcount_max
    resume = args.resume
    genome = args.genome
    program = args.program
    n_neighbors = args.n_neighbors

    program = program.split(",")
    if len(program) == 2:
        fun_dict={'combine':combine_m6anet,'preprocess_gene':preprocess_gene_xpore_m6anet,'preprocess_tx':preprocess_tx_xpore_m6anet}
    else:
        if program[0] == 'xpore':
            fun_dict={'combine':combine_xpore,'preprocess_gene':preprocess_gene_xpore,'preprocess_tx':preprocess_tx_xpore}
        else:
            fun_dict={'combine':combine_m6anet,'preprocess_gene':preprocess_gene_m6anet,'preprocess_tx':preprocess_tx_m6anet} 

    if genome and (None in [args.gtf_path_or_url,args.transcript_fasta_paths_or_urls]):
        print('please provide the following')
        print('- gtf_path_or_url')
        print('- transcript_fasta_paths_or_urls')
    else:
        gtf_path_or_url = args.gtf_path_or_url
        transcript_fasta_paths_or_urls = args.transcript_fasta_paths_or_urls
        
    misc.makedirs(out_dir) #todo: check every level.
    
    # (1) For each read, combine multiple events aligned to the same positions, the results from nanopolish eventalign, into a single event per position.
    if not args.skip_eventalign_indexing:
        parallel_index(eventalign_filepath,summary_filepath,chunk_size,out_dir,n_processes,resume)
    
    # (2) Create a .json file, where the info of all reads are stored per position, for modelling.
    if genome:
        if os.path.exists(transcript_fasta_paths_or_urls+'.pickle') and os.path.exists(gtf_path_or_url+'.pickle'):
            with open(transcript_fasta_paths_or_urls+'.pickle','rb') as fasta_pickle, open(gtf_path_or_url+'.pickle','rb') as gtf_pickle:
                fasta_dict = pickle.load(fasta_pickle)
                gtf_dict = pickle.load(gtf_pickle)
        else:
            fasta_dict = readFasta(transcript_fasta_paths_or_urls)
            gtf_dict = readGTF(gtf_path_or_url)
        parallel_preprocess_gene(program,eventalign_filepath,fasta_dict,gtf_dict,fun_dict,out_dir,n_processes,readcount_min,readcount_max,n_neighbors,resume)

    else:
        parallel_preprocess_tx(program,eventalign_filepath,fun_dict,out_dir,n_processes,readcount_min,readcount_max,n_neighbors,resume)

if __name__ == '__main__':
    main()
