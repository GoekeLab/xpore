import numpy as np
import pandas as pd
import os,re
import multiprocessing 
import h5py
import csv
import ujson
from operator import itemgetter
from collections import defaultdict
from io import StringIO

from . import helper
from ..utils import misc

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

def parallel_index(eventalign_filepath,chunk_size,out_dir,n_processes,resume):
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
    
def t2g(gene_id,fasta_dict,annotation_dict,g2t_mapping,df_eventalign_index,readcount_min):
    tx_ids = []
    t2g_dict = {}
    transcripts = [tx for tx in annotation_dict if tx in g2t_mapping[gene_id]]
    n_reads = sum([len(df_eventalign_index.loc[tx]) for tx in transcripts])
    if n_reads >= readcount_min:
        for tx in transcripts:
            tx_seq = fasta_dict[tx][0]
            tx_contig = annotation_dict[tx]['chr']
            if tx_seq is None:
                continue
            for exon_num in range(len(annotation_dict[tx]['exon'])):
                g_interval=annotation_dict[tx]['exon'][exon_num]
                tx_interval=annotation_dict[tx]['tx_exon'][exon_num]
                for g_pos in range(g_interval[0],g_interval[1]+1): # Exclude the rims of exons.
                    dis_from_start = g_pos - g_interval[0]
                    if annotation_dict[tx]['strand'] == "+":
                        tx_pos = tx_interval[0] + dis_from_start
                    elif annotation_dict[tx]['strand'] == "-":
                        tx_pos = tx_interval[1] - dis_from_start
                    if (g_interval[0] <= g_pos < g_interval[0]+2) or (g_interval[1]-2 < g_pos <= g_interval[1]): # Todo: To improve the mapping
                        kmer = 'XXXXX'
                    else:
                        kmer = tx_seq[tx_pos-2:tx_pos+3]
                    t2g_dict[(tx,tx_pos)] = (tx_contig,gene_id,g_pos,kmer) # tx.contig is chromosome.
            tx_ids += [tx]

                
    return n_reads, tx_ids, t2g_dict

def combine(events_str,kmer_source='reference_kmer'):
    f_string = StringIO(events_str)
    eventalign_result = pd.read_csv(f_string,delimiter='\t',names=['contig','position','reference_kmer','read_index',
                         'strand','event_index','event_level_mean','event_stdv','event_length','model_kmer',
                         'model_mean', 'model_stdv', 'standardized_level', 'start_idx', 'end_idx'])
    f_string.close()

    # Reverse complement helper for matching genome-aligned reads where
    # reference_kmer may be the reverse complement of model_kmer
    def _revcomp(seq):
        comp = str.maketrans('ACGT','TGCA')
        return seq.translate(comp)[::-1]

    # When using model_kmer, accept both direct and reverse-complement matches
    # (genome alignments can have reverse-oriented reads).
    # When using reference_kmer, only accept exact matches (same as xPore v2.1 behaviour)
    if kmer_source == 'model_kmer':
        cond_successfully_eventaligned = (
            (eventalign_result['reference_kmer'] == eventalign_result['model_kmer']) |
            (eventalign_result['reference_kmer'].apply(_revcomp) == eventalign_result['model_kmer'])
        )
    else:
        cond_successfully_eventaligned = (
            eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
        )
    if cond_successfully_eventaligned.sum() != 0:

        eventalign_result = eventalign_result[cond_successfully_eventaligned]

        # Select kmer column based on kmer_source flag: model_kmer for genome alignments,
        # reference_kmer for transcriptome alignments
        kmer_col = 'model_kmer' if kmer_source == 'model_kmer' else 'reference_kmer'
        keys = ['read_index','contig','position',kmer_col] # for groupby
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


#        eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]    #### CHANGE MADE ####
        eventalign_result['transcript_id'] = [contig for contig in eventalign_result['contig']]
        #eventalign_result['transcript_id'] = eventalign_result['contig']

        eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
        # eventalign_result = misc.str_encode(eventalign_result)
#         eventalign_result['read_id'] = [read_name]*len(eventalign_result)

        # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx']
        # features_dtype = np.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'), ('end_idx', '<i8')])
        
#         features = ['transcript_id','transcriptomic_position','reference_kmer','norm_mean']

#         df_events = eventalign_result[['read_index']+features]
#         # print(df_events.head())

        # Use the selected kmer column in the output features
        features = ['transcript_id','transcriptomic_position',kmer_col,'norm_mean']
#        np_events = eventalign_result[features].reset_index().values.ravel().view(dtype=[('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')])
        df_events = eventalign_result[features]
        np_events = np.rec.fromrecords(df_events, names=[*df_events])
        return np_events, kmer_col

def readFasta(transcript_fasta,is_gff):
    fasta=open(transcript_fasta,"r")
    entries,separate_by_pipe="",False
    for ln in fasta:
        entries+=ln
    entries=entries.split(">")
    if len(entries[1].split("|"))>1:
        separate_by_pipe=True
    dict={}
    for entry in entries:
        entry=entry.split("\n")
        if len(entry[0].split()) > 0:
            id=entry[0].split('.')[0]
            seq="".join(entry[1:])
            dict[id]=[seq]
            if is_gff > 0:
                if separate_by_pipe == True:
                    g_id=info[1],split(".")[0]
                else:
                    g_id=entry[0].split("gene:")[1].split(".")[0]
                dict[id].append(g_id)
    return dict

def readAnnotation(gtf_or_gff):
    gtf=open(gtf_or_gff,"r")
    dict,is_gff={},0
    for ln in gtf:
        if not ln.startswith("#"):
            ln=ln.strip("\n").split("\t")
            if is_gff == 0:
                if ln[-1].startswith("ID") or ln[-1].startswith("Parent"):
                    is_gff = 1
                else:
                    is_gff = -1
            if is_gff < 0:
                if ln[2] == "transcript" or ln[2] == "exon":
                    chr,type,start,end=ln[0],ln[2],int(ln[3]),int(ln[4])
                    attrList=ln[-1].split(";")
                    attrDict={}
                    for k in attrList:
                        p=k.strip().split(" ")
                        if len(p) == 2:
                            attrDict[p[0]]=p[1].strip('\"')
                    ##tx_id=ln[-1].split('; transcript_id "')[1].split('";')[0]
                    ##g_id=ln[-1].split('gene_id "')[1].split('";')[0]
                    tx_id = attrDict["transcript_id"]
                    g_id = attrDict["gene_id"]
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
            if is_gff > 0:
                if ln[2] == "exon" or ln[2] == "mRNA":
                    chr,type,start,end=ln[0],ln[2],int(ln[3]),int(ln[4])
                    tx_id=ln[-1].split('transcript:')[1].split(';')[0]
                    if ln[2] == "mRNA":
                        type="transcript"
                    if tx_id not in dict:
                        dict[tx_id]={'chr':chr,'strand':ln[6]}
                        if type == "transcript":
                            dict[tx_id][type]=(start,end)
                        if type == "exon":
                            dict[tx_id][type]=[(start,end)]
                    else:
                        if type == "transcript" and type not in dict[tx_id]:
                            dict[tx_id][type]=(start,end)
                        if type == "exon":
                            if type not in dict[tx_id]:
                                dict[tx_id][type]=[(start,end)]
                            else:
                                dict[tx_id][type].append((start,end))
    #convert genomic positions to tx positions
    if is_gff < 0:
        for id in dict:
            tx_pos,tx_start=[],0
            for pair in dict[id]["exon"]:
                tx_end=pair[1]-pair[0]+tx_start
                tx_pos.append((tx_start,tx_end))
                tx_start=tx_end+1
            dict[id]['tx_exon']=tx_pos
    else:
        for id in dict:
            tx_pos,tx_start=[],0
            if dict[id]["strand"] == "-":
                dict[id]["exon"].sort(key=lambda tup: tup[0], reverse=True)
            for pair in dict[id]["exon"]:
                tx_end=pair[1]-pair[0]+tx_start
                tx_pos.append((tx_start,tx_end))
                tx_start=tx_end+1
            dict[id]['tx_exon']=tx_pos
    return (dict,is_gff)

def parallel_preprocess_gene(eventalign_filepath,fasta_dict,annotation_dict,is_gff,out_dir,n_processes,readcount_min,readcount_max,resume,kmer_source='reference_kmer'):
    
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

    # Get all gene ids and create a dict of eventalign.combine index.
#     gene_ids = set()

    if is_gff > 0: ##add g_id from fasta dict entry if gff annotation is used
        for tx_id in annotation_dict:
            try:
                annotation_dict[tx_id]['g_id']=fasta_dict[tx_id][1]
            except KeyError:
                continue

    df_eventalign_index = pd.read_csv(os.path.join(out_dir,'eventalign.index'))
    df_eventalign_index['transcript_id'] = [tx_id.split('.')[0] for tx_id in  df_eventalign_index['transcript_id']]
#    df_eventalign_index['transcript_id'] = [tx_id for tx_id in  df_eventalign_index['transcript_id']]
    df_eventalign_index.set_index('transcript_id',inplace=True)
    g2t_mapping = defaultdict(list)

    for tx_id in set(df_eventalign_index.index):
        try:
##           g_id = ensembl.transcript_by_id(tx_id).gene_id 
            g_id = annotation_dict[tx_id]['g_id'] 
        except KeyError:
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

    # kmer_col == kmer_source (argparse restricts it to a valid column name); bind it up front so it
    # is always defined even if no read in this gene returns a result from combine() (avoids NameError).
    kmer_col = kmer_source

    with open(eventalign_filepath,'r') as eventalign_result:

        for gene_id in g2t_mapping:
                        
            if resume and (gene_id in gene_ids_done):
                continue
            # mapping a gene <-> transcripts

            n_reads, tx_ids, t2g_mapping = t2g(gene_id,fasta_dict,annotation_dict,g2t_mapping,df_eventalign_index,readcount_min)
            #
            if n_reads >= readcount_min: 
                data_dict = dict()
                readcount = 0
                for tx_id in tx_ids:
                    for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                        read_index,pos_start,pos_end = row['read_index'],row['pos_start'],row['pos_end']
                        eventalign_result.seek(pos_start,0)
                        events_str = eventalign_result.read(pos_end-pos_start)
                        # Pass kmer_source to combine so it selects the correct kmer column
                        result = combine(events_str,kmer_source)
                        #data = np.genfromtxt(f_string,delimiter=',',dtype=np.dtype([('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')]))
                        if (result is not None):
                            data, kmer_col = result
                            if data.size > 1:
                                data_dict[read_index] = data
                        readcount += 1

                if len(data_dict)>=readcount_min:
#                     print(gene_id,len(data_dict)) #len(data_dict) is the number of reads to be processed.
                    task_queue.put((gene_id,data_dict,kmer_col,t2g_mapping,readcount_max,out_paths)) # Blocked if necessary until a free slot is available.
                    gene_ids_processed += [gene_id]


    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(out_paths['log'],'a+') as f:
        f.write('Total %d genes.\n' %len(gene_ids_processed))
        f.write(helper.decor_message('successfully finished'))

def preprocess_gene(gene_id,data_dict,kmer_col,t2g_mapping,readcount_max,out_paths,locks):
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
    condition_labels = []
    run_labels = []
    read_ids = []
    genomic_coordinates = []
    
    # Concatenate
#     if len(data_dict) == 0:
#         return


    for read_index,events_per_read in data_dict.items():
#         if len(events_per_read) > 0:
        # ===== transcript to gene coordinates ===== # TODO: to use gtf.
#        tx_ids = [tx_id.decode('UTF-8').split('.')[0] for tx_id in events_per_read['transcript_id']]
        tx_ids = [tx_id.split('.')[0] for tx_id in events_per_read['transcript_id']] 
        tx_positions = events_per_read['transcriptomic_position']
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- np structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = np.array(genomic_coordinate,dtype=np.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        # ===== 

        # Based on Ensembl, remove transcript version.

        events_per_read['transcript_id'] = tx_ids
        # Use dynamic kmer_col (reference_kmer or model_kmer) in the dtype based on kmer_source
        events_per_read = np.array(events_per_read,dtype=np.dtype([('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), (kmer_col, 'S5'), ('norm_mean', '<f8')]))

        #

        events += [events_per_read]
        genomic_coordinates += [genomic_coordinate]
        n_events_per_read = len(events_per_read)
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
        # Cap reads per site (not per gene) to readcount_max; None means no limit
        if readcount_max is not None and len(y_array) > readcount_max:
            y_array = y_array[:readcount_max]
            g_kmer_array = g_kmer_array[:readcount_max]
            g_positions_array = g_positions_array[:readcount_max]
        if (len(set(g_kmer_array)) == 1) and ('XXXXX' in set(g_kmer_array)) or (len(y_array) == 0):
            continue

        if 'XXXXX' in set(g_kmer_array):
            y_array = y_array[g_kmer_array != 'XXXXX']
            assert len(y_array) == len(g_kmer_array) - (g_kmer_array=='XXXXX').sum()
            g_kmer_array = g_kmer_array[g_kmer_array != 'XXXXX']

        try:
            assert len(set(g_kmer_array)) == 1
            assert list(set(g_kmer_array))[0].count('N') == 0 ##to weed out the mapped kmers from tx_seq that contain 'N', which is not in diffmod's model_kmer
            assert {position} == set(g_positions_array)
        except:
            asserted = False
            continue  # skip this position, don't stop processing remaining positions
        kmer = set(g_kmer_array).pop()

        data[position] = {kmer: list(y_array)} #,'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}
        
    # write to file.
    log_str = '%s: %s' %(gene_id,asserted)

    with locks['json'], open(out_paths['json'],'a') as f:

        pos_start = f.tell()
        f.write('{')
        f.write('"%s":' %gene_id)
        ujson.dump(data, f)
        f.write('}\n')
        pos_end = f.tell()

    with locks['index'], open(out_paths['index'],'a') as f:
        f.write('%s,%d,%d\n' %(gene_id,pos_start,pos_end))
        
    with locks['readcount'], open(out_paths['readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
        n_reads = len(data_dict)
        f.write('%s,%d\n' %(gene_id,n_reads))
        
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(log_str + '\n')


def parallel_preprocess_tx(eventalign_filepath,out_dir,n_processes,readcount_min,readcount_max,resume,kmer_source='reference_kmer'):
    
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
    df_eventalign_index = pd.read_csv(os.path.join(out_dir,'eventalign.index'))
#    df_eventalign_index['transcript_id'] = [tx_id.split('.')[0] for tx_id in  df_eventalign_index['transcript_id']]
#    df_eventalign_index['transcript_id'] = [tx_id for tx_id in  df_eventalign_index['transcript_id']]
    tx_ids = df_eventalign_index['transcript_id'].values.tolist()
    tx_ids = list(dict.fromkeys(tx_ids))
    df_eventalign_index.set_index('transcript_id',inplace=True)
    # kmer_col == kmer_source (argparse restricts it to a valid column name); bind it up front so it
    # is always defined even if no read in this transcript returns a result from combine() (avoids NameError).
    kmer_col = kmer_source
    with open(eventalign_filepath,'r') as eventalign_result:
        for tx_id in tx_ids:
            data_dict = dict()
            readcount = 0
            for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                read_index,pos_start,pos_end = row['read_index'],row['pos_start'],row['pos_end']
                eventalign_result.seek(pos_start,0)
                events_str = eventalign_result.read(pos_end-pos_start)
                # Pass kmer_source to combine so it selects the correct kmer column
                result = combine(events_str,kmer_source)
                if (result is not None):
                    data, kmer_col = result
                    if data.size > 1:
                        data_dict[read_index] = data
                readcount += 1
            if readcount>=readcount_min:
                task_queue.put((tx_id,data_dict,kmer_col,readcount_max,out_paths)) # Blocked if necessary until a free slot is available.
                tx_ids_processed += [tx_id]

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(out_paths['log'],'a+') as f:
        f.write('Total %d transcripts.\n' %len(tx_ids_processed))
        f.write(helper.decor_message('successfully finished'))

def preprocess_tx(tx_id,data_dict,kmer_col,readcount_max,out_paths,locks):
    """
    Reshape one aligned reference's per-read events into per-position signal arrays.

    Transcriptome mode (no --genome) (or genome alignment mode): positions stay in original reference coordinates;
    no conversion to genomic coordinates for tx alignments is done here (that's preprocess_gene). Groups reads at each
    position, caps them at readcount_max, checks k-mer consistency, and appends the
    result to the data.json / .index / .readcount / .log files. Returns None.

    Parameters
    ----------
    tx_id : str                     Reference (for tx alignments, Transcript ID).
    data_dict : {read_index: events_array}   Per-read arrays from combine().
    kmer_col : str                  'reference_kmer' or 'model_kmer'.
    readcount_max : int or None     Per-site read cap; None means no cap.
    out_paths, locks : dict         Output file paths and their write locks.
    
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
    if kmer_col == 'model_kmer':
        # Group by (position, kmer) so forward/reverse reads at the same position
        # produce separate entries with their respective kmers
        idx_sorted = np.lexsort((events[kmer_col], events['transcriptomic_position']))
        positions_sorted = events['transcriptomic_position'][idx_sorted]
        kmers_sorted = events[kmer_col][idx_sorted]
        compound_keys = np.array(list(zip(positions_sorted, kmers_sorted)),
                                 dtype=[('pos', positions_sorted.dtype), ('kmer', kmers_sorted.dtype)])
        unique_pairs, index = np.unique(compound_keys, return_index=True)
        unique_positions = unique_pairs['pos']
    else:
        idx_sorted = np.argsort(events['transcriptomic_position'])
        unique_positions, index = np.unique(events['transcriptomic_position'][idx_sorted], return_index=True)
    y_arrays = np.split(events['norm_mean'][idx_sorted], index[1:])
    kmer_arrays = np.split(events[kmer_col][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
    asserted = True
#     for key_tuple,y_array,kmer_array in zip(key_tuples,y_arrays,kmer_arrays):
    for position,y_array,kmer_array in zip(unique_positions,y_arrays,kmer_arrays):

        position = int(position)
        # Cap reads per site (not per transcript) to readcount_max; None means no limit
        if readcount_max is not None and len(y_array) > readcount_max:
            y_array = y_array[:readcount_max]
            kmer_array = kmer_array[:readcount_max]
        if (len(set(kmer_array)) == 1) and ('XXXXX' in set(kmer_array)) or (len(y_array) == 0):
            continue

        if 'XXXXX' in set(kmer_array):
            y_array = y_array[kmer_array != 'XXXXX']
            assert len(y_array) == len(kmer_array) - (kmer_array=='XXXXX').sum()
            kmer_array = kmer_array[kmer_array != 'XXXXX']

        try:
            assert len(set(kmer_array)) == 1
            assert list(set(kmer_array))[0].count('N') == 0 ##to weed out the mapped kmers from tx_seq that contain 'N', which is not in diffmod's model_kmer
        except:
            asserted = False
            continue  # skip this position, don't stop processing remaining positions
        kmer = set(kmer_array).pop()

        data[position][kmer] = list(np.around(y_array,decimals=2))
        
    # write to file.
    log_str = '%s: %s.' %(tx_id,asserted)
    with locks['json'], open(out_paths['json'],'a') as f:
        pos_start = f.tell()
        f.write('{')
        f.write('"%s":' %tx_id)
        ujson.dump(data, f)
        f.write('}\n')
        pos_end = f.tell()
        
    with locks['index'], open(out_paths['index'],'a') as f:
        f.write('%s,%d,%d\n' %(tx_id,pos_start,pos_end))
        
    with locks['readcount'], open(out_paths['readcount'],'a') as f: #todo: repeats no. of tx >> don't want it.
        n_reads = len(data_dict)
        f.write('%s,%d\n' %(tx_id,n_reads))
        
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(log_str + '\n')
        

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

#def check_gene_tx_id_version(gtf_or_gff):
#    gtf=open(gtf_or_gff,"r")
#    extra_version_fields=0
#    for i in range(25):
#        ln=gtf.readline().split('\t')
#        if not ln[0].startswith('#'):
#            if ln[2] == "transcript" or ln[2] == "exon":
#                check_transcript_version = len(ln[-1].split('transcript_version')) == 2
#                check_gene_version = len(ln[-1].split('gene_version')) == 2
#                if check_transcript_version and check_gene_version:
#                   extra_version_fields+=1
#    if extra_version_fields > 0:
#        return True
#    else:
#        return False
 
#def mergeGTFtxIDversion(gtf_or_gff,out_dir):
#    gtf=open(gtf_or_gff,"r")
#    new_gtf_path=os.path.join(out_dir,'transcript_id_version_merged.gtf')
#    new_gtf=open(new_gtf_path,"w")
#    for ln in gtf:
#        if not ln.startswith("#"):
#            ln=ln.split("\t")
#            if ln[2] == "transcript" or ln[2] == "exon":
#                g_id=ln[-1].split('gene_id "')[1].split('";')[0]
#                g_ver=ln[-1].split('; gene_version "')[1].split('";')[0]
#                tx_id=ln[-1].split('; transcript_id "')[1].split('";')[0]
#                tx_ver=ln[-1].split('; transcript_version "')[1].split('";')[0]
#                new_gtf.write('\t'.join(ln[:-1])+'\t'+'gene_id "'+g_id+'.'+g_ver+'"; transcript_id "'+tx_id+'.'+tx_ver+'";'+'\n')
#    new_gtf.close()
#    return new_gtf_path

def dataprep(args):
    #
    n_processes = args.n_processes        
    eventalign_filepath = args.eventalign
    chunk_size = args.chunk_size
    out_dir = args.out_dir
    readcount_min = args.readcount_min    
    readcount_max = args.readcount_max
    resume = args.resume
    genome = args.genome
    kmer_source = args.kmer_source

    # Check for likely genome alignment when user has selected reference_kmer
    if kmer_source == 'reference_kmer':
        with open(eventalign_filepath,'r') as f:
            first_line = f.readline()
            if first_line:
                contig = first_line.split('\t')[0]
                if contig.startswith('chr'):
                    raise ValueError(
                        "It looks like you are using a genome alignment (reference name '%s' starts with 'chr'), "
                        "but --kmer_source is set to 'reference_kmer'. Genome alignments can contain reverse-oriented "
                        "reads, using reference_kmer will cause all reverse reads to be filtered. Please use '--kmer_source model_kmer' instead." % contig
                    )

    if genome and (None in [args.gtf_or_gff,args.transcript_fasta]):
        print('please provide the following')
        print('- gtf_or_gff')
        print('- transcript_fasta')
    else:
        gtf_or_gff = args.gtf_or_gff
        transcript_fasta = args.transcript_fasta
        
    misc.makedirs(out_dir) #todo: check every level.
    
    # (1) For each read, combine multiple events aligned to the same positions, the results from nanopolish eventalign, into a single event per position.
    if not args.skip_eventalign_indexing:
        parallel_index(eventalign_filepath,chunk_size,out_dir,n_processes,resume)
    
    # (2) Create a .json file, where the info of all reads are stored per position, for modelling.
    if genome:
#        merge_transcript_id_version = check_gene_tx_id_version(gtf_or_gff)
#        if merge_transcript_id_version:
#            gtf_or_gff = mergeGTFtxIDversion(gtf_or_gff,out_dir)
        annotation_dict,is_gff = readAnnotation(gtf_or_gff)
        fasta_dict = readFasta(transcript_fasta,is_gff)
        parallel_preprocess_gene(eventalign_filepath,fasta_dict,annotation_dict,is_gff,out_dir,n_processes,readcount_min,readcount_max,resume,kmer_source)
    else:
        parallel_preprocess_tx(eventalign_filepath,out_dir,n_processes,readcount_min,readcount_max,resume,kmer_source)
