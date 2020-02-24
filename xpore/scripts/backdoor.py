import argparse
import numpy
import pandas
import os
import multiprocessing 
import h5py
import csv
import json
import pysam #0-based leftmost coordinate
from pyensembl import EnsemblRelease
from tqdm import tqdm
from operator import itemgetter
from collections import defaultdict
import numpy.lib.recfunctions as rfn 

from . import helper
from ..utils import misc

def get_args():
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('--eventalign', dest='eventalign', help='Eventalign filepath from nanopolish.',required=True)
    parser.add_argument('--summary', dest='summary', help='Summary filepath from nanopolish.',required=True)
    parser.add_argument('--bamtx', dest='bamtx', help='bamtx filepath.',required=True)
    parser.add_argument('--mapping', dest='mapping', help='gene-transcript mapping directory.',required=True)
    parser.add_argument('--data_dir', dest='data_dir', help='Data directory.',required=True)
    parser.add_argument('--out_dir', dest='out_dir', help='Output directory.',required=True)


    # Optional
    # parser.add_argument('--features', dest='features', help='Signal features to extract.',type=list,default=['norm_mean'])
    parser.add_argument('--n_processes', dest='n_processes', help='Number of processes.',type=int, default=1)
    parser.add_argument('--ensembl', dest='ensembl', help='ensembl version.',type=int, default=91)
    parser.add_argument('--read_count_min', dest='read_count_min', help='Minimum of read counts per gene.',type=int, default=1)
    parser.add_argument('--resume', dest='resume', help='Resume.',default=False,action='store_true') #todo


    
    return parser.parse_args()

def combine(read_name,eventalign_per_read,out_paths,read_ids_processed,locks):
    if read_ids_processed[read_name]:
        return

    eventalign_result = pandas.DataFrame.from_records(eventalign_per_read)

    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    eventalign_result = eventalign_result[cond_successfully_eventaligned]

    keys = ['read_index','contig','position','reference_kmer'] # for groupby
    eventalign_result['length'] = pandas.to_numeric(eventalign_result['end_idx'])-pandas.to_numeric(eventalign_result['start_idx'])
    eventalign_result['sum_norm_mean'] = pandas.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']

    eventalign_result = eventalign_result.groupby(keys)  
    sum_norm_mean = eventalign_result['sum_norm_mean'].sum() 
    start_idx = eventalign_result['start_idx'].min()
    end_idx = eventalign_result['end_idx'].max()
    total_length = eventalign_result['length'].sum()

    eventalign_result = pandas.concat([start_idx,end_idx],axis=1)
    eventalign_result['norm_mean'] = sum_norm_mean/total_length

    eventalign_result.reset_index(inplace=True)

    eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]
    eventalign_result['transcriptomic_position'] = pandas.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
    # eventalign_result = misc.str_encode(eventalign_result)
    eventalign_result['read_id'] = [read_name]*len(eventalign_result)

    # features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean','start_idx','end_idx']
    # features_dtype = numpy.dtype([('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8'), ('start_idx', '<i8'), ('end_idx', '<i8')])
    features = ['read_id','transcript_id','transcriptomic_position','reference_kmer','norm_mean']

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


def parallel_combine(eventalign_filepath,summary_filepath,out_dir,n_processes,resume):
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['hdf5','log']:
        out_paths[out_filetype] = os.path.join(out_dir,'eventalign.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
        
        
    read_ids_processed = defaultdict(lambda: False, read_ids_processed)
    if resume:
        with open(out_paths['log'],'r') as f:
            reader = csv.reader(f)
            read_ids_processed.update({read_id: True for read_id in sum(list(reader),[])})
            
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
    with open(eventalign_filepath,'r') as eventalign_file, open(summary_filepath,'r') as summary_file:

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
                task_queue.put((read_name,eventalign_per_read,out_paths,read_ids_processed))
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

def count_reads(version,bamtx_filepath,out_dir):
    """
    Counts number of aligned reads from bamtx file. Returns a dataframe of ['chr','gene_id','gene_name','transcript_id','n_reads'].
    """
    grch38 = EnsemblRelease(version) # human reference genome GRCh38 release 91 used in the ont mapping.    
    bam_file = pysam.AlignmentFile(bamtx_filepath, "rb")

    rows = []
    for contig in tqdm(bam_file.references):
        n_reads = bam_file.count(contig=contig,read_callback=lambda r: (not r.is_secondary) and (not r.is_qcfail) and ((r.flag==0) | (r.flag==16)) ) #flag: 0(+),16(-) 
        if n_reads == 0:
            continue
        tx_id = contig.split('.')[0]     
        try:
            tx = grch38.transcript_by_id(tx_id)
        except ValueError:
            continue
        else:
            g_id = tx.gene_id
            chromosome_id = tx.contig
            g_name = grch38.gene_name_of_gene_id(g_id)
            rows += [[chromosome_id,g_id,g_name,tx_id,n_reads]]

    df_count = pandas.DataFrame.from_records(rows,columns=['chr','gene_id','gene_name','transcript_id','n_reads'])
    
    # write to file.
    read_count_filepath = os.path.join(out_dir,'read_count.csv')
    df_count.to_csv(read_count_filepath,index=False,header=True,sep=',')
    return df_count
    
def parallel_preprocess(df_count,gt_mapping_dir,out_dir,n_processes,read_count_min):
    
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['json','index','log']:
        out_paths[out_filetype] = os.path.join(out_dir,'data.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
                
    # Writing the starting of the files.
    with open(out_paths['json'],'w') as f:
        f.write('{\n')
        f.write('"genes":{')
    with open(out_paths['index'],'w') as f:
        f.write('gene_id,start,end\n') # header
    open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=preprocess,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Load tasks into task_queue.
    df_count.set_index('gene_id',inplace=True)
    # gene_ids = ['ENSG00000168496','ENSG00000204388','ENSG00000123989','ENSG00000170144'] #todo
    # gene_ids = df_count.index
    
    # Load tasks into task_queue.
    with open('/ploy_ont_data/genes.csv','r') as f:
        reader = csv.reader(f)
        gene_ids = sum(list(reader),[]) # flatten the list.

    with h5py.File(os.path.join(out_dir,'eventalign.hdf5'),'r') as f:
        gene_ids_processed = []
        for gene_id in gene_ids:            
            gt_mapping_filepath = os.path.join(gt_mapping_dir,'%s.csv' %gene_id)
            df_gt = pandas.read_csv(gt_mapping_filepath)

            keys = zip(df_gt['tx_id'], df_gt['tx_pos'])
            values = zip(df_gt['chr'], df_gt['g_id'], df_gt['g_pos'], df_gt['kmer'])
            t2g_mapping = dict(zip(keys,values))
            
            try:
                df = df_count.loc[gene_id]
            except KeyError:
                continue
            n_reads = df['n_reads'].sum()
            read_ids = []
            if n_reads >= read_count_min:
                tx_ids = df_gt['tx_id'].unique()
                data_dict = dict()
                for tx_id in tx_ids:
                    if tx_id not in f: # no eventalign for tx_id
                        continue
                    for read_id in f[tx_id].keys():
                        if read_id not in read_ids:
                            data_dict[read_id] = f[tx_id][read_id]['events'][:]
                            read_ids += [read_id]

            if len(read_ids) > read_count_min:
                task_queue.put((gene_id,data_dict,t2g_mapping,out_paths)) # Blocked if necessary until a free slot is available. 
                gene_ids_processed += [gene_id]

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

    # Write the ending of the json file.
    with open(out_paths['json'],'a+') as f:
        f.seek(0,2)  # end of file
        f.truncate(f.tell()-1) 
        f.write('\n}\n}\n')
    ###
    
    with open(out_paths['log'],'a+') as f:
        f.write('Total %d genes.\n' %len(gene_ids_processed))
        f.write(helper.decor_message('successfully finished'))

def preprocess(gene_id,data_dict,t2g_mapping,out_paths,locks):  
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

        # ===== transcript to gene coordinates ===== #
        tx_ids = [tx_id.decode('UTF-8') for tx_id in events_per_read['transcript_id']]
        tx_positions = events_per_read['transcriptomic_position']
        
        genomic_coordinate = list(itemgetter(*zip(tx_ids,tx_positions))(t2g_mapping)) # genomic_coordinates -- numpy structured array of 'chr','gene_id','genomic_position','kmer'
        genomic_coordinate = numpy.array(genomic_coordinate,dtype=numpy.dtype([('chr','<U2'),('gene_id','<U15'),('genomic_position','<i4'),('g_kmer','<U5')]))
        # ===== 

        events += [events_per_read]
        genomic_coordinates += [genomic_coordinate]
        n_events_per_read = len(events_per_read)
        
    events = numpy.concatenate(events)
    genomic_coordinates = numpy.concatenate(genomic_coordinates)
   
    # Sort and split 
    idx_sorted = numpy.lexsort((events['reference_kmer'],genomic_coordinates['genomic_position'],genomic_coordinates['gene_id']))
    key_tuples, index = numpy.unique(list(zip(genomic_coordinates['gene_id'][idx_sorted],genomic_coordinates['genomic_position'][idx_sorted],events['reference_kmer'][idx_sorted])),return_index = True,axis=0) #'chr',
    y_arrays = numpy.split(events['norm_mean'][idx_sorted], index[1:])
    read_id_arrays = numpy.split(events['read_id'][idx_sorted], index[1:])
    g_kmer_arrays = numpy.split(genomic_coordinates['g_kmer'][idx_sorted], index[1:])

    # Prepare
    # print('Reformating the data for each genomic position ...')
    data = defaultdict(dict)
    # for each position, make it ready for json dump
    for key_tuple,y_array,read_id_array,g_kmer_array in zip(key_tuples,y_arrays,read_id_arrays,g_kmer_arrays):
        gene_id,position,kmer = key_tuple
        if (len(set(g_kmer_array)) == 1) and ('XXXXX' in set(g_kmer_array)) or (len(y_array) == 0):
            continue
                        
        data[position][kmer] = {'norm_means': list(y_array),
        'read_ids': [read_id.decode('UTF-8') for read_id in read_id_array]}
        
    # write to file.
    log_str = '%s: Data preparation ... Done.' %(gene_id)

    with locks['json'], open(out_paths['json'],'a') as f:
        f.write('\n')
        pos_start = f.tell()
        f.write('"%s":' %gene_id)
        json.dump(data, f)
        pos_end = f.tell()
        f.write(',')
    with locks['index'], open(out_paths['index'],'a') as f:
        f.write('%s,%d,%d\n' %(gene_id,pos_start,pos_end))
    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(log_str + '\n')
        
def reformat_eventalign(tx_id,tx_dir,out_paths,locks):
    columns = ['transcript_id','position','reference_kmer','norm_mean']


    for read_filename in os.listdir(tx_dir):        
        with h5py.File(os.path.join(tx_dir,read_filename),'r') as f:
            events = f['nanopolish']['eventalign']['events'][:][columns]
            read_id = read_filename.split('.')[0]
            read_ids = numpy.array([read_id]*len(events),dtype=numpy.dtype([('read_id','S36')]))
            events = rfn.merge_arrays((read_ids, events), asrecarray=True, flatten=True)  
            new_dtype=[('read_id', 'S36'), ('transcript_id', 'S15'), ('transcriptomic_position', '<i8'), ('reference_kmer', 'S5'), ('norm_mean', '<f8')]
            events = events.astype(numpy.dtype(new_dtype))
            events['transcriptomic_position'] += 2  # the middle position of 5-mers.
        with locks['hdf5'], h5py.File(out_paths['hdf5'],'a') as f:
            hf_tx = f.require_group('%s/%s' %(tx_id,read_id))
            if 'events' in hf_tx:
                continue
            else:
                hf_tx['events'] = events

    with locks['log'], open(out_paths['log'],'a') as f:
        f.write(tx_id + '\n')

    
def parallel_reformat_eventalign(data_dir,out_dir,n_processes):
    # Create output paths and locks.
    out_paths,locks = dict(),dict()
    for out_filetype in ['hdf5','log']:
        out_paths[out_filetype] = os.path.join(out_dir,'eventalign.%s' %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
                
    # Writing the starting of the files.
    h5py.File(out_paths['hdf5'],'w').close()
    open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue,task_function=reformat_eventalign,locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    # Load tasks into task_queue.
    with open('/ploy_ont_data/transcripts.csv','r') as f:
        reader = csv.reader(f)
        transcript_ids = sum(list(reader),[]) # flatten the list.
        
    # transcript_ids = os.listdir(data_dir)
    # print(len(transcript_ids),'transcripts')
    
    ###
    transcript_ids_processed = []
    for tx_id in transcript_ids:
        tx_dir = os.path.join(data_dir,tx_id)
        if not os.path.exists(tx_dir):
            continue
        transcript_ids_processed += [tx_id]
        task_queue.put((tx_id,tx_dir,out_paths)) # Blocked if necessary until a free slot is available. 


    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()
    
    with open(out_paths['log'],'a') as f:
        f.write('Total %d transcripts.\n' %len(transcript_ids_processed))
        f.write(helper.decor_message('successfully finished'))

    
    
def main():
    
    args = get_args()
    #
    n_processes = args.n_processes
    resume = args.resume
    
    data_dir = args.data_dir

    eventalign_filepath = args.eventalign
    summary_filepath = args.summary
    bamtx_filepath = args.bamtx
    out_dir = args.out_dir
    ensembl_version = args.ensembl
    gt_mapping_dir = args.mapping
    read_count_min = 30 #args.read_count_min

    misc.makedirs(out_dir) #todo: check every level.
    
    # (0)
    # parallel_reformat_eventalign(data_dir,out_dir,n_processes)
    
    # (1) For each read, combine multiple events aligned to the same positions, the results from nanopolish eventalign, into a single event per position.
    # parallel_combine(eventalign_filepath,summary_filepath,out_dir,n_processes,resume)
    
    # (2) Generate read count from the bamtx file.
    read_count_filepath = os.path.join(out_dir,'read_count.csv')
    if os.path.exists(read_count_filepath):
        df_count = pandas.read_csv(read_count_filepath)
    else:
        df_count = count_reads(ensembl_version,bamtx_filepath,out_dir)

    # (3) Create a .json file, where the info of all reads are stored per position, for modelling.
    parallel_preprocess(df_count,gt_mapping_dir,out_dir,n_processes,read_count_min)
    
if __name__ == '__main__':
    """
    Usage:
        xpore-backdoor --data_dir eventalign_by_transcript --eventalign eventalign.txt --summary summary.txt --mapping /ploy_ont_workspace/out/Release_v1_0/statCompare/data/mapping --bamtx bamtx/aligned.bam --out_dir dataprep --n_processes 8
    """
    main()

    
