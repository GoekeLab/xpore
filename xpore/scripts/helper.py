import gzip
import multiprocessing
import numpy
import os
import pandas
from functools import reduce
from collections import defaultdict


class EventalignFile:
    
    def __init__(self, fn):
        self._fn = fn
        self._open()

    def _open(self):
        fn = self._fn
        if os.path.splitext(fn)[1] == '.gz':
            self._handle = gzip.open(fn)
            self._decode_method = bytes.decode
        else:
            self._handle = open(fn)
            self._decode_method = str

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()

    def close(self):
        self._handle.close()

    def __iter__(self):
        return self

    def __next__(self):
        return self._decode_method(next(self._handle))


def decor_message(text,opt='simple'):
    text = text.upper()
    if opt == 'header':
        return text
    else:
        return '--- ' + text + ' ---\n'

def end_queue(task_queue,n_processes):
    for _ in range(n_processes):
        task_queue.put(None)
    return task_queue
        
def get_ids(f_index,data_info): #todo
    df_list = []
      
    for condition_name, run_names in data_info.items():
        list_of_set_ids = []
        for run_name in run_names:
            list_of_set_ids += [set(f_index[run_name].keys())]
        # ids = reduce(lambda x,y: x.intersection(y), list_of_set_ids)
        ids = reduce(lambda x,y: x.union(y), list_of_set_ids)
        df_list += [pandas.DataFrame({'ids':list(ids),condition_name:[1]*len(ids)})]
    df_merged = reduce(lambda  left,right: pandas.merge(left,right,on=['ids'], how='outer'), df_list).fillna(0).set_index('ids')
    return sorted(list(df_merged[df_merged.sum(axis=1) >= 2].index)) # At least two conditions.
    
## tmp: to remove
# def get_gene_ids(config_filepath): # arguments are not used.
#     import os,pandas
#     from ..diffmod.configurator import Configurator
#     # config_filepath = '/ploy_ont_workspace/github/experiments/Release_v1_0/config_manuscript/gmm_HEK293T-KO_HEK293T-WT_HEPG2-WT_K562-WT_A549-WT_MCF7-WT_reps_v01.ini'
#     config = Configurator(config_filepath) 
#     paths = config.get_paths()
#     info = config.get_info()
#     criteria = config.get_criteria() 
#     df_gt_ids = pandas.read_csv('/ploy_ont_workspace/out/Release_v1_0/statCompare/data/mapping_gt_ids.csv')    
#     gene_ids = set(df_gt_ids['g_id'].unique())
#     read_count_sum_min,read_count_sum_max = criteria['read_count_min'],criteria['read_count_max']
#     df_read_count = {}
#     for run_name in set(info['run_names']):
#         read_count_filepath = os.path.join(paths['data_dir'],run_name,'summary','read_count_per_gene.csv')
#         df_read_count[run_name] = pandas.read_csv(read_count_filepath).set_index('g_id') 

#     for run_name in set(info['run_names']):
#         df_read_count[run_name].reset_index(inplace=True)
#         cond = df_read_count[run_name]['n_reads'] >= criteria['read_count_min']
#         cond &= df_read_count[run_name]['n_reads'] <= criteria['read_count_max']
#         gene_ids = gene_ids.intersection(set(df_read_count[run_name].loc[cond,'g_id'].values))

#     gene_ids = sorted(list(gene_ids))
#     return gene_ids
##

class Consumer(multiprocessing.Process):
    """ For parallelisation """
    
    def __init__(self,task_queue,task_function,locks=None,result_queue=None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.locks = locks
        self.task_function = task_function
        self.result_queue = result_queue
        
    def run(self):
        proc_name = self.name
        while True:
            next_task_args = self.task_queue.get()
            if next_task_args is None:
                self.task_queue.task_done()
                break
            result = self.task_function(*next_task_args,self.locks)    
            self.task_queue.task_done()
            if self.result_queue is not None:
                self.result_queue.put(result)

def read_last_line(filepath): # https://stackoverflow.com/questions/3346430/what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file/3346788
    if not os.path.exists(filepath):
        return
    with open(filepath, "rb") as f:
        first = f.readline()        # Read the first line.
        if first == b'':
            return
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
        last = f.readline()         # Read last line.
    return last

def is_successful(filepath):
    return read_last_line(filepath) == b'--- SUCCESSFULLY FINISHED ---\n'
    