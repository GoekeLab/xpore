import requests
import os
import pytest
import subprocess

## stage nanopolish eventalign.txt files from demo test data
data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
r_KO=requests.get('https://github.com/yuukiiwa/test-datasets/raw/nanoseq/modification_nanopolish/HEK293T-METTL3-KO-rep1_eventalign.txt', allow_redirects=True)
open(os.path.join(data_dir,'HEK293T-METTL3-KO-rep1_eventalign.txt'), 'wb').write(r_KO.content)
r_WT=requests.get('https://github.com/yuukiiwa/test-datasets/raw/nanoseq/modification_nanopolish/HEK293T-WT-rep1_eventalign.txt', allow_redirects=True)
open(os.path.join(data_dir,'HEK293T-WT-rep1_eventalign.txt'), 'wb').write(r_WT.content)

ref_dir = os.path.join(os.path.dirname(__file__), 'test_data','reference')
if not os.path.exists(ref_dir):
    os.makedirs(ref_dir)
r_FASTA=requests.get('https://github.com/yuukiiwa/test-datasets/raw/nanoseq/reference/modification_transcriptome_subset.fa', allow_redirects=True)
open(os.path.join(ref_dir,'modification_transcriptome_subset.fa'), 'wb').write(r_FASTA.content)
r_GTF=requests.get('https://raw.githubusercontent.com/yuukiiwa/test-datasets/nanoseq/reference/modification_transcriptome_subset.gtf', allow_redirects=True)
open(os.path.join(ref_dir,'modification_transcriptome_subset.gtf'), 'wb').write(r_GTF.content)

def _create_diffmod_config(outdir,type):
    yml_path=os.path.join(outdir,type+'_config.yml')
    outfile=open(yml_path,"w")
    outfile.write(
    'notes: Pairwise comparison without replicates with default parameter setting.'+ '\n'+
    '\n'+
    'data:'+'\n'+
    '    KO:'+'\n'+
    '        rep1: '+os.path.join(outdir,'dataprepKO_'+type)+'\n'+
    '    WT:'+'\n'+
    '        rep1: '+os.path.join(outdir,'dataprepWT_'+type)+'\n'+
    '\n'+
    'out: '+os.path.join(outdir,'diffmod_'+type)+'\n'+
    '\n'+
    'method:'+'\n'+
    '    prefiltering:'+'\n'+
    '        method: t-test'+'\n'+
    '        threshold: 0.1'+'\n'+'\n'
    )
    return yml_path

def _readDiffmodTable(path):
    file=open(path,'r')
    file.readline()
    rank,asi=[],[]
    for ln in file:
        ln=ln.strip().split(",")
        rank.append(ln[0].split(".")[0]+"_"+ln[1]+"_"+ln[2])
        asi.append((ln[0].split(".")[0]+"_"+ln[1]+"_"+ln[2],float(ln[4])))
    return (rank,asi)

def _checkDiffmodPval(test_run_out_path,expected_out_path):
    test_run_rank,test_run_asi=_readDiffmodTable(test_run_out_path)
    expected_rank,expected_asi=_readDiffmodTable(expected_out_path)
    test_run_asi_sort=sorted(test_run_asi, key=lambda tup: tup[1])
    expected_asi_sort=sorted(expected_asi, key=lambda tup: tup[1])
    check_n_entries = len(test_run_asi_sort) == len(expected_asi_sort)
    n_same_entries = 0
    for i in range(len(test_run_asi_sort)):
        gene,pval=test_run_asi_sort[i]
        for j in range(len(expected_asi_sort)):
            exp_gene,exp_pval=expected_asi_sort[j]
            if gene == exp_gene and pval == exp_pval:
                n_same_entries += 1
    check_pval =  n_same_entries == len(test_run_asi_sort)
    if check_n_entries and check_pval:
        return True
    else:
        return False

@pytest.mark.parametrize("type",['transcript','gene'])
def test_xpore(type,tmpdir):
    eventalign_KO=os.path.join(os.path.dirname(__file__),'test_data','HEK293T-METTL3-KO-rep1_eventalign.txt')
    eventalign_WT=os.path.join(os.path.dirname(__file__),'test_data','HEK293T-WT-rep1_eventalign.txt')
    fasta=os.path.join(os.path.dirname(__file__),'test_data','reference','modification_transcriptome_subset.fa')
    gtf=os.path.join(os.path.dirname(__file__),'test_data','reference','modification_transcriptome_subset.gtf')

    outdir=str(tmpdir)
    ## run xpore-dataprep for each sample
    if type == 'transcript':
        dataprepKO_args = ['xpore-dataprep',
                           '--eventalign', eventalign_KO,
                           '--out_dir', os.path.join(outdir,'dataprepKO_'+type)]
        subprocess.run(' '.join(dataprepKO_args), shell=True, check=True)

        dataprepWT_args = ['xpore-dataprep',
                           '--eventalign', eventalign_WT,
                           '--out_dir', os.path.join(outdir,'dataprepWT_'+type)]
        subprocess.run(' '.join(dataprepWT_args), shell=True, check=True)
    elif type == 'gene':
        dataprepKO_args = ['xpore-dataprep',
                           '--eventalign', eventalign_KO,
                           '--out_dir', os.path.join(outdir,'dataprepKO_'+type),
                           '--genome',
                           '--gtf_path_or_url', gtf,
                           '--transcript_fasta_paths_or_urls', fasta]
        subprocess.run(' '.join(dataprepKO_args), shell=True, check=True)

        dataprepWT_args = ['xpore-dataprep',
                           '--eventalign', eventalign_WT,
                           '--out_dir', os.path.join(outdir,'dataprepWT_'+type),
                           '--genome',
                           '--gtf_path_or_url', gtf,
                           '--transcript_fasta_paths_or_urls', fasta]
        subprocess.run(' '.join(dataprepWT_args), shell=True, check=True)

    ## run diffmod
    config_yml_path = _create_diffmod_config(outdir,type)
    diffmod_args = ['xpore-diffmod',
                    '--config', config_yml_path]
    subprocess.run(' '.join(diffmod_args), shell=True, check=True)

    ## check diffmod output
    test_run_out_path = os.path.join(outdir,'diffmod_'+type,'diffmod.table')
    expected_out_path = os.path.join(os.path.dirname(__file__),'expected_outputs',type+'_diffmod_1_0.table')
    diffmod_check = _checkDiffmodPval(test_run_out_path,expected_out_path)
    assert(diffmod_check == True)
