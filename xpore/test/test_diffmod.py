import pandas as pd
import numpy as np
import pytest
import os
import shutil
import subprocess
from xpore.scripts import diffmod

#@pytest.fixture
#def dataprep_args():
#    return {
#            'eventalign': os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/eventalign.txt"),
#            'out_dir': 'dataprep',
#            'readcount_min': 1,
#            'readcount_max': 1000,
#            'n_processes': 2
#            }

@pytest.fixture
def diffmod_args():
    class DiffmodArgs:
        #config = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/config.yml")
        outfile=open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'neues_config.yml'),'w')
        for ln in open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/config.yml"),'r'):
            if 'rep1' in ln:
                ln=ln.split(':')
                outfile.write(': '.join([ln[0], os.path.join(os.path.abspath(os.path.dirname(__file__)),ln[-1].strip())]))
            else:
                outfile.write(ln)
        outfile.close()
        config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'neues_config.yml')
        n_processes = 2
        save_models = False
        resume = False
        ids = []
    return DiffmodArgs

def test_diffmod(diffmod_args):
    #cmd = ['xpore diffmod --config',os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/config.yml"),'--n_processes 4']
    #subprocess.run(' '.join(cmd), shell=True, check=True)
    diffmod.diffmod(diffmod_args)
