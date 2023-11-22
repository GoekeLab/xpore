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
    return {
            'config': os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/config.yml"),
            'out_dir': 'diffmod',
            'n_processes': 2,
            'save_models': False,
            'resume': False,
            'ids': []
            }


def test_diffmod(diffmod_args):
    cmd = ['xpore diffmod --config',diffmod_args.config,'--out_dir diffmod']
    subprocess.run(' '.join(cmd), shell=True, check=True)
    #diffmod.diffmod(diffmod_args)
