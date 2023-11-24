import pandas as pd
import numpy as np
import pytest
import os
import shutil
from xpore.scripts import diffmod

@pytest.fixture
def diffmod_args():
    class DiffmodArgs:
        outfile=open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'neues_config.yml'),'w')
        for ln in open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/config.yml"),'r'):
            if 'rep1' in ln:
                ln=ln.split(':')
                outfile.write(': '.join([ln[0], os.path.join(os.path.abspath(os.path.dirname(__file__)),ln[-1].strip())])+'\n')
            else:
                outfile.write(ln)
        outfile.close()
        config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'neues_config.yml')
        n_processes = 2
        save_models = False
        resume = False
        ids = []
    return DiffmodArgs

# the following runs the updated xpore diffmod from the PR and compares its diffmod.table to the diffmod.table generated
# by the original xpore diffmod (as long as the z-scores, DMRs are not too off, it should be fine)
def test_diffmod(diffmod_args):
    diffmod.diffmod(diffmod_args)

    #assert(os.path.exists(os.path.join(os.path.abspath(os.path.dirname(__file__)), "out/diffmod.table")))
    test_diffmod_path=os.path.join(os.getcwd(), "out/diffmod.table")
    original_diffmod_path=os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/original_diffmod.table")
    test_diffmod_table=pd.read_csv(test_diffmod_path).sort_values(["id", "position", "kmer"]).reset_index(drop=True)
    original_diffmod_table=pd.read_csv(original_diffmod_path).sort_values(["id", "position", "kmer"]).reset_index(drop=True)

    assert(np.all(original_diffmod_table["id"] == test_diffmod_table["id"]))
    assert(np.allclose(original_diffmod_table["diff_mod_rate_KO_vs_WT"], test_diffmod_table["diff_mod_rate_KO_vs_WT"]))
    assert(np.allclose(original_diffmod_table["z_score_KO_vs_WT"], test_diffmod_table["z_score_KO_vs_WT"]))
