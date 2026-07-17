"""Tests for diffmod.io.load_data -- assembling per-site read data for modelling.

load_data() collects, per (position, kmer) present in every run, the reads for each
condition, applies the min/max read-count rules, and returns model-ready arrays. These tests
cover the v2.2 change (over-cap runs are truncated, not dropped) alongside the min_count
filter, output shapes, the pooled round-robin selection (when users have multiple samples per condition), and the >=2-conditions requirement.
"""
import numpy as np  # noqa: F401  (imported for parity with the module under test)

from xpore.diffmod.io import load_data

IDX = 'geneX'
POS = 100
KMER = 'GGACT'


def _wrap(values):
    """Nest a read-value list into the {idx: {pos: {kmer: values}}} shape load_data reads."""
    return {IDX: {POS: {KMER: list(values)}}}


def test_over_max_count_truncates_not_drops():
    # A site in a run has MORE than max_count reads. The v2.1 behaviour dropped the site. 
    # the new behaviour truncates each site to max_count.
    data_dict = {
        ('KO', 'r1'): _wrap([0.1, 0.2, 0.3, 0.4, 0.5]),   # 5 > 3
        ('WT', 'r2'): _wrap([1.1, 1.2, 1.3, 1.4]),         # 4 > 3
    }
    data = load_data(IDX, data_dict, min_count=1, max_count=3)

    key = (IDX, POS, KMER)
    assert key in data                                     # site kept, not dropped
    result = data[key]
    assert len(result['y']) == 6                           # 3 (KO) + 3 (WT)
    assert (result['y_condition_names'] == 'KO').sum() == 3
    assert (result['y_condition_names'] == 'WT').sum() == 3


def test_run_below_min_count_dropped_site_survives():
    # KO has one run+site below min_count and one above. The below-min site+run is dropped.
    data_dict = {
        ('KO', 'r1'): _wrap([0.1]),               # 1 < min_count -> dropped
        ('KO', 'r2'): _wrap([0.2, 0.3, 0.4]),     # 3 >= min_count -> kept
        ('WT', 'r3'): _wrap([1.1, 1.2, 1.3]),
    }
    data = load_data(IDX, data_dict, min_count=2, max_count=100)

    result = data[(IDX, POS, KMER)]
    ko_values = result['y'][result['y_condition_names'] == 'KO']
    assert sorted(ko_values) == [0.2, 0.3, 0.4]            # r1's single read excluded
    assert 'r1' not in result['run_names']


def test_basic_two_condition_shapes():
    data_dict = {
        ('KO', 'r1'): _wrap([0.1, 0.2, 0.3]),
        ('WT', 'r2'): _wrap([1.1, 1.2, 1.3, 1.4]),
    }
    data = load_data(IDX, data_dict, min_count=1, max_count=100)

    result = data[(IDX, POS, KMER)]
    assert len(result['y']) == 7
    assert result['x'].shape == (7, 2)                     # one-hot over 2 conditions
    assert result['r'].shape == (7, 2)                     # one-hot over 2 runs
    assert result['condition_names'] == ['KO', 'WT']
    assert result['run_names'] == ['r1', 'r2']


def test_pooling_over_cap_round_robin():
    # KO's two runs total 7 reads > max_count 4. Round-robin selection takes reads evenly
    # across the runs (2 + 2), not the first 4 (which would be 3 from r1 + 1 from r2).
    data_dict = {
        ('KO', 'r1'): _wrap([1.0, 2.0, 3.0]),
        ('KO', 'r2'): _wrap([10.0, 20.0, 30.0, 40.0]),
        ('WT', 'r3'): _wrap([100.0, 200.0]),
    }
    data = load_data(IDX, data_dict, min_count=1, max_count=4, pooling=True)

    result = data[(IDX, POS, KMER)]
    ko = result['y_condition_names'] == 'KO'
    assert ko.sum() == 4                                   # capped at max_count
    ko_runs = result['y_run_names'][ko]
    assert (ko_runs == 'r1').sum() == 2                    # evenly split -> round-robin
    assert (ko_runs == 'r2').sum() == 2


def test_single_condition_site_skipped():
    # Only one condition present -> fewer than 2 conditions -> the site is skipped entirely.
    data_dict = {
        ('KO', 'r1'): _wrap([0.1, 0.2, 0.3]),
        ('KO', 'r2'): _wrap([0.4, 0.5, 0.6]),
    }
    data = load_data(IDX, data_dict, min_count=1, max_count=100)
    assert len(data) == 0
