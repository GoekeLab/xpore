"""Tests for dataprep.preprocess_tx -- per-transcript reshaping into per-position signal.

preprocess_tx() concatenates a transcript's reads, groups them by position (or by
(position, kmer) in model_kmer mode), caps reads per site, runs a kmer-consistency check,
and appends the result to data.json. These tests exercise the v2.2 changes: the per-site
read cap, the break->continue fix, and the model_kmer (position, kmer) grouping.
"""
from conftest import make_events_array, read_json_output

from xpore.scripts.dataprep import preprocess_tx


def test_per_site_readcount_max_caps_reads(out_paths, locks):
    # Five reads at one site; a per-site cap of 3 keeps only 3 values for that site.
    data_dict = {
        i: make_events_array([(102, 'GGACT', float(i))], kmer_col='reference_kmer')
        for i in range(5)
    }
    preprocess_tx('tx1', data_dict, 'reference_kmer', 3, out_paths, locks)

    data = read_json_output(out_paths['json'])['tx1']
    assert len(data['102']['GGACT']) == 3


def test_readcount_max_none_keeps_all_reads(out_paths, locks):
    # readcount_max=None means no per-site cap -> every read is kept.
    data_dict = {
        i: make_events_array([(102, 'GGACT', float(i))], kmer_col='reference_kmer')
        for i in range(5)
    }
    preprocess_tx('tx1', data_dict, 'reference_kmer', None, out_paths, locks)

    data = read_json_output(out_paths['json'])['tx1']
    assert len(data['102']['GGACT']) == 5


def test_model_kmer_mode_keeps_both_orientations(out_paths, locks):
    # Forward and reverse reads at the same position carry different model_kmers. In
    # model_kmer mode reads are grouped by (position, kmer), so both orientations are kept.
    data_dict = {
        0: make_events_array([(102, 'GGACT', 100.0)], kmer_col='model_kmer'),
        1: make_events_array([(102, 'AGTCC', 88.0)], kmer_col='model_kmer'),
    }
    preprocess_tx('tx1', data_dict, 'model_kmer', None, out_paths, locks)

    data = read_json_output(out_paths['json'])['tx1']
    assert set(data['102'].keys()) == {'GGACT', 'AGTCC'}
