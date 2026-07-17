"""Tests for dataprep.combine -- per-read event merging and kmer-column selection.

combine() filters events where the reference and model kmers agree, merges multiple events
at the same position into one length-weighted norm_mean, and (new in v2.2) chooses which kmer
column to keep based on --kmer_source, accepting reverse-complement matches for genome reads.
"""
from conftest import make_eventalign_line, make_eventalign_str

from xpore.scripts.dataprep import combine

def test_transcriptome_basic_match():
    # reference_kmer == model_kmer: the transcriptome case, works with the default source.
    events_str = make_eventalign_line(
        contig='ENST1', position=100, read_index=0,
        reference_kmer='GGACT', model_kmer='GGACT',
        event_level_mean=120.5, start_idx=0, end_idx=5,
    )
    result = combine(events_str, kmer_source='reference_kmer')

    assert result is not None
    np_events, kmer_col = result
    assert kmer_col == 'reference_kmer'
    assert np_events.size == 1
    assert np_events['transcript_id'][0] == 'ENST1'
    # combine reports the middle base of the 5-mer: position + 2.
    assert np_events['transcriptomic_position'][0] == 102
    assert np_events['reference_kmer'][0] == 'GGACT'
    assert np_events['norm_mean'][0] == 120.5

def test_reverse_complement_dropped_with_reference_kmer():
    # The same reverse read under reference_kmer source: no exact match, so it is filtered out
    # and combine returns None. This is the behaviour in v2.1. 
    events_str = make_eventalign_line(reference_kmer='GGACT', model_kmer='AGTCC')
    assert combine(events_str, kmer_source='reference_kmer') is None

def test_reverse_complement_kept_with_model_kmer():
    # A reverse-oriented read: reference_kmer is the reverse complement of model_kmer (from the read)
    # (revcomp('GGACT') == 'AGTCC'). With --kmer_source model_kmer this read must be KEPT.
    events_str = make_eventalign_line(
        contig='chr1', position=100, read_index=5,
        reference_kmer='GGACT', model_kmer='AGTCC',
        event_level_mean=88.0, start_idx=2000, end_idx=2005,
    )
    result = combine(events_str, kmer_source='model_kmer')

    assert result is not None
    np_events, kmer_col = result
    assert kmer_col == 'model_kmer'
    assert np_events.size == 1
    # the kmer recorded is the model_kmer (the kmer from the read).
    assert np_events['model_kmer'][0] == 'AGTCC'
    assert np_events['norm_mean'][0] == 88.0


def test_direct_match_with_model_kmer():
    # Test forward reads (reference == model) are kept when kmer_source='model_kmer'.
    events_str = make_eventalign_line(reference_kmer='GGACT', model_kmer='GGACT')
    result = combine(events_str, kmer_source='model_kmer')

    assert result is not None
    np_events, kmer_col = result
    assert kmer_col == 'model_kmer'
    assert np_events['model_kmer'][0] == 'GGACT'


def test_no_match_returns_none():
    # Neither an exact nor a reverse-complement match -> nothing survives the filter.
    events_str = make_eventalign_line(reference_kmer='GGACT', model_kmer='TTTTT')
    assert combine(events_str, kmer_source='model_kmer') is None


def test_norm_mean_is_length_weighted_ref_kmer():
    # when kmer_source='reference_kmer'
    # Two events at the same position within a **forward** read  (read, contig, position, kmer) group with different (start,end)
    # widths -> norm_mean is the length-weighted mean, not the plain mean:
    #   (100*1 + 130*3) / (1 + 3) == 122.5     (the plain mean would be 115.0)
    events_str = make_eventalign_str([
        dict(read_index=0, position=100, reference_kmer='GGACT', model_kmer='GGACT',
             event_level_mean=100.0, start_idx=0, end_idx=1),
        dict(read_index=0, position=100, reference_kmer='GGACT', model_kmer='GGACT',
             event_level_mean=130.0, start_idx=1, end_idx=4),
    ])
    np_events, _ = combine(events_str, kmer_source='reference_kmer')

    assert np_events.size == 1
    assert np_events['norm_mean'][0] == 122.5

def test_norm_mean_is_length_weighted_model_kmer():
    # when kmer_source='model_kmer'
    # Two events at the same position within a **reverse** read (same read, contig, position, kmer) group with different (start,end)
    # widths -> norm_mean is the length-weighted mean, not the plain mean:
    #   (100*1 + 130*3) / (1 + 3) == 122.5     (the plain mean would be 115.0)
    events_str = make_eventalign_str([
        dict(contig='chr1', read_index=0, position=100, reference_kmer='GGACT', model_kmer='AGTCC',
             event_level_mean=100.0, start_idx=0, end_idx=1),
        dict(contig='chr1', read_index=0, position=100, reference_kmer='GGACT', model_kmer='AGTCC',
             event_level_mean=130.0, start_idx=1, end_idx=4),
    ])
    np_events, _ = combine(events_str, kmer_source='model_kmer')

    assert np_events.size == 1
    assert np_events['norm_mean'][0] == 122.5
