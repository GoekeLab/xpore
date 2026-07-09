"""Shared helpers and fixtures for the xpore test suite.

The functions here build the small synthetic inputs the tests need:

* ``make_eventalign_line`` / ``make_eventalign_str`` -- eventalign text, the input to
  ``dataprep.combine``.
* ``make_events_array`` -- a single read's event array (what ``combine`` returns and what
  ``preprocess_tx`` consumes), built with fixed string widths so several reads can be
  ``np.concatenate``-d together.
* ``read_json_output`` -- parse the ``data.json`` that ``preprocess_tx`` writes.

The ``locks`` and ``out_paths`` fixtures provide the two arguments ``preprocess_tx`` writes
through, so the tests never touch multiprocessing.
"""
import json
import threading

import numpy as np
import pytest


# Column order of a nanopolish/f5c eventalign.txt, exactly as combine() parses it.
EVENTALIGN_COLUMNS = [
    'contig', 'position', 'reference_kmer', 'read_index', 'strand',
    'event_index', 'event_level_mean', 'event_stdv', 'event_length',
    'model_kmer', 'model_mean', 'model_stdv', 'standardized_level',
    'start_idx', 'end_idx',
]

# Sensible defaults for one eventalign row; override any field per call. With start_idx=0,
# end_idx=5 the per-row length is 5, so a lone row's norm_mean == its event_level_mean.
_EVENTALIGN_DEFAULTS = {
    'contig': 'ENST1',
    'position': 100,
    'reference_kmer': 'GGACT',
    'read_index': 0,
    'strand': 't',
    'event_index': 1,
    'event_level_mean': 120.0,
    'event_stdv': 2.0,
    'event_length': 0.005,
    'model_kmer': 'GGACT',
    'model_mean': 120.0,
    'model_stdv': 1.5,
    'standardized_level': 0.5,
    'start_idx': 0,
    'end_idx': 5,
}



#helper functions for the combine() tests (in test_combine.py)
def make_eventalign_line(**overrides):
    """Return one tab-separated eventalign row (no trailing newline)."""
    row = dict(_EVENTALIGN_DEFAULTS, **overrides)
    return '\t'.join(str(row[col]) for col in EVENTALIGN_COLUMNS)

def make_eventalign_str(rows):
    """Join a list of row-override dicts into one eventalign string (what combine expects)."""
    return '\n'.join(make_eventalign_line(**row) for row in rows)



#helper functions for the preprocess_tx() tests (in test_preprocess_tx.py)
def _events_dtype(kmer_col):
    # Matches combine()'s output fields. Fixed string widths let arrays for different reads
    # be np.concatenate-d inside preprocess_tx.
    return np.dtype([
        ('transcript_id', '<U15'),
        ('transcriptomic_position', '<i8'),
        (kmer_col, '<U5'),
        ('norm_mean', '<f8'),
    ])

def make_events_array(rows, kmer_col='reference_kmer', transcript_id='tx1'):
    """Build one read's event array.

    ``rows`` is a list of ``(transcriptomic_position, kmer, norm_mean)`` tuples.
    """
    records = [(transcript_id, pos, kmer, norm_mean) for pos, kmer, norm_mean in rows]
    return np.array(records, dtype=_events_dtype(kmer_col))

def read_json_output(json_path):
    """Parse preprocess_tx's data.json into ``{tx_id: {pos: {kmer: [values]}}}``.

    NOTE: JSON object keys are always strings, so genomic positions come back as ``str``
    (e.g. ``'102'``), not ``int``.
    """
    result = {}
    with open(json_path) as f:
        for line in f:
            line = line.strip()
            if line:
                result.update(json.loads(line))
    return result

#Fixtures for the preprocess_tx() tests (test_preprocess_tx.py, test_xpore_v2_1_back_compatibility.py).
@pytest.fixture
def locks():
    """The four locks preprocess_tx uses as context managers."""
    return {name: threading.Lock() for name in ('json', 'index', 'readcount', 'log')}


@pytest.fixture
def out_paths(tmp_path):
    """The four output files preprocess_tx appends to (created lazily in append mode)."""
    return {name: str(tmp_path / f'data.{name}')
            for name in ('json', 'index', 'readcount', 'log')}
