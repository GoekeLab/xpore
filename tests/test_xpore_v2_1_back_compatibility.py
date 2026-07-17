"""Back-compatibility tests: the transcriptome path still matches xPore v2.1.

Freezes outputs captured from xPore v2.1 and asserts the current (v2.2)
code reproduces them -- i.e. the v2.2 changes did not alter transcriptome (reference_kmer)
behaviour. 

Two levels are pinned:
  * combine()       -- the per-read event merge          (fixtures/_regenerate_combine_v2_1.py)
  * preprocess_tx() -- the per-position reshape to data.json
                       (combine + preprocess_tx; fixtures/_regenerate_data_json_v2_1.py)

The data.json comparison is order-insensitive within each site: preprocess_tx orders the
reads at a position via np.argsort, whose tie-breaking is not guaranteed stable across numpy
versions, and read order at a site is not meaningful to diffmod. What must be preserved is the
SET of values at each (position, kmer).
"""
import json
from pathlib import Path

from xpore.scripts.dataprep import combine, preprocess_tx

FIXTURES = Path(__file__).parent / "fixtures"


def test_combine_transcriptome_matches_v2_1():
    reference = json.loads((FIXTURES / "combine_transcriptome_v2_1.json").read_text())

    np_events, kmer_col = combine(reference["input_eventalign"], kmer_source="reference_kmer")

    assert kmer_col == reference["kmer_col"]
    actual = [
        [str(row["transcript_id"]), int(row["transcriptomic_position"]),
         str(row[kmer_col]), float(row["norm_mean"])]
        for row in np_events
    ]
    assert actual == reference["expected"]


def _sort_site_values(data_json):
    """Sort the read values at each (position, kmer) so comparison ignores within-site order."""
    return {
        tx_id: {pos: {kmer: sorted(vals) for kmer, vals in kmers.items()}
                for pos, kmers in positions.items()}
        for tx_id, positions in data_json.items()
    }


def test_preprocess_tx_data_json_matches_v2_1(out_paths, locks):
    reference = json.loads((FIXTURES / "data_json_transcriptome_v2_1.json").read_text())

    # Rebuild the per-read arrays the pipeline feeds preprocess_tx, using the current combine().
    data_dict, kmer_col = {}, "reference_kmer"
    for read in reference["reads"]:
        np_events, kmer_col = combine(read["events"], kmer_source="reference_kmer")
        if np_events.size > 1:                 # mirror the caller's inclusion rule
            data_dict[read["read_index"]] = np_events

    # readcount_max=None matches v2.1 (no per-site cap), so both reshape the same reads.
    preprocess_tx(reference["tx_id"], data_dict, kmer_col, None, out_paths, locks)

    actual = {}
    with open(out_paths["json"]) as f:
        for line in f:
            if line.strip():
                actual.update(json.loads(line))

    assert _sort_site_values(actual) == _sort_site_values(reference["expected"])
