"""Regenerate the xPore v2.1 data.json reference fixture (transcriptome path).

`data_json_transcriptome_v2_1.json` freezes the data.json that combine() + preprocess_tx()
produce on a fixed multi-read transcriptome input, captured from xPore v2.1 (the pre-PR
upstream version). The test
(../test_xpore_v2_1_back_compatibility.py::test_preprocess_tx_data_json_matches_v2_1)
asserts the current (v2.2) code reproduces it.

The committed fixture was captured from upstream/master (GoekeLab/xpore, xPore v2.1) via a
git worktree:

    git worktree add /tmp/xpore-v2.1 upstream/master
    PYTHONPATH=/tmp/xpore-v2.1 python tests/fixtures/_regenerate_data_json_v2_1.py
    git worktree remove /tmp/xpore-v2.1

Handles both the v2.1 and v2.2 combine()/preprocess_tx() signatures, so it can be pointed at
either version; it prints the xpore it imported so you can confirm the source. (Underscore
prefix keeps pytest from collecting this as a test module.)
"""
import inspect
import json
import os
import tempfile
import threading

import xpore
from xpore.scripts.dataprep import combine, preprocess_tx

TX_ID = "ENST1"


def _line(pos, kmer, read_index, event_index, mean, start, end):
    # One eventalign row (reference_kmer == model_kmer, i.e. forward transcriptome alignment).
    return "\t".join(str(x) for x in [
        "ENST1", pos, kmer, read_index, "t", event_index, mean, 2.0, 0.01,
        kmer, mean, 1.5, 0.5, start, end,
    ])


# Three reads with overlapping positions, so multiple reads land on the same site and
# preprocess_tx's per-position grouping is exercised. Each read spans >1 position (size > 1)
# so the caller's inclusion rule keeps it.
READS = [
    {"read_index": 0, "events": "\n".join([
        _line(100, "GGACT", 0, 1, 120.0, 1000, 1005),
        _line(101, "GACTA", 0, 2, 95.0, 1005, 1010),
        _line(102, "ACTAG", 0, 3, 104.0, 1010, 1015),
    ])},
    {"read_index": 1, "events": "\n".join([
        _line(100, "GGACT", 1, 1, 121.0, 2000, 2005),
        _line(101, "GACTA", 1, 2, 96.0, 2005, 2010),
        _line(102, "ACTAG", 1, 3, 105.0, 2010, 2015),
    ])},
    {"read_index": 2, "events": "\n".join([
        _line(101, "GACTA", 2, 1, 94.0, 3000, 3005),
        _line(102, "ACTAG", 2, 2, 106.0, 3005, 3010),
        _line(103, "CTAGT", 2, 3, 110.0, 3010, 3015),
    ])},
]


def _run_combine(events):
    result = combine(events)               # v2.1 returns array; v2.2 returns (array, kmer_col)
    return result[0] if isinstance(result, tuple) else result


def _call_preprocess_tx(tx_id, data_dict, kmer_col, out_paths, locks):
    # v2.2 added kmer_col + readcount_max params; v2.1 has neither. readcount_max=None on v2.2
    # matches v2.1 (no per-site cap), so both reshape the same reads.
    if "kmer_col" in inspect.signature(preprocess_tx).parameters:
        preprocess_tx(tx_id, data_dict, kmer_col, None, out_paths, locks)
    else:
        preprocess_tx(tx_id, data_dict, out_paths, locks)


def main():
    data_dict, kmer_col = {}, "reference_kmer"
    for read in READS:
        np_events = _run_combine(read["events"])
        kmer_col = np_events.dtype.names[2]
        if np_events.size > 1:             # mirror the caller's inclusion rule
            data_dict[read["read_index"]] = np_events

    with tempfile.TemporaryDirectory() as tmp:
        out_paths = {name: os.path.join(tmp, "data.%s" % name)
                     for name in ("json", "index", "readcount", "log")}
        locks = {name: threading.Lock() for name in out_paths}
        _call_preprocess_tx(TX_ID, data_dict, kmer_col, out_paths, locks)
        expected = {}
        with open(out_paths["json"]) as f:
            for line in f:
                if line.strip():
                    expected.update(json.loads(line))

    fixture = {
        "_comment": "xPore v2.1 data.json reference captured from upstream/master; regenerate via _regenerate_data_json_v2_1.py.",
        "tx_id": TX_ID,
        "reads": READS,
        "expected": expected,
    }
    out_path = os.path.join(os.path.dirname(__file__), "data_json_transcriptome_v2_1.json")
    with open(out_path, "w") as f:
        json.dump(fixture, f, indent=2)
        f.write("\n")

    print("imported xpore from:", xpore.__file__)
    print("wrote:", out_path)
    print(json.dumps(expected, indent=2))


if __name__ == "__main__":
    main()
