"""Regenerate the xPore v2.1 reference fixture for combine() on the transcriptome path.

`combine_transcriptome_v2_1.json` freezes combine()'s output on a fixed transcriptome
eventalign input, captured from xPore v2.1 (the pre-PR upstream version). The test
(../test_xpore_v2_1_back_compatibility.py) asserts the current (v2.2) code reproduces it
byte-for-byte -- a regression guard that the v2.2 changes did not alter transcriptome
(reference_kmer) behaviour.

The committed fixture was captured from upstream/master (GoekeLab/xpore, xPore v2.1) via a
git worktree:

    git worktree add /tmp/xpore-v2.1 upstream/master
    PYTHONPATH=/tmp/xpore-v2.1 python tests/fixtures/_regenerate_combine_v2_1.py
    git worktree remove /tmp/xpore-v2.1

Point PYTHONPATH at whichever checkout you want as the reference; the script prints the actual
xpore it imported so you can confirm the source. (Underscore prefix keeps pytest from
collecting this as a test module.)
"""
import json
import os

import xpore
from xpore.scripts.dataprep import combine

# A small but non-trivial transcriptome eventalign snippet for ONE read (read_index 0):
# reference_kmer == model_kmer throughout (forward transcriptome alignment), with position
# 100 carrying two events so the length-weighted merge is exercised.
INPUT_EVENTALIGN = "\n".join([
    "ENST1\t100\tGGACT\t0\tt\t1\t120.5\t2.0\t0.01\tGGACT\t120.0\t1.5\t0.5\t1000\t1005",
    "ENST1\t100\tGGACT\t0\tt\t2\t118.0\t2.1\t0.01\tGGACT\t120.0\t1.5\t0.4\t1005\t1008",
    "ENST1\t101\tGACTA\t0\tt\t3\t95.0\t1.8\t0.01\tGACTA\t95.5\t1.4\t0.3\t1008\t1013",
    "ENST1\t102\tACTAG\t0\tt\t4\t104.5\t1.9\t0.01\tACTAG\t104.0\t1.6\t0.2\t1013\t1018",
])


def main():
    result = combine(INPUT_EVENTALIGN)
    # xPore v2.1 combine returns just the array; v2.2 returns (array, kmer_col).
    np_events = result[0] if isinstance(result, tuple) else result
    kmer_field = np_events.dtype.names[2]  # 'reference_kmer' on the transcriptome path

    expected = [
        [str(row["transcript_id"]), int(row["transcriptomic_position"]),
         str(row[kmer_field]), float(row["norm_mean"])]
        for row in np_events
    ]

    fixture = {
        "_comment": "xPore v2.1 reference output captured from upstream/master; regenerate via _regenerate_combine_v2_1.py.",
        "input_eventalign": INPUT_EVENTALIGN,
        "kmer_col": kmer_field,
        "expected": expected,
    }

    out_path = os.path.join(os.path.dirname(__file__), "combine_transcriptome_v2_1.json")
    with open(out_path, "w") as f:
        json.dump(fixture, f, indent=2)
        f.write("\n")

    print("imported xpore from:", xpore.__file__)
    print("wrote:", out_path)
    print(json.dumps(expected, indent=2))


if __name__ == "__main__":
    main()
