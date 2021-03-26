import numpy as np
from itertools import product

NUM_NEIGHBORING_FEATURES = 1
CENTER_MOTIFS = [['A', 'G', 'T'], ['G', 'A'], ['A'], ['C'], ['A', 'C', 'T']]
FLANKING_MOTIFS = [['G', 'A', 'C', 'T'] for i in range(NUM_NEIGHBORING_FEATURES)]
ALL_KMERS = list(["".join(x) for x in product(*(FLANKING_MOTIFS + CENTER_MOTIFS + FLANKING_MOTIFS))])
ALL_KMERS = np.unique(np.array(list(map(lambda x: [x[i:i+5] for i in range(len(x) -4)], 
                                    ALL_KMERS))).flatten())
KMER_TO_INT = {ALL_KMERS[i]: i for i in range(len(ALL_KMERS))}
INT_TO_KMER =  {i: ALL_KMERS[i] for i in range(len(ALL_KMERS))}
M6A_KMERS = ["".join(x) for x in product(*CENTER_MOTIFS)]
