"""
Given the following sequences:
    AAAAGG
    AATAGG
    ATTAAG
    AAAAAG
    CAACAG
    TATAAG
Create a position-specific scoring matrix (PSSM) with pseudocount 1.

Use log base 2 and round the final matrix to 1 decimal place.

Use background distribution of A: 0.23, C: 0.27, G: 0.27, T: 0.23 instead of an
equiprobable distribution.

Calculate the score of sequence AAAAGG based on your PSSM.
"""

import numpy as np

sequences = [
    "AAAAGG",
    "AATAGG",
    "ATTAAG",
    "AAAAAG",
    "CAACAG",
    "TATAAG",
]
seq_to_score = "AAAAGG"
seq_len = len(seq_to_score)
bases = "ACGT"
n_bases = len(bases)
base_map = {base: idx for idx, base in enumerate(bases)}
counts = np.ones(shape=(n_bases, seq_len))
for seq_idx in range(seq_len):
    for seq in sequences:
        base_idx = base_map[seq[seq_idx]]
        counts[base_idx, seq_idx] += 1
print(counts)
n_sequences = len(sequences)
probabilites = counts / (n_sequences + n_bases + np.zeros(shape=(1, seq_len)))
print(probabilites)
background_probs = np.matrix([[0.23, 0.27, 0.27, 0.23]]).T
print(background_probs)
pssm = np.log2(probabilites / background_probs)
print(pssm.round(1))

score = 0
for seq_idx, base in enumerate(seq_to_score):
    base_idx = base_map[base]
    score += pssm[base_idx, seq_idx]
print(score)
