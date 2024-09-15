"""
Calculating global and local alignment optimizations for DNA sequences.
Using Needleman-Wunsch global alignment.
Using Smith-Waterman local alignment.
"""

import numpy as np
from typing import (
    List,
    Optional,
    Tuple,
)

GAP = "-"
LEFT = 0
UP = 1
DIAGONAL = 2

def global_alignment(
        seq1: str,
        seq2: str,
        match_score: Optional[int]=3,
        mismatch_penalty: Optional[int]=-1,
        gap_penalty: Optional[int]=-7,
        gap_extension_penalty: Optional[int]=-1,
    ) -> Tuple[int, List[Tuple[str]]]:
    """
    Given two DNA sequences, compute the maximum global alignment score via the
    Needleman-Wunsch method. Return the score and all alignments that yield
    that score.

    Parameters
    ----------
    seq1: str
        Represents one of the two DNA sequences that need to be aligned
    seq2: str
        Represents the other DNA sequence that needs to be aligned
    match_score: int
        How much to reward the alignment for finding a match
    mismatch_penalty: int
        How much to penalize the alignment for mismatches
    gap_penalty: int
        How much to penalize the alignment for introducing a new gap
    gap_extension_penalty: int
        How much to penalize the alignment for extending an existing gap

    Returns
    -------
    [0]: int
        The maximum global alignment score achievable
    [1]: List[Tuple[str]]
        All possible alignments that yield the maximum global alignment score
        [0]: The alignment of the seq1 DNA sequence
        [1]: The alignment of the seq2 DNA sequence
    """

    # NOTE for now I am using only affine gap penalty for ease of
    # implementation; that means gap_extension_penalty is ignored and
    # only gap_penalty is used for both opening and extending gaps

    len1 = len(seq1)
    len2 = len(seq2)
    scores = np.zeros([len1 + 1, len2 + 1])

    scores[0, :] = np.arange(len2 + 1) * gap_penalty
    # TODO: introduce gap_extension_penalty
    # scores[0, :] = np.arange(len2 + 1) * gap_extension_penalty
    # scores[0, 1] = gap_penalty - gap_extension_penalty

    scores[:, 0] = np.arange(len1 + 1) * gap_penalty
    # TODO: introduce gap_extension_penalty
    # scores[:, 0] = np.arange(len1 + 1) * gap_extension_penalty
    # scores[1, 0] = gap_penalty - gap_extension_penalty

    paths = [[None for _ in range(len2)] for _ in range(len1)]

    for iii in range(len1):
        for jjj in range(len2):
            # TODO: introduce gap_extension_penalty
            left_score = scores[iii, jjj + 1] + gap_penalty
            up_score = scores[iii + 1, jjj] + gap_penalty
            diag_score = scores[iii, jjj]
            if seq1[iii] == seq2[jjj]:
                diag_score += match_score
            else:
                diag_score += mismatch_penalty
            best_score = -np.inf
            # TODO when you introduce gap_extension_penalty; somehow the best
            # paths for an alignment needs to be able to consider if the a
            # prior path is introducing a new gap or extending an old one
            for score, path in zip(
                (left_score, up_score, diag_score), (UP, LEFT, DIAGONAL),
            ):
                if score > best_score:
                    best_score = score
                    best_paths = [path]
                elif score == best_score:
                    best_paths.append(path)
            scores[iii + 1, jjj + 1] = best_score
            paths[iii][jjj] = best_paths

    alignment_score = int(best_score)
    # TODO need some recursive programming to get all possible alignments from
    # the paths
    alignment1 = seq1[-1]
    alignment2 = seq2[-1]
    while iii and jjj:
        path = paths[iii][jjj][0]
        if path == UP:
            jjj -= 1
            alignment1 = GAP + alignment1
            alignment2 = seq2[jjj] + alignment2
        elif path == LEFT:
            iii -= 1
            alignment1 = seq1[iii] + alignment1
            alignment2 = GAP + alignment2
        elif path == DIAGONAL:
            iii -= 1
            jjj -= 1
            alignment1 = seq1[iii] + alignment1
            alignment2 = seq2[jjj] + alignment2
    if iii:
        alignment1 = seq1[:iii] + alignment1
        alignment2 = GAP * iii + alignment2
    elif jjj:
        alignment1 = GAP * jjj + alignment1
        alignment2 = seq2[:jjj] + alignment2
    return alignment_score, [(alignment1, alignment2)]

seq1="AACTT"
seq2="ACAT"
match_score=3
mismatch_penalty=-1
gap_penalty=-7

if __name__ == "__main__":
    # TODO add some argparse to make this a proper script
    alignment_score, alignments = global_alignment(
        seq1="AACTT",
        seq2="ACAT",
        match_score=3,
        mismatch_penalty=-1,
        gap_penalty=-7,
    )
    print(alignment_score)
    alignment1, alignment2 = alignments[0]
    print(alignment1)
    print(alignment2)
