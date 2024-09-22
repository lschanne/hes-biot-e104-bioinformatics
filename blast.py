"""
Perform your own manual ungapped BLAST search for the DNA sequence below
against the following "database" of 4 sequences. In other words, do not use an
existing BLAST program, but perform the three steps of the BLAST method
manually (or write your own BLAST program).

Use word length w=10 and word score threshold T=45, with match = 5,
mismatch = -2. Find a maximal segment pair that produce an alignment score of
at least S=55, with an extension termination threshold of 40. (20%)

1) First, show the list of words you used to search for in the database.
Then, provide 2) the name, 3) start and end positions and 4) (raw) alignment
score of your hit (e.g. “SEQ2 position 4 to 10, score=60”). Assume that the
position of the first base of each sequence in the database is 1.

Query sequence: CAAATCTTCAGGAC

Database:
SEQ1: CGCTATCGGCATCGACGATGCATCGTAGCAGTCTACTAGCATCGA
SEQ2: TACGGCGCACGCATCGAGTCATGCATCGACAAATCACACTCAGTCTG
SEQ3: CCAGCGAATCTTCAGGAAAAAATCTTGCCCATTATACCGCGCGATAT
SEQ4: CCATATTATCTTCTGCGCGCTATATAACTCGCAGTCGATGAGACTAGC

Hint: be careful calculating the score for 1 mismatch (which means 9 matches!)
"""

import re

WORD_LENGTH = 10
WORD_SCORE_THRESHOLD = 45
MATCH = 5
MISMATCH = -2
MIN_ALIGNMENT_SCORE = 55
TERMINATION_THRESHOLD = 40

QUERY_SEQUENCE = "CAAATCTTCAGGAC"
DATABASE = [
    "CGCTATCGGCATCGACGATGCATCGTAGCAGTCTACTAGCATCGA",
    "TACGGCGCACGCATCGAGTCATGCATCGACAAATCACACTCAGTCTG",
    "CCAGCGAATCTTCAGGAAAAAATCTTGCCCATTATACCGCGCGATAT",
    "CCATATTATCTTCTGCGCGCTATATAACTCGCAGTCGATGAGACTAGC",
]

# Step 1: compile a list of words of length WORD_LENGTH that would produce
# a score of at least T with the query sequence
max_mismatches = 0
score = MATCH * WORD_LENGTH
while score > WORD_SCORE_THRESHOLD:
    score = score - MATCH + MISMATCH
    if score > WORD_SCORE_THRESHOLD:
        max_mismatches += 1
print(f"{max_mismatches=}") # max_mismatches = 0

# since there can be no mismatches, the word list will just be a sliding
# window of the QUERY_SEQUENCE
word_list = []
print("\nword list:")
for idx in range(len(QUERY_SEQUENCE) - WORD_LENGTH + 1):
    word = QUERY_SEQUENCE[idx:idx + WORD_LENGTH]
    print(f"    {word}")
    word_list.append((idx, word))

# Step 2: Scan the sequences in the database to determine if they contain any
# word in the list
print("\nmatches:")
for query_start_idx, word in word_list:
    query_end_idx = query_start_idx + WORD_LENGTH
    for seq_num, seq in enumerate(DATABASE, 1):
        for match in re.finditer(word, seq):
            start_idx = match.start()
            print(f"    Trying {seq_num=} {query_start_idx=} {start_idx=}")

            # Step 3: If a sequence in the database contains a word in the
            # list, extend the alignment on both sides to determine if it can
            # produce a score of at least MIN_ALIGNMNET_SCORE
            end_idx = start_idx + WORD_LENGTH
            best_start_idx = start_idx
            score = best_score = WORD_LENGTH * MATCH
            for iii in range(min(start_idx, query_start_idx)):
                if (
                    seq[start_idx - iii - 1] ==
                    QUERY_SEQUENCE[query_start_idx - iii - 1]
                ):
                    score += MATCH
                else:
                    score += MISMATCH
                if score < TERMINATION_THRESHOLD:
                    break
                if score > best_score:
                    best_score = score
                    best_start_idx = start_idx - iii - 1
            score = best_score
            best_end_idx = end_idx
            max_range = min(
                len(seq) - end_idx,
                len(QUERY_SEQUENCE) - query_end_idx,
            )
            if max_range:
                for iii in range(max_range - 1):
                    if (
                        seq[end_idx + iii + 1] ==
                        QUERY_SEQUENCE[query_end_idx + iii + 1]
                    ):
                        score += MATCH
                    else:
                        score += MISMATCH
                    if score < TERMINATION_THRESHOLD:
                        break
                    if score > best_score:
                        best_score = score
                        best_end_idx = end_idx + iii + 1

            if best_score >= MIN_ALIGNMENT_SCORE:
                print(f"    SEQ{seq_num} position {best_start_idx + 1} to {best_end_idx + 1}, score={best_score}")
