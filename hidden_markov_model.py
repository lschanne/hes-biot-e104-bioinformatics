"""
Find the most likely sequence of states for the following observed DNA sequence given the
three-state hidden Markov model (HMM) described below. Assume equal initial probabilities of
the 3 states. Provide 1) your Viterbi dynamic programming matrix including traceback and 2)
final state path (sequence of states). You can do it either manually, write your own program or
use Microsoft Excel. 3) Show a photo or a screenshot showing that you have done it yourself.
(30%)
Observed DNA sequence: TTTACGGT
"""

import numpy as np

states = {"T repeat": 0, "A linker": 1, "GC-rich": 2}
emissions = {"A": 0, "C": 1, "G": 2, "T": 3}
transition_probabilities = np.matrix([
    [0.8, 0.2, 0],
    [0, 0, 1.0],
    [0.1, 0, 0.9],
])
emission_probabilities = np.matrix([
    [0.2, 0, 0, 0.8],
    [1.0, 0, 0, 0],
    [0.1, 0.4, 0.4, 0.1],
])
n_states = len(states)
initial_probabilities = np.ones(shape=(n_states,)) / n_states
log2_init = np.log2(initial_probabilities)
log2_transition = np.log2(transition_probabilities)
log2_emit = np.log2(emission_probabilities)

observed_sequence = "TTTACGGT"
scores = np.zeros(shape=(n_states, len(observed_sequence)))
path = scores.copy()
first_emit = emissions[observed_sequence[0]]
scores[:,0] = log2_emit[:,first_emit].T + log2_init
for seq_idx, emission in enumerate(observed_sequence[1:], 1):
    emission_idx = emissions[emission]
    for state_idx in range(n_states):
        best_score = -np.inf
        best_path = np.nan
        for prev_state_idx in range(n_states):
            score = scores[prev_state_idx, seq_idx - 1] + log2_transition[prev_state_idx, state_idx]
            if score > best_score:
                best_score = score
                best_path = prev_state_idx
            path[state_idx, seq_idx] = best_path
            scores[state_idx, seq_idx] = best_score + log2_emit[state_idx, emission_idx]

print(scores.round(4))
print(path)
