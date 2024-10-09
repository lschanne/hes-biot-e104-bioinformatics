"""
Train a first-order Markov model from the following DNA sequence.
1) Provide a transition probability matrix rounded to 2 decimal places.
2) calculate the log2 probability of sequence GCACACA given your transition
probability matrix. Assume that the initial probabilities are equal for all
four states. Round to 2 decimal places. (20%)

DNA sequence:
TTTAAAATCCGGCGGCGATCCAATCCCCGTTAGTATAATTAATTTTTCAATCCCCCATTTAATAATACGGCCCGAAATTATTAACGGGCCCCCGGCCGCCCC
"""

import numpy as np

sequence = "TTTAAAATCCGGCGGCGATCCAATCCCCGTTAGTATAATTAATTTTTCAATCCCCCATTTAATAATACGGCCCGAAATTATTAACGGGCCCCCGGCCGCCCC"

states = "AGCT"
state_map = {state: idx for idx, state in enumerate(states)}
print(state_map)
n_states = len(states)
transition_counts = np.zeros(shape=(n_states, n_states))

prev_idx = state_map[sequence[0]]
for curr_state in sequence[1:]:
    curr_idx = state_map[curr_state]
    transition_counts[prev_idx][curr_idx] += 1
    prev_idx = curr_idx

print(transition_counts)
transition_probabilities = transition_counts / transition_counts.sum(axis=1)
print(transition_probabilities)
log_probabilities = np.log2(transition_probabilities)

new_sequence = "GCACACA"
log_prob = np.log2(1 / n_states)
prev_idx = state_map[new_sequence[0]]
for curr_state in new_sequence[1:]:
    curr_idx = state_map[curr_state]
    log_prob += log_probabilities[prev_idx][curr_idx]
    prev_idx = curr_idx
print(log_prob)

