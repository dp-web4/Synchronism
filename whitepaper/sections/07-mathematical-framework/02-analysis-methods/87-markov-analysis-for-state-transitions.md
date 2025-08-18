# Markov Analysis for State Transitions

State transitions of intent between ticks are modeled using Markov
chains. The transition probability matrix is given by:

P(r, t) = \[p_11(r, t) \... p_1n(r, t)\]

\[ \... \... \... \]

\[p_n1(r, t) \... p_nn(r, t)\]

where p_ij(r, t) is the probability of transitioning from state i to
state j at position r and time t.