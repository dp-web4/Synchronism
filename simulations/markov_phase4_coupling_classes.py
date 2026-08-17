#!/usr/bin/env python3
"""Phase 4 toy model for the Markov coherence/emergence arc.

Two persistent binary subsystems A and B are updated synchronously. Each retains memory
of its own current state and is also coupled to the other's current state by parameter J.

J > 0 favors alignment, J < 0 favors anti-alignment, J = 0 factorizes exactly.

Measures:
- stationary correlation and mutual information,
- relation-state residence times,
- observational transfer information I(A_{t+1};B_t|A_t),
- intervention-like causal sensitivity of A's next-state law to B.

Formalism toy only; not a physics model.
"""

import math
import numpy as np

STATES = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
INDEX = {s: i for i, s in enumerate(STATES)}


def sigmoid(x):
    return 1.0 / (1.0 + math.exp(-x))


def make_chain(memory=1.2, coupling=0.0, beta=2.0):
    P = np.zeros((4, 4), dtype=float)
    for i, (a, b) in enumerate(STATES):
        p_a_plus = sigmoid(2.0 * beta * (memory * a + coupling * b))
        p_b_plus = sigmoid(2.0 * beta * (memory * b + coupling * a))
        for j, (ap, bp) in enumerate(STATES):
            pa = p_a_plus if ap == 1 else (1.0 - p_a_plus)
            pb = p_b_plus if bp == 1 else (1.0 - p_b_plus)
            P[i, j] = pa * pb
    return P


def stationary(P):
    vals, vecs = np.linalg.eig(P.T)
    v = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    if v.sum() < 0:
        v = -v
    return v / v.sum()


def stationary_corr(pi):
    return float(sum(p * a * b for p, (a, b) in zip(pi, STATES)))


def mutual_information(pi):
    pa = {a: 0.0 for a in (-1, 1)}
    pb = {b: 0.0 for b in (-1, 1)}
    for p, (a, b) in zip(pi, STATES):
        pa[a] += p
        pb[b] += p
    out = 0.0
    for p, (a, b) in zip(pi, STATES):
        if p > 0:
            out += p * math.log2(p / (pa[a] * pb[b]))
    return out


def transition_prob_a_plus(P, a, b):
    i = INDEX[(a, b)]
    return float(sum(P[i, j] for j, (ap, bp) in enumerate(STATES) if ap == 1))


def transfer_information(P, pi):
    """Stationary observational I(A_{t+1}; B_t | A_t)."""
    p_a = {a: sum(pi[i] for i, (aa, bb) in enumerate(STATES) if aa == a) for a in (-1, 1)}
    p_b_given_a = {(b, a): pi[INDEX[(a, b)]] / p_a[a] for a in (-1, 1) for b in (-1, 1)}
    p_ap_given_ab = {}
    for a in (-1, 1):
        for b in (-1, 1):
            p_plus = transition_prob_a_plus(P, a, b)
            p_ap_given_ab[(1, a, b)] = p_plus
            p_ap_given_ab[(-1, a, b)] = 1.0 - p_plus
    p_ap_given_a = {}
    for a in (-1, 1):
        for ap in (-1, 1):
            p_ap_given_a[(ap, a)] = sum(
                p_b_given_a[(b, a)] * p_ap_given_ab[(ap, a, b)] for b in (-1, 1)
            )
    cmi = 0.0
    for i, (a, b) in enumerate(STATES):
        for ap in (-1, 1):
            p_cond = p_ap_given_ab[(ap, a, b)]
            p_base = p_ap_given_a[(ap, a)]
            p_joint = pi[i] * p_cond
            if p_joint > 0 and p_cond > 0 and p_base > 0:
                cmi += p_joint * math.log2(p_cond / p_base)
    return cmi


def causal_sensitivity(P):
    signed = []
    for a in (-1, 1):
        signed.append(transition_prob_a_plus(P, a, 1) - transition_prob_a_plus(P, a, -1))
    mean_signed = sum(signed) / len(signed)
    return abs(mean_signed), mean_signed


def relation_residence(P, pi, relation):
    inds = [i for i, (a, b) in enumerate(STATES) if a * b == relation]
    mass = pi[inds].sum()
    weights = pi[inds] / mass
    p_stay = sum(w * P[i, inds].sum() for w, i in zip(weights, inds))
    return 1.0 / (1.0 - p_stay)


def main():
    print("J corr MI transfer causal tau_anti tau_align")
    for J in np.linspace(-1.5, 1.5, 13):
        P = make_chain(coupling=float(J))
        pi = stationary(P)
        corr = stationary_corr(pi)
        mi = mutual_information(pi)
        transfer = transfer_information(P, pi)
        causal_mag, causal_signed = causal_sensitivity(P)
        tau_anti = relation_residence(P, pi, -1)
        tau_align = relation_residence(P, pi, 1)
        print(f"{J:+.2f} {corr:+.6f} {mi:.6f} {transfer:.6f} {causal_mag:.6f} {tau_anti:.3f} {tau_align:.3f}")


if __name__ == "__main__":
    main()
