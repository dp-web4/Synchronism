#!/usr/bin/env python3
"""Phase 6: discover slow relational variables from raw microstate trajectories.

The hidden construction contains the same rough hierarchy as Phase 5, but removes the
special analytic privilege that parity observables were exact eigenfunctions.

Changes from Phase 5:
- add a quenched random microstate energy;
- use Metropolis acceptance, making transition rates state-dependent/nonlinear;
- simulate a finite raw trajectory;
- estimate a reversible transition operator from transition counts;
- diagonalize that empirical operator with NO relational candidate features supplied.

Only after the slow modes are learned do we compare them to hidden diagnostic
relations Q, q1, q2, r1..r4.

This is a formalism toy, not a physics model or evidence for Synchronism.
"""

import math
import numpy as np

N = 8
NSTATES = 2 ** N
PAIRS = [(0, 1), (2, 3), (4, 5), (6, 7)]


def all_states():
    return np.array([
        [1 if (s >> i) & 1 else -1 for i in range(N)]
        for s in range(NSTATES)
    ], dtype=int)


def nuisance_energy(states, seed=7):
    """Weak random fields and couplings break exact translation symmetry."""
    rng = np.random.default_rng(seed)
    h = rng.normal(0.0, 0.15, N)
    J = rng.normal(0.0, 0.06, (N, N))
    J = np.triu(J, 1)
    J = J + J.T
    out = np.empty(len(states))
    for k, x in enumerate(states):
        out[k] = -h @ x - 0.5 * np.sum(J * np.outer(x, x))
    return out


def proposal_moves():
    """Symmetric proposal kernel containing a hidden relational hierarchy."""
    moves = []
    add = lambda mask, p, name: moves.append((tuple(sorted(mask)), p, name))

    add((), 0.09, "identity")

    # Fast component/orientation replacement inside each base pair.
    for pair in PAIRS:
        add(pair, 0.60 / 4.0, "pair_gauge")

    # Change pair relations inside quartet 1 while preserving q1.
    for i in PAIRS[0]:
        for j in PAIRS[1]:
            add((i, j), 0.14 / 4.0, "within_q1")

    # Same for quartet 2, at a slightly different rate to break degeneracy.
    for i in PAIRS[2]:
        for j in PAIRS[3]:
            add((i, j), 0.10 / 4.0, "within_q2")

    # Change both quartet relations while preserving their product Q.
    for i in range(4):
        for j in range(4, 8):
            add((i, j), 0.05 / 16.0, "across_quartet")

    # Rare moves destroy Q; asymmetric rates further break exact degeneracies.
    for i in range(4):
        add((i,), 0.006 / 4.0, "single_q1")
    for i in range(4, 8):
        add((i,), 0.014 / 4.0, "single_q2")

    total = sum(p for _, p, _ in moves)
    if abs(total - 1.0) > 1e-12:
        raise RuntimeError(f"proposal probabilities sum to {total}")
    return moves


def build_transition(states, energy, beta=0.7):
    """Metropolis chain with symmetric proposals and state-dependent acceptance."""
    P = np.zeros((NSTATES, NSTATES), dtype=float)
    for i in range(NSTATES):
        for mask, prob, _ in proposal_moves():
            j = i
            for bit in mask:
                j ^= 1 << bit
            if i == j:
                P[i, i] += prob
                continue
            accept = min(1.0, math.exp(-beta * (energy[j] - energy[i])))
            P[i, j] += prob * accept
            P[i, i] += prob * (1.0 - accept)
    return P


def simulate_counts(P, steps=700_000, seed=123):
    rng = np.random.default_rng(seed)
    cumulative = np.cumsum(P, axis=1)
    counts = np.zeros_like(P, dtype=np.int64)
    state = int(rng.integers(NSTATES))
    for _ in range(steps):
        u = rng.random()
        nxt = int(np.searchsorted(cumulative[state], u, side="right"))
        counts[state, nxt] += 1
        state = nxt
    return counts


def empirical_slow_modes(counts):
    """Reversible count estimator followed by symmetric eigendecomposition."""
    C = counts + counts.T
    degree = C.sum(axis=1).astype(float)
    S = C / np.sqrt(degree[:, None] * degree[None, :])
    eigenvalues, U = np.linalg.eigh(S)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    U = U[:, order]
    # Convert normalized symmetric eigenvectors into state functions.
    phi = U / np.sqrt(degree[:, None])
    weights = degree / degree.sum()
    return eigenvalues, phi, weights


def relation(states, indices):
    return np.prod(states[:, list(indices)], axis=1)


def weighted_corr(a, b, w):
    am = np.sum(w * a)
    bm = np.sum(w * b)
    aa = a - am
    bb = b - bm
    return np.sum(w * aa * bb) / math.sqrt(
        np.sum(w * aa * aa) * np.sum(w * bb * bb)
    )


def main():
    states = all_states()
    energy = nuisance_energy(states)
    P = build_transition(states, energy)
    counts = simulate_counts(P)
    eigenvalues, phi, weights = empirical_slow_modes(counts)

    # These hidden variables are diagnostics ONLY. They are never passed into the
    # transition estimator or eigensolver.
    hidden = {
        "Q": relation(states, range(8)),
        "q1": relation(states, range(4)),
        "q2": relation(states, range(4, 8)),
        "r1": relation(states, (0, 1)),
        "r2": relation(states, (2, 3)),
        "r3": relation(states, (4, 5)),
        "r4": relation(states, (6, 7)),
    }

    print(f"raw transitions observed: {counts.sum()}")
    print(f"microstates visited: {(counts.sum(axis=1) > 0).sum()} / {NSTATES}")
    print("\nLeading learned nonconstant slow modes:")
    print("rank   lambda      tau_corr    best hidden diagnostic correlation")

    for rank in range(1, 12):
        lam = float(eigenvalues[rank])
        tau = -1.0 / math.log(abs(lam))
        scores = []
        for name, target in hidden.items():
            c = weighted_corr(phi[:, rank], target, weights)
            scores.append((abs(c), c, name))
        scores.sort(reverse=True)
        _, corr, name = scores[0]
        print(
            f"{rank:>4}   {lam:.6f}   {tau:9.4f}    "
            f"{name:>2s}: corr={corr:+.6f}"
        )


if __name__ == "__main__":
    main()
