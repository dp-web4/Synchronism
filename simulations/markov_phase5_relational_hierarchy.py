#!/usr/bin/env python3
"""Phase 5: recursive relational coarse-graining.

Eight binary microvariables are grouped only for construction into four base pairs,
two quartets, and one full system. The transition process contains four classes of
bit-flip moves:

1. fast pair-gauge moves: flip both members of a base pair;
2. intermediate within-quartet moves: flip one member in each of two base pairs
   in the same quartet;
3. slower across-quartet moves: flip one member in each half;
4. rare single flips.

The intended relational observables are therefore:
  r1=s1*s2, r2=s3*s4, r3=s5*s6, r4=s7*s8
  q1=r1*r2=s1*s2*s3*s4
  q2=r3*r4=s5*s6*s7*s8
  Q=q1*q2=product(s1..s8)

But the analysis is NOT given those observables. It exhaustively evaluates every
nonconstant parity observable (255 total), computes its exact one-step eigenvalue
under the translation-invariant Markov process, and ranks the slow modes.

This is a formalism toy, not a physics model or evidence for Synchronism.
"""

from itertools import combinations
from math import log

N = 8


def build_moves():
    """Return (flip_mask, probability, class_name) tuples.

    Probabilities sum to one. A flip mask is represented as a frozenset of zero-based
    variable indices. The hierarchy is deliberately present in the transition-channel
    structure; the question is whether blind spectral search recovers it.
    """
    moves = [(frozenset(), 0.07, "identity")]

    pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]

    # Relation-preserving orientation changes at the base-pair level.
    for pair in pairs:
        moves.append((frozenset(pair), 0.48 / 4.0, "pair_gauge"))

    # Flip both base relations inside one quartet, preserving the quartet relation.
    for pair_a, pair_b in ((pairs[0], pairs[1]), (pairs[2], pairs[3])):
        for i in pair_a:
            for j in pair_b:
                moves.append((frozenset((i, j)), 0.24 / 8.0, "within_quartet"))

    # Flip one base relation in each quartet, preserving the full-system relation.
    for i in range(4):
        for j in range(4, 8):
            moves.append((frozenset((i, j)), 0.16 / 16.0, "across_quartet"))

    # Rare move that breaks every higher parity invariant containing that variable.
    for i in range(N):
        moves.append((frozenset((i,)), 0.05 / 8.0, "single"))

    total = sum(prob for _, prob, _ in moves)
    if abs(total - 1.0) > 1e-12:
        raise RuntimeError(f"move probabilities sum to {total}, not 1")
    return moves


def parity_eigenvalue(observable, moves):
    """Exact eigenvalue for parity observable prod_{i in observable} s_i.

    A move reverses the parity iff its flip mask overlaps the observable in an odd
    number of variables. For this random-walk process on Z_2^N, every parity
    observable is an exact eigenfunction.
    """
    p_reverse = sum(
        prob
        for mask, prob, _ in moves
        if len(observable.intersection(mask)) % 2 == 1
    )
    return 1.0 - 2.0 * p_reverse


def relaxation_time(lam):
    """Exponential autocorrelation time -1/log(|lambda|)."""
    a = abs(lam)
    if a == 0.0:
        return 0.0
    if a >= 1.0:
        return float("inf")
    return -1.0 / log(a)


def label(obs):
    return "*".join(f"s{i + 1}" for i in sorted(obs))


def main():
    moves = build_moves()
    rows = []

    for size in range(1, N + 1):
        for combo in combinations(range(N), size):
            obs = frozenset(combo)
            lam = parity_eigenvalue(obs, moves)
            rows.append((lam, relaxation_time(lam), obs))

    rows.sort(key=lambda row: row[0], reverse=True)

    print(f"nonconstant parity observables searched: {len(rows)}")
    print("\nTop 20 slow positive modes:")
    print("rank  lambda      tau_corr    observable")
    for rank, (lam, tau, obs) in enumerate(rows[:20], start=1):
        print(f"{rank:>4}  {lam: .6f}   {tau:9.5f}   {label(obs)}")

    expected = {
        "Q (relation among quartet-relations)": frozenset(range(8)),
        "q1 (relation among r1,r2)": frozenset(range(4)),
        "q2 (relation among r3,r4)": frozenset(range(4, 8)),
        "r1": frozenset((0, 1)),
        "r2": frozenset((2, 3)),
        "r3": frozenset((4, 5)),
        "r4": frozenset((6, 7)),
        "s1 (representative microvariable)": frozenset((0,)),
    }

    print("\nConstructed relational hierarchy (reported only after blind ranking):")
    for name, obs in expected.items():
        lam = parity_eigenvalue(obs, moves)
        print(
            f"{name:38s} lambda={lam:.6f}  "
            f"tau_corr={relaxation_time(lam):.6f}  observable={label(obs)}"
        )


if __name__ == "__main__":
    main()
