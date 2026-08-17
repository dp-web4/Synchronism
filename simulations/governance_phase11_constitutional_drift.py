#!/usr/bin/env python3
"""Phase 11: locally valid amendments can accumulate into large constitutional drift.

A toy constitution is a binary decision vector over 64 governance contexts. An amendment
changes exactly one context, so every adjacent constitution differs by only 1/64 = 1.5625%.

We compare:
- local amendment distance,
- endpoint distance from the constitution a relying party originally trusted,
- MRH/stakes-weighted distance,
- cumulative path length,
- maximum historical excursion.

The point is to show that witnessed/authorized lineage does not by itself justify indefinite
trust transfer. Constitutional compatibility is a separate, relying-party-relative question.

This is a governance formalism toy, not a proposed production trust formula.
"""

import numpy as np

N = 64

# Four contexts are treated as high-stakes in this relying party's MRH.
WEIGHTS = np.full(N, 0.60 / 60.0)
WEIGHTS[:4] = 0.10
assert abs(WEIGHTS.sum() - 1.0) < 1e-12


def unweighted_distance(a, b):
    return float(np.mean(a != b))


def weighted_distance(a, b):
    return float(np.sum(WEIGHTS[a != b]))


def chain_diagnostics(chain):
    origin = chain[0]
    local = [weighted_distance(chain[i], chain[i + 1]) for i in range(len(chain) - 1)]
    endpoint = weighted_distance(origin, chain[-1])
    path_length = sum(local)
    max_excursion = max(weighted_distance(origin, c) for c in chain)
    return endpoint, path_length, max_excursion, local


def monotone_drift(steps=32):
    chain = [np.zeros(N, dtype=np.int8)]
    for i in range(steps):
        nxt = chain[-1].copy()
        nxt[i] ^= 1
        chain.append(nxt)
    return chain


def out_and_back(steps=32):
    chain = monotone_drift(steps)
    for i in reversed(range(steps)):
        nxt = chain[-1].copy()
        nxt[i] ^= 1
        chain.append(nxt)
    return chain


def main():
    monotone = monotone_drift(32)
    endpoint_w, path_w, excursion_w, local_w = chain_diagnostics(monotone)

    print("Monotone 32-amendment drift")
    print(f"adjacent unweighted change: {1/N:.6f}")
    print(f"largest adjacent weighted change: {max(local_w):.6f}")
    print(f"endpoint unweighted distance: {unweighted_distance(monotone[0], monotone[-1]):.6f}")
    print(f"endpoint MRH-weighted distance: {endpoint_w:.6f}")
    print(f"cumulative weighted path length: {path_w:.6f}")
    print(f"maximum weighted excursion: {excursion_w:.6f}")

    returned = out_and_back(32)
    endpoint_w, path_w, excursion_w, _ = chain_diagnostics(returned)
    print("\nOut-and-back amendment history")
    print(f"endpoint unweighted distance: {unweighted_distance(returned[0], returned[-1]):.6f}")
    print(f"endpoint MRH-weighted distance: {endpoint_w:.6f}")
    print(f"cumulative weighted path length: {path_w:.6f}")
    print(f"maximum weighted excursion: {excursion_w:.6f}")

    origin = np.zeros(N, dtype=np.int8)
    high = origin.copy()
    high[:4] = 1
    low = origin.copy()
    low[4:8] = 1

    print("\nSame raw edit count, different relying-party relevance")
    print(f"4 high-stakes edits: unweighted={unweighted_distance(origin, high):.6f}, weighted={weighted_distance(origin, high):.6f}")
    print(f"4 low-stakes edits:  unweighted={unweighted_distance(origin, low):.6f}, weighted={weighted_distance(origin, low):.6f}")


if __name__ == "__main__":
    main()
