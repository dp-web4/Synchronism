#!/usr/bin/env python3
"""Phase 8: relational identity through membership and topology turnover.

A society-like toy consists of 12 current members connected in a cycle graph. Members
carry binary states and locally repair toward an anti-aligned relational pattern:
neighbors prefer opposite states.

Three event classes occur:
1. relational repair/update of one current member;
2. membership replacement: remove one member identity, insert a brand-new one with
   random initial state;
3. topology rewiring: a 2-opt move changes current adjacency while preserving a single
   even cycle.

There are no semantically meaningful role labels in the measured macro-observable. The
identity diagnostic is permutation-invariant edge satisfaction: the fraction of current
relations that realize the anti-aligned organization.

For each episode we track the original member IDs and original edges until BOTH have
vanished. We then ask whether the relational organization still exists.

This is a formalism toy, not a physics model or evidence for Synchronism.
"""

import itertools
import math
import random
import numpy as np

N = 12
BETA = 3.0


def graph_edges(order):
    return {
        frozenset((order[i], order[(i + 1) % N]))
        for i in range(N)
    }


def edge_satisfaction(order, spins):
    """Fraction of current edges whose endpoint states are opposite."""
    return sum(
        spins[order[i]] != spins[order[(i + 1) % N]]
        for i in range(N)
    ) / N


def run_episode(p_replace, p_rewire, seed, max_steps=100_000):
    rng = random.Random(seed)

    # Initial organization is perfectly anti-aligned on the current graph.
    order = list(range(N))
    spins = {i: (1 if i % 2 == 0 else -1) for i in range(N)}

    initial_members = set(order)
    initial_edges = graph_edges(order)
    next_member_id = N

    member_turnover_time = None
    edge_turnover_time = None

    for t in range(1, max_steps + 1):
        u = rng.random()

        if u < p_replace:
            # Remove a component identity and insert a brand-new one. The newcomer
            # initially has no inherited state memory.
            idx = rng.randrange(N)
            old = order[idx]
            new = next_member_id
            next_member_id += 1
            order[idx] = new
            spins.pop(old, None)
            spins[new] = 1 if rng.random() < 0.5 else -1

        elif u < p_replace + p_rewire:
            # 2-opt rewiring: change two graph relations while keeping one even cycle.
            while True:
                i, j = sorted(rng.sample(range(N), 2))
                if j - i >= 2 and not (i == 0 and j == N - 1):
                    break
            order[i + 1:j + 1] = reversed(order[i + 1:j + 1])

        else:
            # Local relational repair: heat-bath update for antiferromagnetic coupling.
            idx = rng.randrange(N)
            node = order[idx]
            left = order[(idx - 1) % N]
            right = order[(idx + 1) % N]
            neighbor_sum = spins[left] + spins[right]
            field = -neighbor_sum
            p_plus = 1.0 / (1.0 + math.exp(-2.0 * BETA * field))
            spins[node] = 1 if rng.random() < p_plus else -1

        if member_turnover_time is None:
            if not (set(order) & initial_members):
                member_turnover_time = t

        if edge_turnover_time is None:
            if not (graph_edges(order) & initial_edges):
                edge_turnover_time = t

        if member_turnover_time is not None and edge_turnover_time is not None:
            return {
                "joint_turnover_time": t,
                "member_turnover_time": member_turnover_time,
                "edge_turnover_time": edge_turnover_time,
                "satisfaction": edge_satisfaction(order, spins),
            }

    raise RuntimeError("joint turnover not reached within max_steps")


def random_baseline(threshold=10 / 12):
    sats = []
    for bits in itertools.product((-1, 1), repeat=N):
        sat = sum(
            bits[i] != bits[(i + 1) % N]
            for i in range(N)
        ) / N
        sats.append(sat)
    sats = np.array(sats)
    return float(sats.mean()), float((sats >= threshold).mean())


def summarize(p_replace, p_rewire, episodes=1000, seed0=10_000):
    rows = [
        run_episode(p_replace, p_rewire, seed0 + i)
        for i in range(episodes)
    ]
    joint = np.array([r["joint_turnover_time"] for r in rows], dtype=float)
    sat = np.array([r["satisfaction"] for r in rows], dtype=float)
    return {
        "p_replace": p_replace,
        "p_rewire": p_rewire,
        "mean_joint_turnover": joint.mean(),
        "mean_satisfaction": sat.mean(),
        "fraction_ge_10_of_12": (sat >= 10 / 12).mean(),
        "fraction_perfect": (sat == 1.0).mean(),
    }


def main():
    baseline_mean, baseline_high = random_baseline()
    print("random-state baseline on a 12-cycle:")
    print(f"  mean edge satisfaction: {baseline_mean:.6f}")
    print(f"  P(satisfaction >= 10/12): {baseline_high:.6f}")

    print("\nturnover sweep (1000 episodes each):")
    print(
        "p_replace p_rewire mean_joint mean_sat "
        "P(sat>=10/12) P(perfect)"
    )
    for pr, pw in ((0.02, 0.001), (0.03, 0.003), (0.05, 0.003)):
        r = summarize(pr, pw)
        print(
            f"{pr:0.3f}     {pw:0.3f}    "
            f"{r['mean_joint_turnover']:9.3f}  "
            f"{r['mean_satisfaction']:.6f}    "
            f"{r['fraction_ge_10_of_12']:.4f}         "
            f"{r['fraction_perfect']:.4f}"
        )


if __name__ == "__main__":
    main()
