#!/usr/bin/env python3
"""Phase 2 toy model for the Markov coherence/emergence arc.

Build an 8-state hierarchical Markov chain with:
- four strong 2-state pairs,
- two weaker 4-state macro-groups,
- weak cross-group leakage,
then add a small deterministic perturbation.

Exhaustively evaluates all partitions with 2-4 blocks and no singletons.
Reports:
- geometric-mean dynamic coherence (exit/internal time-scale ratio),
- approximate strong-lumpability error,
- Markov Stability across observation times.

This is a formalism toy, not a physics model.
"""

import math
import numpy as np


def set_partitions(seq):
    if not seq:
        yield []
        return
    first = seq[0]
    for rest in set_partitions(seq[1:]):
        yield [[first]] + [b[:] for b in rest]
        for i in range(len(rest)):
            nr = [b[:] for b in rest]
            nr[i] = [first] + nr[i]
            yield nr


def canonical_partitions(n):
    seen = set()
    out = []
    for part in set_partitions(list(range(n))):
        blocks = tuple(sorted((tuple(sorted(b)) for b in part), key=lambda x: (x[0], len(x), x)))
        if blocks in seen:
            continue
        seen.add(blocks)
        if 2 <= len(blocks) <= 4 and all(len(b) >= 2 for b in blocks):
            out.append(blocks)
    return out


def make_chain(seed=42):
    n = 8
    pairs = [{0, 1}, {2, 3}, {4, 5}, {6, 7}]
    macro = [{0, 1, 2, 3}, {4, 5, 6, 7}]
    P = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            if i == j:
                P[i, j] = 0.50
            elif any(i in pair and j in pair for pair in pairs):
                P[i, j] = 0.30
            elif any(i in group and j in group for group in macro):
                P[i, j] = 0.09
            else:
                P[i, j] = 0.005

    rng = np.random.default_rng(seed)
    P = np.clip(P + rng.normal(0.0, 0.003, size=P.shape), 1e-5, None)
    return P / P.sum(axis=1, keepdims=True)


def stationary(P):
    vals, vecs = np.linalg.eig(P.T)
    pi = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    return pi / pi.sum()


def block_metrics(P, block):
    idx = list(block)
    Q = P[np.ix_(idx, idx)]

    rho = max(abs(np.linalg.eigvals(Q)))
    tau_exit = 1.0 / (1.0 - rho)

    row_stay = Q.sum(axis=1)
    W = Q / row_stay[:, None]
    eigs = sorted((abs(x) for x in np.linalg.eigvals(W)), reverse=True)
    lam2 = eigs[1]
    tau_internal = 1.0 / (1.0 - lam2)

    return tau_exit, tau_internal, tau_exit / tau_internal


def lumpability_error(P, part):
    err = 0.0
    for source in part:
        for target in part:
            probs = [P[i, list(target)].sum() for i in source]
            err = max(err, max(probs) - min(probs))
    return err


def markov_stability(P, pi, part, t):
    Pt = np.linalg.matrix_power(P, t)
    M = np.diag(pi) @ Pt - np.outer(pi, pi)
    return float(sum(M[np.ix_(list(block), list(block))].sum() for block in part))


def main():
    P = make_chain()
    pi = stationary(P)
    parts = canonical_partitions(8)

    rows = []
    for part in parts:
        metrics = [block_metrics(P, block) for block in part]
        coherence = math.exp(sum(math.log(m[2]) for m in metrics) / len(metrics))
        rows.append((part, coherence, lumpability_error(P, part)))

    print(f"candidate partitions: {len(rows)}")
    print("\nBest partition by Markov Stability at observation time t:")
    for t in [1, 2, 3, 5, 10, 20, 50]:
        best = max(rows, key=lambda row: markov_stability(P, pi, row[0], t))
        print(
            f"t={t:>2}: {best[0]}  "
            f"stability={markov_stability(P, pi, best[0], t):.6f}  "
            f"coherence={best[1]:.3f}  lump_err={best[2]:.5f}"
        )

    print("\nConstructed hierarchy:")
    targets = [
        ((0, 1), (2, 3), (4, 5), (6, 7)),
        ((0, 1, 2, 3), (4, 5, 6, 7)),
    ]
    for target in targets:
        row = next(row for row in rows if row[0] == target)
        print(f"{target}: coherence={row[1]:.6f}, lump_err={row[2]:.6f}")
        for block in target:
            te, ti, c = block_metrics(P, block)
            print(f"  {block}: tau_exit={te:.6f}, tau_internal={ti:.6f}, C={c:.6f}")


if __name__ == "__main__":
    main()
