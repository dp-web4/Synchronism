#!/usr/bin/env python3
"""Phase 9: constitutional continuity under governance-mechanism replacement.

A six-member ring supports one macro relational organization: neighboring states prefer
anti-alignment. Three deliberately different update kernels all stabilize that same
macro relation:

1. local heat-bath repair,
2. edge arbitration,
3. noisy best response.

The active governance kernel itself is a Markov state and can be replaced at any step.
We build the exact joint transition matrix over (microstate, active kernel), then measure:
- stationary relational satisfaction,
- lifetime of the slow macro mode,
- correlation of that learned slow mode with hidden staggered identity,
- how many kernel tenures fit inside one macro-identity lifetime.

A fourth, constitution-incompatible kernel favors alignment rather than anti-alignment.
Replacing one member of the admissible kernel family with that incompatible mechanism is
used as a falsification/control arm.

The point is not that anti-alignment is governance. It is a transparent stand-in for a
higher-scale invariant implemented by multiple lower-scale procedures.

This is a formalism/governance toy, not a physics model or evidence for Synchronism.
"""

import math
import numpy as np

N = 6
NSTATES = 2 ** N
STATES = np.array(
    [[1 if (s >> i) & 1 else -1 for i in range(N)] for s in range(NSTATES)],
    dtype=int,
)
ROLE_SIGN = np.array([1 if i % 2 == 0 else -1 for i in range(N)], dtype=int)


def state_index(x):
    return int(sum((1 if x[i] > 0 else 0) << i for i in range(N)))


def satisfaction(x):
    return sum(x[i] != x[(i + 1) % N] for i in range(N)) / N


def staggered_identity(x):
    return float(np.mean(ROLE_SIGN * x))


def heatbath_kernel(beta=2.4):
    P = np.zeros((NSTATES, NSTATES), dtype=float)
    for si, x in enumerate(STATES):
        for i in range(N):
            left = x[(i - 1) % N]
            right = x[(i + 1) % N]
            field = -(left + right) / 2.0
            p_plus = 1.0 / (1.0 + math.exp(-2.0 * beta * field))
            for value, p in ((1, p_plus), (-1, 1.0 - p_plus)):
                y = x.copy()
                y[i] = value
                P[si, state_index(y)] += p / N
    return P


def edge_arbitration_kernel(fix=0.96, noise=0.02):
    P = np.zeros((NSTATES, NSTATES), dtype=float)
    for si, x in enumerate(STATES):
        for i in range(N):
            u, v = i, (i + 1) % N
            if x[u] == x[v]:
                for target in (u, v):
                    y = x.copy()
                    y[target] *= -1
                    P[si, state_index(y)] += (0.5 * fix) / N
                P[si, si] += (1.0 - fix) / N
            else:
                for target in (u, v):
                    y = x.copy()
                    y[target] *= -1
                    P[si, state_index(y)] += (0.5 * noise) / N
                P[si, si] += (1.0 - noise) / N
    return P


def best_response_kernel(p_best=0.93, align=False):
    P = np.zeros((NSTATES, NSTATES), dtype=float)
    for si, x in enumerate(STATES):
        for i in range(N):
            left = x[(i - 1) % N]
            right = x[(i + 1) % N]
            scores = {}
            for value in (-1, 1):
                if align:
                    scores[value] = int(value == left) + int(value == right)
                else:
                    scores[value] = int(value != left) + int(value != right)

            if scores[1] > scores[-1]:
                probs = {1: p_best, -1: 1.0 - p_best}
            elif scores[-1] > scores[1]:
                probs = {-1: p_best, 1: 1.0 - p_best}
            else:
                probs = {1: 0.5, -1: 0.5}

            for value, p in probs.items():
                y = x.copy()
                y[i] = value
                P[si, state_index(y)] += p / N
    return P


def joint_chain(kernels, switch_probability):
    """Exact chain over (active governance kernel, microstate)."""
    m = len(kernels)
    P = np.zeros((m * NSTATES, m * NSTATES), dtype=float)

    for current_kernel, K in enumerate(kernels):
        for si in range(NSTATES):
            for sj, p_state in enumerate(K[si]):
                if p_state == 0.0:
                    continue
                for next_kernel in range(m):
                    if next_kernel == current_kernel:
                        p_kernel = 1.0 - switch_probability
                    else:
                        p_kernel = switch_probability / (m - 1)
                    P[
                        current_kernel * NSTATES + si,
                        next_kernel * NSTATES + sj,
                    ] += p_state * p_kernel
    return P


def stationary(P):
    values, vectors = np.linalg.eig(P.T)
    k = int(np.argmin(np.abs(values - 1.0)))
    pi = np.real(vectors[:, k])
    if pi.sum() < 0:
        pi = -pi
    pi = np.maximum(pi, 0.0)
    return pi / pi.sum()


def weighted_corr(a, b, weights):
    am = float(weights @ a)
    bm = float(weights @ b)
    aa = a - am
    bb = b - bm
    return float(
        (weights @ (aa * bb))
        / math.sqrt((weights @ (aa * aa)) * (weights @ (bb * bb)))
    )


def analyze(kernels, switch_probability):
    P = joint_chain(kernels, switch_probability)
    pi = stationary(P)
    m = len(kernels)

    sat = np.tile(np.array([satisfaction(x) for x in STATES]), m)
    hidden_identity = np.tile(
        np.array([staggered_identity(x) for x in STATES]), m
    )

    values, vectors = np.linalg.eig(P)
    order = np.argsort(np.abs(values))[::-1]
    values = values[order]
    vectors = vectors[:, order]

    lam = values[1]
    tau_identity = -1.0 / math.log(abs(lam))
    slow_mode = np.real(vectors[:, 1])
    corr = abs(weighted_corr(slow_mode, hidden_identity, pi))

    mean_tenure = 1.0 / switch_probability
    return {
        "mean_satisfaction": float(pi @ sat),
        "lambda_identity": float(np.real(lam)),
        "tau_identity": tau_identity,
        "hidden_corr": corr,
        "mean_kernel_tenure": mean_tenure,
        "identity_per_tenure": tau_identity / mean_tenure,
    }


def mean_satisfaction_drift(kernel):
    """Uniform-state average E[S_{t+1}-S_t | S_t=s] by satisfaction level."""
    current = np.array([satisfaction(x) for x in STATES])
    expected_next = kernel @ current
    out = {}
    for s in sorted(set(current)):
        inds = np.where(np.isclose(current, s))[0]
        out[s] = float(np.mean(expected_next[inds] - s))
    return out


def main():
    heat = heatbath_kernel()
    edge = edge_arbitration_kernel()
    best = best_response_kernel()
    incompatible = best_response_kernel(align=True)

    good = [heat, edge, best]
    mixed = [heat, edge, incompatible]

    print("Constitution-compatible kernel-family sweep")
    print(
        "switch   mean_S   tau_identity   mean_tenure   identity/tenure   "
        "corr(slow,hidden)"
    )
    for switch in (0.01, 0.03, 0.08, 0.20, 0.50):
        r = analyze(good, switch)
        print(
            f"{switch:0.2f}     {r['mean_satisfaction']:.6f}   "
            f"{r['tau_identity']:11.3f}   {r['mean_kernel_tenure']:11.3f}   "
            f"{r['identity_per_tenure']:15.3f}   {r['hidden_corr']:.6f}"
        )

    print("\nControl: one constitution-incompatible kernel in the family")
    print("switch   mean_S   tau_identity   corr(slow,hidden)")
    for switch in (0.03, 0.08, 0.20):
        r = analyze(mixed, switch)
        print(
            f"{switch:0.2f}     {r['mean_satisfaction']:.6f}   "
            f"{r['tau_identity']:11.3f}   {r['hidden_corr']:.6f}"
        )

    print("\nOne-step constitutional drift E[Delta satisfaction | current satisfaction]")
    for name, kernel in (
        ("heatbath", heat),
        ("edge_arbitration", edge),
        ("best_response", best),
        ("incompatible_alignment", incompatible),
    ):
        drift = mean_satisfaction_drift(kernel)
        formatted = ", ".join(f"S={s:.3f}:{d:+.4f}" for s, d in drift.items())
        print(f"{name:24s} {formatted}")


if __name__ == "__main__":
    main()
