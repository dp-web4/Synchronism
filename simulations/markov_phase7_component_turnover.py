#!/usr/bin/env python3
"""Phase 7: component turnover while larger relational identity persists.

Eight role positions carry replaceable components. The component at each role has a
binary state x_i. The role topology supports a staggered collective relation: adjacent
roles prefer opposite states, equivalently the transformed variables y_i=r_i*x_i
prefer common alignment for fixed role signs r_i in {+1,-1}.

Two kinds of events occur:
1. repair/update events: one occupant updates from the current collective relation;
2. turnover events: one occupant is replaced by a brand-new component identity and
   its state is initialized randomly.

The learner never receives component IDs or the staggered-order diagnostic. It sees
only raw 8-bit state transitions, learns the slowest empirical Markov mode, and only
then is that mode compared with the hidden relational order parameter.

Separate turnover episodes track explicit component identities to measure how long it
takes until every original component has been replaced.

This is a formalism toy, not a physics model or evidence for Synchronism.
"""

import math
import numpy as np

N = 8
NSTATES = 2 ** N
ROLE_SIGN = np.array([1, -1, 1, -1, 1, -1, 1, -1], dtype=int)
J = 1.0
BETA = 3.0
TURNOVER = 0.08

RNG_FIELDS = np.random.default_rng(42)
LOCAL_FIELDS = RNG_FIELDS.normal(0.0, 0.08, N)

STATES = np.array(
    [[1 if (s >> i) & 1 else -1 for i in range(N)] for s in range(NSTATES)],
    dtype=int,
)


def state_index(x):
    return int(sum((1 if x[i] > 0 else 0) << i for i in range(N)))


def staggered_order(x):
    return float(np.mean(ROLE_SIGN * x))


def heatbath_probability_plus(x, i, beta=BETA):
    stagger_sum = float(np.sum(ROLE_SIGN * x))
    collective = ROLE_SIGN[i] * J * (
        stagger_sum - ROLE_SIGN[i] * x[i]
    ) / (N - 1)
    local = collective + LOCAL_FIELDS[i]
    return 1.0 / (1.0 + math.exp(-2.0 * beta * local))


def build_exact_transition(turnover=TURNOVER, beta=BETA):
    P = np.zeros((NSTATES, NSTATES), dtype=float)
    for si, x in enumerate(STATES):
        # Replacement: brand-new component state is initially random.
        for i in range(N):
            for new_spin in (-1, 1):
                y = x.copy()
                y[i] = new_spin
                P[si, state_index(y)] += turnover * (1.0 / N) * 0.5

        # Relational repair/update.
        for i in range(N):
            p_plus = heatbath_probability_plus(x, i, beta=beta)
            for new_spin, prob in ((1, p_plus), (-1, 1.0 - p_plus)):
                y = x.copy()
                y[i] = new_spin
                P[si, state_index(y)] += (1.0 - turnover) * (1.0 / N) * prob
    return P


def stochastic_step(x, rng, turnover=TURNOVER, beta=BETA):
    x = x.copy()
    if rng.random() < turnover:
        i = int(rng.integers(N))
        x[i] = 1 if rng.random() < 0.5 else -1
        return x, True, i

    i = int(rng.integers(N))
    p_plus = heatbath_probability_plus(x, i, beta=beta)
    x[i] = 1 if rng.random() < p_plus else -1
    return x, False, i


def learn_from_raw_transitions(steps=500_000, seed=123):
    rng = np.random.default_rng(seed)
    counts = np.zeros((NSTATES, NSTATES), dtype=np.int64)
    x = STATES[int(rng.integers(NSTATES))].copy()

    for _ in range(steps):
        i = state_index(x)
        y, _, _ = stochastic_step(x, rng)
        j = state_index(y)
        counts[i, j] += 1
        x = y

    # Reversible count estimator. No relational features are supplied.
    C = counts + counts.T
    degree = C.sum(axis=1).astype(float)
    active = np.where(degree > 0)[0]
    S = C[np.ix_(active, active)] / np.sqrt(
        degree[active, None] * degree[None, active]
    )

    eigenvalues, U = np.linalg.eigh(S)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    U = U[:, order]

    phi = np.full((NSTATES, len(active)), np.nan)
    phi[active, :] = U / np.sqrt(degree[active, None])

    weights = np.zeros(NSTATES)
    weights[active] = degree[active] / degree[active].sum()
    return counts, eigenvalues, phi, weights


def weighted_corr(a, b, w):
    mask = np.isfinite(a) & (w > 0)
    a = a[mask]
    b = b[mask]
    w = w[mask]
    w = w / w.sum()
    am = np.sum(w * a)
    bm = np.sum(w * b)
    aa = a - am
    bb = b - bm
    return float(
        np.sum(w * aa * bb)
        / math.sqrt(np.sum(w * aa * aa) * np.sum(w * bb * bb))
    )


def complete_turnover_episode(mode, seed, max_steps=200_000):
    rng = np.random.default_rng(seed)

    # Start in one relational basin. Component IDs are arbitrary occupants of roles.
    x = ROLE_SIGN.copy()
    ids = np.arange(N, dtype=np.int64)
    originals = set(ids.tolist())
    next_id = N

    start_sign = 1 if mode[state_index(x)] >= 0 else -1
    first_mode_flip = None
    replacements = 0

    for t in range(1, max_steps + 1):
        if rng.random() < TURNOVER:
            i = int(rng.integers(N))
            ids[i] = next_id
            next_id += 1
            replacements += 1
            x[i] = 1 if rng.random() < 0.5 else -1
        else:
            i = int(rng.integers(N))
            p_plus = heatbath_probability_plus(x, i)
            x[i] = 1 if rng.random() < p_plus else -1

        current_sign = 1 if mode[state_index(x)] >= 0 else -1
        if first_mode_flip is None and current_sign != start_sign:
            first_mode_flip = t

        if not any(component_id in originals for component_id in ids):
            return {
                "steps": t,
                "same_learned_identity": current_sign == start_sign,
                "first_mode_flip": first_mode_flip,
                "replacements": replacements,
                "staggered_order": staggered_order(x),
            }

    raise RuntimeError("complete turnover did not occur within max_steps")


def exact_turnover_sweep():
    H_N = sum(1.0 / k for k in range(1, N + 1))
    out = []
    for q in (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.25, 0.30, 0.40):
        P = build_exact_transition(turnover=q)
        eigenvalues = np.linalg.eigvals(P)
        eigenvalues = sorted(np.real(eigenvalues), reverse=True)
        lam = eigenvalues[1]
        tau_identity = -1.0 / math.log(lam)

        # Coupon-collector mean: replacement events arrive with probability q per
        # microstep and choose one of N role occupants uniformly.
        tau_ancestry = N * H_N / q
        out.append((q, tau_identity, tau_ancestry, tau_identity / tau_ancestry))
    return out


def main():
    counts, eigenvalues, phi, weights = learn_from_raw_transitions()
    diagnostic = np.array([staggered_order(x) for x in STATES])

    learned = phi[:, 1].copy()
    corr = weighted_corr(learned, diagnostic, weights)
    if corr < 0:
        learned = -learned
        corr = -corr

    lam = float(eigenvalues[1])
    tau_learned = -1.0 / math.log(lam)

    print(f"raw transitions observed: {counts.sum()}")
    print(f"microstates visited: {(counts.sum(axis=1) > 0).sum()} / {NSTATES}")
    print(f"learned slow-mode eigenvalue: {lam:.8f}")
    print(f"learned slow-mode correlation time: {tau_learned:.3f}")
    print(f"post-hoc corr(mode, hidden staggered relation): {corr:.6f}")

    episodes = [complete_turnover_episode(learned, 1000 + k) for k in range(2000)]
    times = np.array([e["steps"] for e in episodes], dtype=float)
    same = np.array([e["same_learned_identity"] for e in episodes], dtype=float)
    replacements = np.array([e["replacements"] for e in episodes], dtype=float)

    print("\nComplete component-turnover episodes:")
    print(f"mean time until all original occupants are gone: {times.mean():.3f}")
    print(f"median time until all original occupants are gone: {np.median(times):.3f}")
    print(f"mean replacement events required: {replacements.mean():.3f}")
    print(f"fraction retaining same learned macro-identity at full turnover: {same.mean():.4f}")
    print(f"tau_identity / mean full-turnover time: {tau_learned / times.mean():.3f}")

    print("\nExact turnover-rate sweep:")
    print("q      tau_identity   tau_ancestry   R_turn=tau_identity/tau_ancestry")
    for q, ti, ta, ratio in exact_turnover_sweep():
        print(f"{q:0.2f}   {ti:11.3f}   {ta:12.3f}   {ratio:8.3f}")


if __name__ == "__main__":
    main()
