#!/usr/bin/env python3
"""SA-3A exact finite-world controls (Codex, 2026-09-08).

Registration: explorations/2026-09-08-codex-sa3a-loss-and-identifiability-charter.md
Run from the repo root: python3 simulations/mrh_sa3a_loss_identifiability.py

Stdlib only. No sampling, fitting, network access, or file writes. JSON goes to
stdout; a failed acceptance check raises and exits nonzero. All information is
in nats. Worlds are joint tables over (retained X, excluded Z, target Y).
"""

import hashlib
import itertools
import json
import math
from pathlib import Path


TOL = 1e-12


def validate(table, width):
    """Reject malformed probability tables; do not repair or normalize them."""
    if not table:
        raise ValueError("empty probability table")
    for state, probability in table.items():
        if len(state) != width:
            raise ValueError("wrong state dimension")
        if not math.isfinite(probability) or probability < 0:
            raise ValueError("invalid probability")
    if abs(math.fsum(table.values()) - 1) > TOL:
        raise ValueError("probability mass must equal one")


def marginal(table, axes):
    result = {}
    for state, probability in table.items():
        key = tuple(state[i] for i in axes)
        result[key] = result.get(key, 0.0) + probability
    return result


def log_risk(table, observed):
    """Optimal expected -log P(Y | observed), computed from predictive rows."""
    joint = marginal(table, tuple(observed) + (2,))
    conditions = marginal(table, observed)
    return math.fsum(
        -p * math.log(p / conditions[key[:-1]])
        for key, p in joint.items() if p > 0
    )


def square_risk(table, observed):
    """Optimal expected square loss around each conditional mean."""
    conditions = marginal(table, observed)
    weighted = {}
    for state, p in table.items():
        key = tuple(state[i] for i in observed)
        weighted[key] = weighted.get(key, 0.0) + p * state[2]
    means = {key: weighted[key] / p for key, p in conditions.items() if p > 0}
    return math.fsum(
        p * (state[2] - means[tuple(state[i] for i in observed)]) ** 2
        for state, p in table.items() if p > 0
    )


def cmi(table):
    """I(Y;Z|X), directly via the conditional probability ratio (not risk gap)."""
    px = marginal(table, (0,))
    pxz = marginal(table, (0, 1))
    pxy = marginal(table, (0, 2))
    return math.fsum(
        p * math.log(p * px[(x,)] / (pxz[(x, z)] * pxy[(x, y)]))
        for (x, z, y), p in table.items() if p > 0
    )


def scores(table):
    validate(table, 3)
    return {
        "cmi_nats": cmi(table),
        "log_risk_restricted": log_risk(table, (0,)),
        "log_risk_extended": log_risk(table, (0, 1)),
        "square_risk_restricted": square_risk(table, (0,)),
        "square_risk_extended": square_risk(table, (0, 1)),
    }


def observer_summary(pxy):
    """Only sees P(X,Y). Does not accept Z, a world ID, or hidden truth.

    For unrestricted finite hidden extensions, [0,H(Y|X)] is the sharp CMI
    interval: independent Z attains zero; Z=Y attains the upper endpoint.
    This is a population identifiability bound, not a learned classifier or a
    finite-data confidence interval.
    """
    validate(pxy, 2)
    px = marginal(pxy, (0,))
    entropy = math.fsum(
        -p * math.log(p / px[(x,)]) for (x, _), p in pxy.items() if p > 0
    )
    return {
        "status": "not_identifiable" if entropy > TOL else "identified_zero",
        "exclusion_log_loss_interval_nats": [0.0, entropy],
        "condition": "known P(X,Y); unrestricted finite hidden extensions",
    }


def add(table, state, mass):
    table[state] = table.get(state, 0.0) + mass


def binary_entropy(p):
    return -p * math.log(p) - (1 - p) * math.log(1 - p)


def suite():
    checks = []

    def near(label, actual, expected):
        if not math.isfinite(actual) or abs(actual - expected) > TOL:
            raise AssertionError(f"{label}: {actual!r} != {expected!r}")
        checks.append(label)

    def require(label, condition):
        if not condition:
            raise AssertionError(label)
        checks.append(label)

    def measured(label, table):
        result = scores(table)
        near(label + ": log-loss gap equals CMI",
             result["log_risk_restricted"] - result["log_risk_extended"],
             result["cmi_nats"])
        require(label + ": CMI nonnegative", result["cmi_nats"] >= -TOL)
        return result

    # Probability-instrument negative controls are separate from scientific cases.
    for label, bad in (
        ("empty", {}),
        ("unnormalized", {(0, 0, 0): 0.5}),
        ("negative", {(0, 0, 0): -0.1, (0, 0, 1): 1.1}),
        ("nonfinite", {(0, 0, 0): float("nan")}),
        ("wrong_dimension", {(0, 0): 1.0}),
    ):
        try:
            validate(bad, 3)
        except ValueError:
            checks.append("table rejects " + label)
        else:
            raise AssertionError("accepted malformed table: " + label)

    p = 0.1
    hp = binary_entropy(p)
    gain = math.log(2) - hp
    results = {}

    closed = {}
    for x, z, noise in itertools.product((0, 1), repeat=3):
        add(closed, (x, z, x ^ noise), 0.25 * (p if noise else 1 - p))
    r = measured("T1", closed)
    near("T1: no exclusion information", r["cmi_nats"], 0.0)
    near("T1: irreducible log risk", r["log_risk_restricted"], hp)
    near("T1: restricted MSE", r["square_risk_restricted"], 0.09)
    near("T1: extended MSE", r["square_risk_extended"], 0.09)
    results["T1_closed_but_noisy"] = r

    variance_only = {(0, z, sign * (z + 1)): 0.25
                     for z in (0, 1) for sign in (-1, 1)}
    r = measured("T2", variance_only)
    near("T2: magnitude reveals bit", r["cmi_nats"], math.log(2))
    near("T2: restricted MSE", r["square_risk_restricted"], 2.5)
    near("T2: extended MSE", r["square_risk_extended"], 2.5)
    r["mse_ratio"] = r["square_risk_restricted"] / r["square_risk_extended"]
    r["gaussian_formula_if_misapplied"] = math.exp(2 * r["cmi_nats"])
    near("T2: MSE ratio one", r["mse_ratio"], 1.0)
    near("T2: Gaussian ratio would be four", r["gaussian_formula_if_misapplied"], 4.0)
    results["T2_distribution_not_mean"] = r

    mean_channel = {}
    for z, noise in itertools.product((0, 1), repeat=2):
        add(mean_channel, (0, z, z ^ noise), 0.5 * (p if noise else 1 - p))
    r = measured("T3", mean_channel)
    scaled = {(x, z, 10 * y): mass for (x, z, y), mass in mean_channel.items()}
    rs = measured("T3 scaled", scaled)
    near("T3: CMI", r["cmi_nats"], gain)
    near("T3: restricted MSE", r["square_risk_restricted"], 0.25)
    near("T3: extended MSE", r["square_risk_extended"], 0.09)
    near("T3: scaled information invariant", rs["cmi_nats"], r["cmi_nats"])
    for key in ("square_risk_restricted", "square_risk_extended"):
        near("T3: scaling " + key, rs[key], 100 * r[key])
    results["T3_mean_and_scale"] = {"original": r, "scaled_by_10": rs}

    pairs = list(itertools.product((0, 1), repeat=2))
    transition = [[0.0 for _ in pairs] for _ in pairs]
    history = {}
    for i, (previous, current) in enumerate(pairs):
        for noise, probability in ((0, 1 - p), (1, p)):
            future = previous ^ noise
            j = pairs.index((current, future))
            transition[i][j] += probability
            add(history, (current, previous, future), 0.25 * probability)
        near(f"T4: transition row {i}", sum(transition[i]), 1.0)
    for j in range(4):
        near(f"T4: stationary pair {j}", sum(0.25 * row[j] for row in transition), 0.25)
    r = measured("T4", history)
    near("T4: memory information", r["cmi_nats"], gain)
    near("T4: present-only MSE", r["square_risk_restricted"], 0.25)
    near("T4: history-aware MSE", r["square_risk_extended"], 0.09)
    results["T4_memory"] = {"scores": r, "pair_transition_matrix": transition}

    direct, hidden = {}, {}
    for x, z in itertools.product((0, 1), repeat=2):
        pz = p if z else 1 - p
        add(hidden, (x, z, x ^ z), 0.5 * pz)
        for noise in (0, 1):
            add(direct, (x, z, x ^ noise), 0.5 * pz * (p if noise else 1 - p))
    rd, rh = measured("T5 direct", direct), measured("T5 hidden", hidden)
    od, oh = marginal(direct, (0, 2)), marginal(hidden, (0, 2))
    require("T5: same observable support", od.keys() == oh.keys())
    for key in od:
        near(f"T5: same observable probability {key}", od[key], oh[key])
    near("T5: direct CMI zero", rd["cmi_nats"], 0.0)
    near("T5: hidden CMI endpoint", rh["cmi_nats"], hp)
    sd, sh = observer_summary(od), observer_summary(oh)
    require("T5: observer refuses unique answer", sd["status"] == sh["status"] == "not_identifiable")
    for a, b in zip(sd["exclusion_log_loss_interval_nats"], sh["exclusion_log_loss_interval_nats"]):
        near("T5: same observable-only interval", a, b)
    near("T5: interval lower endpoint", sd["exclusion_log_loss_interval_nats"][0], 0.0)
    near("T5: interval upper endpoint", sd["exclusion_log_loss_interval_nats"][1], hp)
    deterministic = observer_summary({(0, 0): 0.5, (1, 1): 0.5})
    require("observer identifies deterministic zero", deterministic["status"] == "identified_zero")
    results["T5_observational_equivalence"] = {
        "direct_scores": rd, "hidden_scores": rh, "observer": sd,
        "observable_rows": [{"x": x, "y": y, "p": mass} for (x, y), mass in sorted(od.items())],
    }

    targets = (-1.0, 1.0)
    residuals = [y - (y + 10.0) for y in targets]
    mean_residual = sum(residuals) / len(residuals)
    residual_variance = sum((r - mean_residual) ** 2 for r in residuals) / len(residuals)
    mse = sum(r * r for r in residuals) / len(residuals)
    target_variance = sum(y * y for y in targets) / len(targets)
    e_style = 1 - residual_variance / target_variance
    near("T6: zero centered residual variance", residual_variance, 0.0)
    near("T6: perfect variance index", e_style, 1.0)
    near("T6: MSE remains 100", mse, 100.0)
    results["T6_bias_blind_variance"] = {
        "residual_mean": mean_residual, "residual_variance": residual_variance,
        "e_style_index": e_style, "mse": mse,
    }

    return {
        "schema": "mrh-sa3a-exact-controls-v1", "status": "pass",
        "scope": "finite population controls; no physical validation or fitted estimator",
        "registration_commit": "5e82a9f4", "absolute_tolerance": TOL,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks_passed": len(checks), "checks": checks, "cases": results,
    }


if __name__ == "__main__":
    print(json.dumps(suite(), indent=2, sort_keys=True, allow_nan=False))
