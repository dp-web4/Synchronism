"""
03 — Nonlocal-grid CHSH: the frontier variant (PREDICTIONS.md bet B1, open arm).

The local construction (02) gave S = 1.98 <= 2: a faithful *local* single-observer model is
local-realist and cannot reproduce quantum nonlocality. The only way the single-observer
ontology could still buy a novel result is if **the one shared substrate ("the grid") does
work across the separated regions** — a NONLOCAL channel — WITHOUT producing classical
signaling (Alice's marginal must not depend on Bob's setting).

This harness instruments exactly that. A tunable nonlocal coupling g lets the grid mix
region B's state-and-setting into Alice's measurement (and symmetrically). At g=0 it reduces
to the local model (02). We sweep g and report, at each g, BOTH:
  - S         : the CHSH value (does it exceed the classical bound 2?)
  - signaling : |P(A=+ | b0) - P(A=+ | b1)| (does Alice's marginal leak Bob's setting?)

THE QUESTION: is there ANY g where S > 2 with signaling ~ 0 (a no-signaling nonlocal
violation, i.e. genuinely quantum-like)? Or does the single-grid coupling only ever exceed
the bound by simultaneously producing signaling (an unphysical, invalid 'violation')?

EXPECTED: S > 2 only appears together with signaling > 0. That would be the productive
boundary: the single shared substrate cannot have it both ways — couple weakly and stay
local-realist (S<=2), or couple strongly and signal. A no-signaling violation would require
a quantum-entanglement-like primitive the ontology does NOT yet contain (and would have to
DERIVE, not assume). If instead some g gives S>2 with signaling~0, that is the program's
first genuinely novel result and demands independent replication + a hidden-channel audit.

numpy only. Headless. Writes results/nonlocal_chsh_result.json.
"""
import json
import os
import numpy as np

N_TRIALS = 400_000
SEED = 23
CLASSICAL_BOUND = 2.0
TSIRELSON_BOUND = 2.0 * np.sqrt(2.0)
G_VALUES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]

# Canonical angle sets; we take the max |S| over them at each g.
ANGLE_SETS = [
    (0.0, np.pi / 2, np.pi / 4, 3 * np.pi / 4),
    (0.0, np.pi / 2, np.pi / 4, -np.pi / 4),
    (0.0, np.pi / 4, np.pi / 8, 3 * np.pi / 8),
]


def outcomes(lam, a, b, g):
    """
    Nonlocal-grid measurement. Alice's effective measurement phase is a g-weighted mix of
    her own target (thetaA = lam) and — carried by the shared grid — Bob's target rotated by
    Bob's setting (thetaB - b). Symmetric for Bob. g=0 -> purely local (reduces to sim 02).
    Outcome = sign(cos(effective_phase - own_setting)).
    """
    thetaA = lam
    thetaB = lam
    effA = np.angle((1.0 - g) * np.exp(1j * thetaA) + g * np.exp(1j * (thetaB - b)))
    effB = np.angle((1.0 - g) * np.exp(1j * thetaB) + g * np.exp(1j * (thetaA - a)))
    oa = np.sign(np.cos(effA - a)); oa[oa == 0] = 1.0
    ob = np.sign(np.cos(effB - b)); ob[ob == 0] = 1.0
    return oa, ob


def chsh_at(rng, g, a0, a1, b0, b1):
    lam = rng.uniform(-np.pi, np.pi, N_TRIALS)
    oa00, ob00 = outcomes(lam, a0, b0, g)
    oa01, ob01 = outcomes(lam, a0, b1, g)
    oa10, ob10 = outcomes(lam, a1, b0, g)
    oa11, ob11 = outcomes(lam, a1, b1, g)
    E00 = float(np.mean(oa00 * ob00))
    E01 = float(np.mean(oa01 * ob01))
    E10 = float(np.mean(oa10 * ob10))
    E11 = float(np.mean(oa11 * ob11))
    S = E00 - E01 + E10 + E11
    # signaling: Alice at a0, does her marginal depend on Bob's setting (b0 vs b1)?
    signaling = abs(float(np.mean(oa00 > 0)) - float(np.mean(oa01 > 0)))
    return S, signaling


def main():
    rng = np.random.default_rng(SEED)
    sweep = []
    for g in G_VALUES:
        best_S, sig_at_best = 0.0, 0.0
        for (a0, a1, b0, b1) in ANGLE_SETS:
            S, sig = chsh_at(rng, g, a0, a1, b0, b1)
            if abs(S) > abs(best_S):
                best_S, sig_at_best = S, sig
        sweep.append({"g": g, "S": round(abs(best_S), 4), "signaling": round(sig_at_best, 4)})

    # Is there a no-signaling violation anywhere on the sweep?
    no_sig_violation = [r for r in sweep if r["S"] > CLASSICAL_BOUND + 0.02 and r["signaling"] < 0.02]
    any_violation = [r for r in sweep if r["S"] > CLASSICAL_BOUND + 0.02]

    if no_sig_violation:
        verdict = (
            "FOUND a no-signaling nonlocal violation (S>2, signaling~0) at "
            f"{no_sig_violation}. EXTRAORDINARY and PROVISIONAL — the single-grid channel "
            "would be reproducing quantum-like nonlocality without signaling. Demands "
            "independent replication and an audit for a hidden classical channel before any "
            "claim. If it survives, this is the program's first novel result."
        )
    elif any_violation:
        verdict = (
            "Violation (S>2) occurs ONLY together with signaling > 0. The single shared "
            "substrate cannot have it both ways: weak coupling stays local-realist (S<=2); "
            "strong enough coupling to exceed the bound also makes Alice's marginal leak "
            "Bob's setting (unphysical). PRODUCTIVE BOUNDARY: a no-signaling violation would "
            "require a quantum-entanglement-like primitive the ontology does not contain and "
            "would have to DERIVE, not assume. Bet B1 stays refuted; the gap is now precisely "
            "named."
        )
    else:
        verdict = (
            "No violation at any g on the sweep (S stays <= 2). The nonlocal-grid mixing as "
            "constructed does not exceed the classical bound. Bet B1 stays refuted; consider "
            "richer grid-mediation forms before concluding the channel is exhausted."
        )

    result = {
        "n_trials": N_TRIALS, "seed": SEED,
        "classical_bound": CLASSICAL_BOUND, "tsirelson_bound": round(TSIRELSON_BOUND, 4),
        "sweep": sweep,
        "no_signaling_violation_found": bool(no_sig_violation),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    out = os.path.join(os.path.dirname(__file__), "results", "nonlocal_chsh_result.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out}")
    print("\n g     S      signaling")
    for r in sweep:
        print(f"  {r['g']:.2f}  {r['S']:.3f}   {r['signaling']:.3f}")


if __name__ == "__main__":
    main()
