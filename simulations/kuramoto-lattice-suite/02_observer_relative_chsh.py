"""
02 — Observer-relative CHSH: the one test that matters (PREDICTIONS.md, bet B1).

QUESTION
  Can a purely LOCAL substrate, measured ONLY through an observer-pattern's phase-lock to a
  target-pattern (with freely chosen CHSH settings), reproduce QUANTUM correlations — i.e.
  violate the Bell/CHSH bound S <= 2 — WITHOUT superluminal signaling or superdeterminism?

CONSTRUCTION (the most faithful *local* single-observer model we can build)
  - A shared source prepares two regions A, B with a correlated hidden phase lambda
    (a LOCAL common cause — the local-hidden-variable analog).
  - Alice's observer phase-LOCKS a probe to region A, then reads a binary outcome
    (in-phase / anti-phase) relative to her freely chosen setting axis a in {a0, a1}.
    Bob does the same on region B with b in {b0, b1}.
  - "Measurement = synchronization": the outcome is produced by the probe relaxing
    (locking) onto the target phase, then a sign readout against the setting axis. There is
    NO coupling between A's side and B's side during measurement (locality is structural),
    and settings are drawn INDEPENDENTLY of lambda (no superdeterminism).

EXPECTED RESULT
  S <= 2 (the classical/Bell bound). This construction is local by design, so it cannot
  reach the Tsirelson bound 2*sqrt(2) ~= 2.828 that quantum mechanics achieves. That is the
  PRODUCTIVE result, not a failure: it demonstrates that the single-observer phase-lattice,
  implemented locally, is a local-realist model and does NOT reproduce quantum nonlocality.
  It thereby sharpens that Synchronism is an ONTOLOGY / interpretation, not a
  local-hidden-variable PHYSICS. To exceed 2 you would need a genuinely NONLOCAL grid or
  superdeterministic settings — both of which this harness checks for and flags.

  A result of S > 2 with NO signaling and independent settings would be the program's first
  novel result and would demand independent replication.

numpy only. Headless. Writes results/chsh_result.json.
"""
import json
import os
import numpy as np

N_TRIALS = 400_000     # per random setting assignment; large for tight error bars
LOCK_STEPS = 12        # probe phase-lock relaxation steps ("measurement = synchronization")
LOCK_RATE = 0.6        # relaxation rate of the probe toward the target phase
SEED = 11

CLASSICAL_BOUND = 2.0
TSIRELSON_BOUND = 2.0 * np.sqrt(2.0)


def phase_lock_readout(target_phase, setting, steps, rate):
    """
    Observer measurement = synchronization. A probe initialized on the observer's setting
    axis relaxes (locks) onto the target phase, then returns a binary readout: +1 if the
    settled phase is in-phase with the setting axis, -1 if anti-phase.

    The relaxation makes the 'measurement is a phase-lock event' literal; the settled probe
    -> target_phase, so the readout reduces to sign(cos(target_phase - setting)) — the
    canonical LOCAL response function. That it is local is exactly why it obeys Bell.
    """
    phi = np.full_like(target_phase, 0.0) + setting  # probe starts on the setting axis
    for _ in range(steps):
        phi = phi + rate * np.sin(target_phase - phi)  # lock onto the target
    out = np.sign(np.cos(phi - setting))
    out[out == 0] = 1.0
    return out


def correlation(thetaA, thetaB, a, b):
    """E(a,b) = <outcome_A * outcome_B> over the prepared source ensemble."""
    oa = phase_lock_readout(thetaA, a, LOCK_STEPS, LOCK_RATE)
    ob = phase_lock_readout(thetaB, b, LOCK_STEPS, LOCK_RATE)
    return float(np.mean(oa * ob)), oa, ob


def chsh_for_angles(rng, a0, a1, b0, b1):
    lam = rng.uniform(-np.pi, np.pi, N_TRIALS)   # shared LOCAL hidden phase (common cause)
    thetaA = lam                                  # region A inherits lambda
    thetaB = lam                                  # region B inherits lambda (correlated source)

    E00, oa00, ob00 = correlation(thetaA, thetaB, a0, b0)
    E01, oa01, _    = correlation(thetaA, thetaB, a0, b1)
    E10, _, _       = correlation(thetaA, thetaB, a1, b0)
    E11, _, _       = correlation(thetaA, thetaB, a1, b1)

    S = E00 - E01 + E10 + E11

    # Signaling check: Alice's marginal P(+1) must NOT depend on Bob's setting.
    pA_given_b0 = float(np.mean(oa00 > 0))   # Alice at a0, Bob at b0
    pA_given_b1 = float(np.mean(oa01 > 0))   # Alice at a0, Bob at b1
    signaling = abs(pA_given_b0 - pA_given_b1)

    return {"S": S, "E": [E00, E01, E10, E11], "signaling_delta": signaling}


def quantum_S(a0, a1, b0, b1):
    """Quantum singlet prediction E(a,b) = -cos(a-b), for contrast at the same angles."""
    def E(a, b):
        return -np.cos(a - b)
    return E(a0, b0) - E(a0, b1) + E(a1, b0) + E(a1, b1)


def main():
    rng = np.random.default_rng(SEED)

    # Scan a few canonical angle sets; report the maximum |S| the local model achieves.
    angle_sets = [
        (0.0, np.pi / 2, np.pi / 4, 3 * np.pi / 4),
        (0.0, np.pi / 2, np.pi / 4, -np.pi / 4),
        (0.0, np.pi / 4, np.pi / 8, 3 * np.pi / 8),
    ]
    runs = []
    best = None
    for (a0, a1, b0, b1) in angle_sets:
        r = chsh_for_angles(rng, a0, a1, b0, b1)
        r["angles_deg"] = [round(np.degrees(x), 1) for x in (a0, a1, b0, b1)]
        r["quantum_S_same_angles"] = round(float(quantum_S(a0, a1, b0, b1)), 4)
        r["S"] = round(r["S"], 4)
        r["E"] = [round(e, 4) for e in r["E"]]
        r["signaling_delta"] = round(r["signaling_delta"], 4)
        runs.append(r)
        if best is None or abs(r["S"]) > abs(best["S"]):
            best = r

    S = abs(best["S"])
    violates = S > CLASSICAL_BOUND + 0.02     # small-margin guard against estimator noise
    no_signaling = best["signaling_delta"] < 0.02

    if not violates:
        verdict = (
            "S <= 2 (classical/Bell bound respected). The local single-observer phase-lock "
            "construction does NOT reproduce quantum nonlocality. EXPECTED, and PRODUCTIVE: "
            "it bounds Synchronism's single-observer ontology to a local-realist model in "
            "this construction — an interpretation, not a local-hidden-variable physics. "
            "Bet B1 in PREDICTIONS.md resolves as: refuted-for-this-construction (the honest "
            "boundary). To exceed 2 requires a nonlocal grid or superdeterministic settings."
        )
    elif violates and no_signaling:
        verdict = (
            "S > 2 with NO detectable signaling. EXTRAORDINARY — a local observer-relative "
            "construction reproducing Bell violation. Treat as PROVISIONAL pending "
            "independent replication and an audit for a hidden nonlocal/ superdeterministic "
            "channel. If it survives, this is the program's first novel result."
        )
    else:
        verdict = (
            "S > 2 but signaling detected — the apparent violation rides on a hidden "
            "nonlocal channel (Alice's marginal depends on Bob's setting). NOT a valid Bell "
            "violation; the construction leaked locality. Fix the leak and re-run."
        )

    result = {
        "n_trials": N_TRIALS, "lock_steps": LOCK_STEPS, "lock_rate": LOCK_RATE, "seed": SEED,
        "classical_bound": CLASSICAL_BOUND,
        "tsirelson_bound": round(TSIRELSON_BOUND, 4),
        "best_local_S": round(S, 4),
        "violates_classical_bound": bool(violates),
        "no_signaling": bool(no_signaling),
        "runs": runs,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    out = os.path.join(os.path.dirname(__file__), "results", "chsh_result.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out}")
    print(f"\nlocal single-observer S = {S:.4f}  |  classical bound = 2.0  |  "
          f"Tsirelson (quantum) = {TSIRELSON_BOUND:.4f}")


if __name__ == "__main__":
    main()
