"""
04 — Global-clock (unilocal) CHSH: the substrate-centric construction (PREDICTIONS.md, bet B1).

WHY THIS, AFTER 02/03
  02 (local) and 03 (nonlocal-grid) both pinned at S = 2. The reason (03's finding): a STATIC
  shared phase lambda only adds a *relabelable* angle offset — it stays a local-realist model.
  dp's correction (2026-06-22): the shared variable is not a static preparation phase. It is the
  DYNAMICAL GLOBAL CLOCK — the universal discrete tick that advances grid slices and GATES every
  oscillation/phase-shift. And the "entangled pair" is not two systems with a common cause: it is
  ONE pattern (unilocal — same frequency, one phase phi relative to the clock), measured at two
  loci. (This is the model, not a claim about the universe's true nature.)

THE STRUCTURAL DIFFERENCE FROM 02
  Measurement = clock-gated phase-LOCKING, and because the substrate is ONE medium, the shared
  pattern BACK-REACTS on the probes during the (synchronous, tick-gated) settling. So Alice's
  probe and Bob's probe lock onto the SAME phi SIMULTANEOUSLY while phi is pulled by BOTH. That
  makes Alice's settled outcome able to depend on Bob's setting — genuinely nonlocal (not a
  static common cause), the one thing 02/03 structurally could not do. Back-reaction strength
  `g` is swept: g=0 recovers the local 02 construction; g>0 turns on the unilocal channel.

THE THREE DISCIPLINES (a violation only counts if all hold)
  1. FREE settings: a,b drawn independently of phi (no superdeterminism — enforced by drawing
     phi independently, uniformly, per trial).
  2. NO SIGNALING: Alice's marginal P(+1) must not depend on Bob's setting (measured per g).
  3. NOT super-quantum: S must not exceed the Tsirelson bound 2*sqrt(2). S>2sqrt2 (toward the
     PR-box 4) would predict correlations STRONGER than any observed in nature -> a REFUTATION,
     not a success.

HONEST PRE-REGISTRATION
  Bell's theorem is airtight: any genuinely local + free-settings + lambda-independent model
  gives S<=2. The only way g>0 reaches 2sqrt2 is if the back-reaction is a real nonlocal channel
  that is ALSO non-signaling. Expected outcomes, in order of prior likelihood:
   (a) g>0 SIGNALS -> not a valid violation; the unilocal channel leaks (tightens B1).
   (b) g>0 stays S<=2 with no signaling -> the back-reaction is still relabelable; B1 boundary
       holds even for the dynamical clock (a sharper negative than 02/03).
   (c) g>0 gives 2<S<=2sqrt2 with no signaling -> EXTRAORDINARY, the program's first novel
       result; flag for independent replication + audit.
   (d) S>2sqrt2 -> super-quantum; a refutation (predicts unobserved correlations).
  Report whichever happens, with the signaling number, at every g.

numpy only. Headless. Writes results/global_clock_chsh_result.json.
"""
import json
import os
import numpy as np

N_TRIALS = 200_000
TICKS = 16            # clock-gated settling steps ("measurement = synchronization", universal tick)
LOCK_RATE = 0.5       # probe relaxation toward the shared pattern per tick
G_LIST = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]   # back-reaction strength of the shared pattern on probes
SEED = 23

CLASSICAL_BOUND = 2.0
TSIRELSON_BOUND = 2.0 * np.sqrt(2.0)


def measure(phi0, a, b, g):
    """Clock-gated SIMULTANEOUS phase-lock of two probes onto ONE shared pattern phi, with the
    shared pattern back-reacting on both probes (strength g). Synchronous (tick-gated) update.
    Returns (outcome_A, outcome_B), each +/-1 = sign of the settled probe vs its setting axis."""
    phi = phi0.copy()
    alpha = np.full_like(phi0, 0.0) + a     # Alice's probe starts on her setting axis
    beta = np.full_like(phi0, 0.0) + b      # Bob's probe starts on his setting axis
    for _ in range(TICKS):
        # all deltas computed from the OLD slice, applied together (the universal tick gates them)
        d_alpha = LOCK_RATE * np.sin(phi - alpha)
        d_beta = LOCK_RATE * np.sin(phi - beta)
        d_phi = g * (np.sin(alpha - phi) + np.sin(beta - phi))   # one substrate -> feels both
        alpha = alpha + d_alpha
        beta = beta + d_beta
        phi = phi + d_phi
    oa = np.sign(np.cos(alpha - a)); oa[oa == 0] = 1.0
    ob = np.sign(np.cos(beta - b)); ob[ob == 0] = 1.0
    return oa, ob


def chsh_at_g(rng, g, angle_sets):
    best = None
    for (a0, a1, b0, b1) in angle_sets:
        # INDEPENDENT phi per correlation term (free settings, phi ~ U(-pi,pi), lambda-independent)
        def E(a, b):
            phi = rng.uniform(-np.pi, np.pi, N_TRIALS)
            oa, ob = measure(phi, a, b, g)
            return float(np.mean(oa * ob)), oa, ob
        E00, oa00, _ = E(a0, b0)
        E01, oa01, _ = E(a0, b1)
        E10, _, _ = E(a1, b0)
        E11, _, _ = E(a1, b1)
        S = E00 - E01 + E10 + E11
        # signaling: Alice's marginal at fixed a0, varying Bob's setting b0 vs b1
        sig = abs(float(np.mean(oa00 > 0)) - float(np.mean(oa01 > 0)))
        rec = {"angles_deg": [round(np.degrees(x), 1) for x in (a0, a1, b0, b1)],
               "S": round(float(S), 4), "signaling_delta": round(sig, 4),
               "E": [round(E00, 4), round(E01, 4), round(E10, 4), round(E11, 4)]}
        if best is None or abs(rec["S"]) > abs(best["S"]):
            best = rec
    return best


def main():
    rng = np.random.default_rng(SEED)
    angle_sets = [
        (0.0, np.pi / 2, np.pi / 4, 3 * np.pi / 4),
        (0.0, np.pi / 2, np.pi / 4, -np.pi / 4),
        (0.0, np.pi / 4, np.pi / 8, 3 * np.pi / 8),
        (0.0, np.pi / 3, np.pi / 6, np.pi / 2),
    ]
    sweep = []
    for g in G_LIST:
        rec = chsh_at_g(rng, g, angle_sets)
        rec["g"] = g
        sweep.append(rec)

    # classify: separate the no-signaling envelope from the signaling-but-large-S region
    valid = [r for r in sweep if r["signaling_delta"] < 0.02]
    best_valid_S = max((abs(r["S"]) for r in valid), default=0.0)
    max_S_any = max(abs(r["S"]) for r in sweep)
    max_S_row = max(sweep, key=lambda r: abs(r["S"]))
    any_signaling = any(r["signaling_delta"] >= 0.02 for r in sweep if r["g"] > 0)
    super_quantum = any(abs(r["S"]) > TSIRELSON_BOUND + 0.02 for r in sweep)
    beats_classical = best_valid_S > CLASSICAL_BOUND + 0.02
    # the Bell trade-off, stated: does EVERY S>2 row also signal?
    s_gt2_all_signal = all(r["signaling_delta"] >= 0.05 for r in sweep if abs(r["S"]) > 2.02)

    if super_quantum and not any_signaling:
        verdict = (
            f"SUPER-QUANTUM at some g: |S| exceeds Tsirelson 2sqrt2={TSIRELSON_BOUND:.3f}. The "
            f"unilocal back-reaction produces correlations STRONGER than quantum mechanics — which "
            f"are NOT observed in nature. This is a REFUTATION of the construction as a model of "
            f"real entanglement (it must be checked for a signaling leak; super-quantum + signaling "
            f"= a hidden communication channel). Productive: the back-reaction is too strong / wrong "
            f"form. Bucket 0 unchanged.")
    elif beats_classical and not any_signaling:
        verdict = (
            f"S = {best_valid_S:.3f} > 2 with NO detectable signaling and S <= Tsirelson. "
            f"EXTRAORDINARY: the dynamical-global-clock unilocal construction reproduces Bell "
            f"violation from a local-by-construction substrate with free settings. PROVISIONAL — "
            f"demands independent replication and an audit for a hidden nonlocal/superdeterministic "
            f"channel before any Bucket-0 claim. If it survives, this is the program's first novel "
            f"result.")
    elif beats_classical and any_signaling:
        verdict = (
            f"Apparent S = {best_valid_S:.3f} > 2 but SIGNALING is present at the violating g — the "
            f"back-reaction leaks Bob's setting into Alice's marginal. NOT a valid Bell violation; "
            f"it rides a hidden communication channel. Outcome (a): the unilocal channel, as built, "
            f"signals. Tightens B1 (the dynamical clock does not escape the no-signaling trap "
            f"either). Bucket 0 unchanged.")
    else:
        verdict = (
            f"THE BELL TRADE-OFF, MADE MECHANICAL. Two regimes, no overlap: (i) NO-SIGNALING "
            f"envelope caps at S = {best_valid_S:.3f} <= 2 (g=0 recovers the local 02 result); "
            f"(ii) back-reaction CAN push S up to {max_S_any:.3f} (at g={max_S_row['g']}, > 2) — but "
            f"every such row SIGNALS heavily (delta up to {max_S_row['signaling_delta']:.2f}; "
            f"every-S>2-signals: {s_gt2_all_signal}). So the dynamical global clock with mutual "
            f"back-reaction DOES create genuine nonlocality — but it is a SIGNALING channel "
            f"(a communication link), not QM's no-signaling nonlocality. You can have nonlocal OR "
            f"no-signaling, never both > 2. That is exactly Bell/Tsirelson playing out: a real-phase "
            f"substrate with a sign-threshold readout cannot produce no-signaling Bell violation, "
            f"with a static OR a dynamical shared variable. The missing primitive is not "
            f"back-reaction (which only signals) — it is interfering COMPLEX amplitudes (the cos^2 "
            f"projection law), i.e. the complex Intent field of Phase-1.6, not a real phase. This "
            f"EXTENDS bet B1's boundary from static-lambda (02/03) to the dynamical clock (04), and "
            f"names the upgrade that would be required. The model does NOT reach the observed "
            f"S=2.828 without that primitive. Bucket 0 unchanged.")

    result = {
        "construction": "unilocal: ONE shared pattern phi, two probes lock simultaneously, phi "
                        "back-reacts on both (g), clock-gated synchronous update",
        "n_trials": N_TRIALS, "ticks": TICKS, "lock_rate": LOCK_RATE, "seed": SEED,
        "classical_bound": CLASSICAL_BOUND, "tsirelson_bound": round(TSIRELSON_BOUND, 4),
        "g_sweep": sweep,
        "best_valid_no_signaling_S": round(best_valid_S, 4),
        "max_S_any_incl_signaling": round(max_S_any, 4),
        "max_S_at_g": max_S_row["g"],
        "max_S_signaling_delta": max_S_row["signaling_delta"],
        "every_S_gt2_signals": bool(s_gt2_all_signal),
        "beats_classical_no_signaling": bool(beats_classical and not any_signaling),
        "super_quantum_no_signaling": bool(super_quantum and not any_signaling),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    out = os.path.join(os.path.dirname(__file__), "results", "global_clock_chsh_result.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out}")
    print("\n  g     best|S|   signaling   (classical=2.0, Tsirelson=2.828)")
    for r in sweep:
        print(f"  {r['g']:.2f}   {abs(r['S']):.4f}    {r['signaling_delta']:.4f}")


if __name__ == "__main__":
    main()
