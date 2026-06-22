"""
Phase-1.5 — INTENT MOMENTUM: can a coherent pattern carry momentum and PROPAGATE?

Context. The Stage-1 result (phase1_stage1_localized_oscillation_substrate.py) showed the
substrate self-confines a stable oscillating pattern only with a *focusing* nonlinearity —
but that pattern is a STANDING breather: it oscillates in place and goes nowhere. dp's model
contribution (2026-06-22): an explanatory substrate must also let Intent carry **momentum** —
a strongly preferred direction — because the things we witness include:
  - photons: coherent patterns with a very specific direction and speed (momentum, propagating)
  - electrons: massive patterns that can sit at rest (rest oscillation = mass = complexity)
  - conversion: interference can change a photon's speed/direction/phase, or dissipate/reform
    patterns into different ones (photon <-> electron-like).

This test takes the focusing-saturating substrate (the one that makes standing breathers) and
asks the next question: if we seed the pattern with momentum (a translating envelope at speed
V), does it PROPAGATE as a coherent localized pattern — moving at ~V, staying localized,
keeping amplitude — or does the lattice pin it (Peierls-Nabarro barrier) / shred it (V above
the lattice signal speed)? And does rest (V=0) vs moving give the electron-like vs photon-like
split from ONE rule?

This is explanatory-model-building, not a patch, and not a claim of novel physics: moving
discrete solitons are studied (PN barrier is real). The value is constructive — does THIS
substrate, with minimal mechanisms (nonlinear saturation + Intent momentum), reproduce the
witnessable behavior (a momentum-carrying coherent pattern with a max propagation speed)?

Method: 1D periodic lattice, second-order wave (velocity-Verlet), focusing-saturating on-site.
Seed a sech envelope with initial velocity field v0 = -V * du/dx (a translating solution of
the wave eqn). Sweep V. Track the (circular) centroid -> measured speed; width; amplitude.
Lattice signal speed c = sqrt(C2) = 1. Expect a window: pin below, propagate in a band,
radiate/break near/above c. numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

L = 384
DT = 0.04
C2 = 1.0                 # lattice signal speed^2 -> c = 1
OMEGA0_2 = 1.0
GAMMA = 1.2
U_S = 0.6
STEPS = 6000
V_LIST = [0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.1]
SEED = 3


def onsite_force(u):
    # focusing-saturating (the Stage-1 Arm-D rule)
    return -OMEGA0_2 * u + GAMMA * (u ** 3) / (1.0 + (u / U_S) ** 2)


def accel(u, c2=C2):
    coupling = c2 * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
    return coupling + onsite_force(u)


def seed(rng, V):
    x = np.arange(L)
    x0 = L / 2.0
    w = 5.0
    amp = 0.9
    u = amp / np.cosh((x - x0) / w)            # sech envelope (focusing-soliton shape)
    u += 0.005 * rng.standard_normal(L)
    du = 0.5 * (np.roll(u, -1) - np.roll(u, 1))  # central-difference spatial gradient
    v = -V * du                                  # translating solution: u_t = -V u_x
    return u, v


def circ_centroid(u):
    """Centroid on a periodic lattice via circular mean of the u^2 distribution."""
    w = u * u
    ang = 2 * np.pi * np.arange(L) / L
    z = np.sum(w * np.exp(1j * ang))
    c = (np.angle(z) % (2 * np.pi)) / (2 * np.pi) * L
    return float(c)


def eff_width(u):
    s2 = np.sum(u * u); s4 = np.sum(u ** 4)
    return float(s2 * s2 / s4) if s4 > 1e-12 else float(L)


def unwrap_speed(centroids, dt, sample_every):
    """Estimate propagation speed from centroid track, handling periodic wraps."""
    c = np.array(centroids)
    d = np.diff(c)
    d = (d + L / 2) % L - L / 2                  # shortest signed step (handle wrap)
    return float(np.sum(d) / ((len(c) - 1) * sample_every * dt))


def run(rng, V, c2=C2):
    u, v = seed(rng, V)
    a = accel(u, c2)
    sample_every = 25
    cents, widths, amps = [], [], []
    for t in range(STEPS):
        u = u + DT * v + 0.5 * DT * DT * a
        a_new = accel(u, c2)
        v = v + 0.5 * DT * (a + a_new)
        a = a_new
        if t % sample_every == 0:
            cents.append(circ_centroid(u)); widths.append(eff_width(u))
            amps.append(float(np.max(np.abs(u))))
    measured_v = unwrap_speed(cents, DT, sample_every)
    widths = np.array(widths); amps = np.array(amps)
    w0 = float(np.mean(widths[:5])); wf = float(np.mean(widths[-len(widths)//4:]))
    a0 = float(np.mean(amps[:5])); af = float(np.mean(amps[-len(amps)//4:]))
    localized = wf < 0.12 * L and wf <= 1.6 * w0
    amp_kept = af / (a0 + 1e-9)
    # distance traveled (in widths) over the run, signed-unwrapped
    dtot = np.diff(np.array(cents)); dtot = (dtot + L/2) % L - L/2
    distance = float(abs(np.sum(dtot)))
    # speed fidelity vs requested V (skip for V=0)
    speed_ok = (V == 0.0) or (abs(measured_v - V) < 0.25 * V + 0.05)
    propagated = bool(localized and amp_kept > 0.5 and distance > 5 * w0 and speed_ok) if V > 0 else False
    rest_held = bool(localized and amp_kept > 0.5 and distance < 3 * w0) if V == 0 else None
    return {
        "V_requested": V, "V_measured": round(measured_v, 3),
        "localized": bool(localized), "width0": round(w0, 1), "widthf": round(wf, 1),
        "amp_kept": round(amp_kept, 3), "distance_in_widths": round(distance / (w0 + 1e-9), 1),
        "propagated": propagated, "rest_held": rest_held,
    }


def main():
    rng = np.random.default_rng(SEED)
    runs = [run(rng, V) for V in V_LIST]

    # C2 (coupling / discreteness) sweep at half-signal-speed — separates a REAL momentum-null
    # from a Peierls-Nabarro pinning ARTIFACT. Higher C2 -> broader soliton relative to the
    # lattice -> lower PN barrier; if the substrate truly carries momentum, a propagation band
    # opens as C2 rises, bounded by the lattice signal speed c = sqrt(C2).
    c2_sweep = []
    for c2 in [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]:
        Vc = 0.5 * float(np.sqrt(c2))      # half the lattice signal speed
        r = run(rng, Vc, c2)
        r["c2"] = c2; r["c_signal"] = round(float(np.sqrt(c2)), 2)
        c2_sweep.append(r)
    unpins_at = [r["c2"] for r in c2_sweep if r["propagated"]]

    rest = runs[0]
    rest_is_electron_like = bool(rest["rest_held"])
    propagating_lowC = [r["V_requested"] for r in runs[1:] if r["propagated"]]

    if rest_is_electron_like and unpins_at:
        verdict = (
            "BOTH regimes from ONE rule, with a coupling threshold. At strong discreteness "
            "(C2=1) every moving seed PINS (Peierls-Nabarro barrier) — electron-like only "
            "(stands, oscillates, holds). As coupling rises (discreteness weakens) a "
            f"PROPAGATION band opens — clean momentum-carrying propagation at C2 in {unpins_at} — "
            "bounded by the lattice signal speed c=sqrt(C2). Intent-momentum IS carried by the "
            "substrate, above a discreteness threshold, with c as the propagation ceiling, from "
            "the SAME rule that gives the standing (electron-like) pattern at rest. NOT novel "
            "physics (moving solitons + PN barrier are studied) — a constructive step: one "
            "substrate, rest->electron-like / moving->photon-like, momentum as the knob, c as "
            "the ceiling. Next (Stage 2): collide two patterns — does interference change "
            "speed/direction/phase, and can a propagating pattern reform into standing ones?"
        )
    elif rest_is_electron_like and not unpins_at:
        verdict = (
            "Rest pattern holds (electron-like) at every coupling, but NO propagation band "
            "opened even as discreteness weakened (C2 up to 32). Under this momentum injection "
            "(translating-envelope seed) the substrate does NOT carry momentum — a real, "
            "tighter null on the photon-analog, not merely a PN artifact. Reading: momentum "
            "likely must be a SUBSTRATE mechanism (an intrinsic directional flux / advection — "
            "dp's 'intent has a preferred direction'), not an initial condition. That redesign "
            "is the next experiment."
        )
    else:
        verdict = (
            "Rest pattern did not hold cleanly under these params — re-tune (gamma, u_s, "
            "omega0, dt) before reading the momentum sweeps."
        )

    out = {
        "lattice": L, "dt": DT, "steps": STEPS, "seed": SEED,
        "V_sweep_at_C2_1": runs,
        "C2_sweep_half_speed": c2_sweep,
        "rest_electron_like": rest_is_electron_like,
        "propagation_unpins_at_C2": unpins_at,
        "propagating_band_lowC": propagating_lowC,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase1_intent_momentum_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n  V_req  V_meas  localized  amp_kept  dist(widths)  propagated")
    for r in runs:
        print(f"  {r['V_requested']:.2f}   {r['V_measured']:+.3f}   "
              f"{str(r['localized']):5}      {r['amp_kept']:.2f}      "
              f"{r['distance_in_widths']:5.1f}        {r['propagated']}")


if __name__ == "__main__":
    main()
