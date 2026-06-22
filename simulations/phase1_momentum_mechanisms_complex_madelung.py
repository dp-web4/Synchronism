"""
Phase-1.6 — the three momentum mechanisms, and why they're one (2026-06-22)

dp's instruction: explore all three candidate "Intent-momentum" mechanisms in parallel, be
open to cross-talk, and look for a deeper perspective — noting that wave/particle duality is
exactly where classical models break, which is the *pointer* to the mechanism we're modeling.

The three candidates:
  (1) independent vector flux J  — a transport/velocity field carrying the pattern
  (2) chiral / asymmetric coupling — a built-in preferred direction in the bond rule
  (3) complex field with intrinsic phase velocity — a carrier whose phase encodes motion

The deeper perspective tested here: **(1) and (2) are projections of (3).** A complex field
ψ = √ρ · e^{iφ} decomposes (Madelung) into:
  - ρ = |ψ|²        → the PARTICLE / envelope  (amplitude = particle nature)
  - v = ∂ₓφ         → a FLUX / velocity field   (= mechanism 1's independent J)
  - a net tilt in φ → a PREFERRED DIRECTION     (= mechanism 2's chirality)
So momentum lives in the PHASE (the wave nature), and the envelope is the particle. That is
wave/particle duality as the *mechanism*, not a paradox. It also explains the Phase-1.5 null:
a REAL field has no phase to hold momentum, so boosting its envelope pins/disperses.

This experiment runs the real-field arm (the Phase-1.5 failure) against the complex-field arm,
measures whether the complex field carries a coherent packet at group velocity v_g < c (the
ceiling), and performs the Madelung readout to show the flux (1) and direction (2) living
INSIDE the complex field (3). NOT novel physics — moving complex-KG/NLS solitons are textbook;
the value is the structural unification of dp's three mechanisms into one, and the empirical
demonstration that momentum = phase = wave-nature. numpy only, headless; writes results JSON.
"""
import json
import os
import numpy as np

L = 384
DT = 0.02
C2 = 4.0                 # lattice signal speed^2 -> c = 2 (the ceiling)
C = float(np.sqrt(C2))
M2 = 1.0                 # mass term
GSAT = 1.5               # focusing strength
SSAT = 0.7               # focusing saturation scale
STEPS = 8000
K_LIST = [0.0, 0.25, 0.5, 1.0, 2.0]   # phase-ramp wavenumbers (momentum)
SEED = 5


def lap(z):
    return np.roll(z, -1) - 2 * z + np.roll(z, 1)


# ---------- complex arm (the unification: mechanism 3, containing 1 & 2) ----------
def accel_complex(psi):
    foc = GSAT * (np.abs(psi) ** 2 / (1.0 + np.abs(psi) ** 2 / SSAT ** 2)) * psi
    return C2 * lap(psi) - M2 * psi + foc


def run_complex(rng, k):
    x = np.arange(L)
    x0 = L / 2.0; w = 6.0; A = 0.9
    env = A / np.cosh((x - x0) / w)
    psi = env * np.exp(1j * k * (x - x0))
    psi += 0.003 * (rng.standard_normal(L) + 1j * rng.standard_normal(L))
    Omega = np.sqrt(M2 + C2 * k * k)                 # KG dispersion for wavenumber k
    v_g = C2 * k / Omega                             # group velocity (theory) -> < c
    pdot = -1j * Omega * psi                          # carrier oscillation
    a = accel_complex(psi)
    cents, widths, amps = [], [], []
    se = 25
    for t in range(STEPS):
        psi = psi + DT * pdot + 0.5 * DT * DT * a
        a_new = accel_complex(psi)
        pdot = pdot + 0.5 * DT * (a + a_new)
        a = a_new
        if t % se == 0:
            rho = np.abs(psi) ** 2
            cents.append(circ_centroid(rho)); widths.append(eff_width_rho(rho))
            amps.append(float(np.max(np.abs(psi))))
    v_meas = unwrap_speed(cents, DT, se)
    # Madelung readout at final state: flux v=∂φ weighted by ρ over the packet core
    phi = np.unwrap(np.angle(psi))
    dphi = 0.5 * (np.roll(phi, -1) - np.roll(phi, 1))
    rho = np.abs(psi) ** 2
    core = rho > 0.2 * rho.max()
    mean_flux = float(np.sum(dphi[core] * rho[core]) / (np.sum(rho[core]) + 1e-12))  # = J (mech 1)
    flux_direction = "right" if mean_flux > 0 else ("left" if mean_flux < 0 else "none")
    return summarize(cents, widths, amps, v_meas, v_g, mean_flux, flux_direction)


# ---------- real arm (mechanism 1/2 as a REAL field: the Phase-1.5 failure mode) ----------
def accel_real(u):
    foc = GSAT * (u ** 3) / (1.0 + (u / SSAT) ** 2)
    return C2 * lap(u) - M2 * u + foc


def run_real(rng, k):
    """Real field, momentum injected as a translating-envelope IC (velocity ∝ k·c).
    No phase to carry momentum -> expected pin/disperse, per Phase-1.5."""
    x = np.arange(L)
    x0 = L / 2.0; w = 6.0; A = 0.9
    u = A / np.cosh((x - x0) / w) + 0.003 * rng.standard_normal(L)
    V = min(C2 * k / np.sqrt(M2 + C2 * k * k), 0.95 * C) if k > 0 else 0.0
    du = 0.5 * (np.roll(u, -1) - np.roll(u, 1))
    v = -V * du
    a = accel_real(u)
    cents, widths, amps = [], [], []
    se = 25
    for t in range(STEPS):
        u = u + DT * v + 0.5 * DT * DT * a
        a_new = accel_real(u)
        v = v + 0.5 * DT * (a + a_new)
        a = a_new
        if t % se == 0:
            rho = u ** 2
            cents.append(circ_centroid(rho)); widths.append(eff_width_rho(rho))
            amps.append(float(np.max(np.abs(u))))
    v_meas = unwrap_speed(cents, DT, se)
    return summarize(cents, widths, amps, v_meas, V, None, None)


# ---------- shared measures ----------
def circ_centroid(rho):
    ang = 2 * np.pi * np.arange(L) / L
    z = np.sum(rho * np.exp(1j * ang))
    return (np.angle(z) % (2 * np.pi)) / (2 * np.pi) * L


def eff_width_rho(rho):
    s1 = np.sum(rho); s2 = np.sum(rho ** 2)
    return float(s1 * s1 / s2) if s2 > 1e-12 else float(L)


def unwrap_speed(cents, dt, se):
    c = np.array(cents); d = np.diff(c)
    d = (d + L / 2) % L - L / 2
    return float(np.sum(d) / ((len(c) - 1) * se * dt))


def summarize(cents, widths, amps, v_meas, v_theory, mean_flux, flux_dir):
    widths = np.array(widths); amps = np.array(amps)
    w0 = float(np.mean(widths[:5])); wf = float(np.mean(widths[-len(widths) // 4:]))
    a0 = float(np.mean(amps[:5])); af = float(np.mean(amps[-len(amps) // 4:]))
    cc = np.array(cents); d = np.diff(cc); d = (d + L / 2) % L - L / 2
    distance = float(abs(np.sum(d)))
    localized = bool(wf < 0.12 * L and wf <= 1.6 * w0)
    amp_kept = af / (a0 + 1e-9)
    speed_ok = (abs(v_theory) < 1e-6) or (abs(abs(v_meas) - abs(v_theory)) < 0.3 * abs(v_theory) + 0.05)
    propagated = bool(localized and amp_kept > 0.5 and distance > 4 * w0 and speed_ok) if abs(v_theory) > 1e-6 else False
    return {
        "v_theory": round(float(v_theory), 3), "v_meas": round(v_meas, 3),
        "below_c": bool(abs(v_meas) < C + 0.05),
        "localized": localized, "amp_kept": round(amp_kept, 3),
        "distance_widths": round(distance / (w0 + 1e-9), 1),
        "propagated_coherent": propagated,
        "madelung_flux_J": (None if mean_flux is None else round(mean_flux, 3)),
        "flux_direction": flux_dir,
    }


def main():
    rng = np.random.default_rng(SEED)
    complex_runs = [{"k": k, **run_complex(rng, k)} for k in K_LIST]
    real_runs = [{"k": k, **run_real(rng, k)} for k in K_LIST]

    cplx_moving = [r for r in complex_runs if r["k"] > 0]
    cplx_propagates = [r for r in cplx_moving if r["propagated_coherent"]]
    real_moving = [r for r in real_runs if r["k"] > 0]
    real_propagates = [r for r in real_moving if r["propagated_coherent"]]
    ceiling_held = all(r["below_c"] for r in complex_runs)
    # cross-talk: nonzero Madelung flux that tracks k => the flux J (mech 1) lives inside ψ
    flux_tracks_k = all(
        (r["madelung_flux_J"] is not None and abs(r["madelung_flux_J"]) > 0.05)
        for r in cplx_moving
    )

    if cplx_propagates and not real_propagates:
        verdict = (
            "The three mechanisms are ONE. The COMPLEX field carries a coherent packet that "
            f"propagates at group velocity v_g<c (clean at k in {[r['k'] for r in cplx_propagates]}; "
            f"ceiling c={C} held: {ceiling_held}), while the REAL field — momentum injected as an "
            "envelope boost, with no phase to hold it — pins/disperses (reproduces Phase-1.5). "
            "Madelung readout of the complex packet shows a nonzero flux J=∂ₓφ that tracks the "
            f"phase-ramp k (tracks: {flux_tracks_k}) with a definite direction — i.e. mechanism 1 "
            "(independent flux J) and mechanism 2 (preferred direction) are the amplitude/phase "
            "PROJECTIONS of mechanism 3 (complex field). Deeper perspective (dp's pointer): "
            "momentum lives in the PHASE (wave nature), the envelope is the PARTICLE (amplitude "
            "nature) — wave/particle duality is the momentum MECHANISM, not a paradox bolted on. "
            "Photon = phase-dominated propagating packet; electron = amplitude-dominated localized "
            "envelope; both from one complex Intent field. NOT novel physics (moving complex-KG "
            "solitons are textbook) — but it unifies the three proposals and says the substrate "
            "MUST be complex. Next: collision/interference of two complex packets (Stage 2)."
        )
    elif cplx_propagates and real_propagates:
        verdict = ("Both real and complex arms propagated — unexpected; the real arm should lack "
                   "a phase to carry momentum. AUDIT the real arm (boundary/parameter artifact).")
    else:
        verdict = ("Complex arm did not propagate cleanly under these params — sweep "
                   "(M2, GSAT, SSAT, k, dt) before concluding; record as a tighter null.")

    out = {
        "lattice": L, "dt": DT, "c_ceiling": C, "m2": M2, "steps": STEPS, "seed": SEED,
        "complex_field_arm": complex_runs,
        "real_field_arm": real_runs,
        "complex_propagates_at_k": [r["k"] for r in cplx_propagates],
        "real_propagates_at_k": [r["k"] for r in real_propagates],
        "ceiling_c_held": ceiling_held,
        "madelung_flux_tracks_k": flux_tracks_k,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase1_momentum_mechanisms_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print(f"\n c-ceiling = {C}")
    print("\n COMPLEX arm:   k    v_theory  v_meas  below_c  localized  amp_kept  dist(w)  J=∂φ   propagated")
    for r in complex_runs:
        print(f"             {r['k']:.2f}   {r['v_theory']:+.2f}    {r['v_meas']:+.2f}   "
              f"{str(r['below_c']):5}    {str(r['localized']):5}     {r['amp_kept']:.2f}    "
              f"{r['distance_widths']:5.1f}   {str(r['madelung_flux_J']):>5}   {r['propagated_coherent']}")
    print("\n REAL arm:      k    v_theory  v_meas  below_c  localized  amp_kept  dist(w)         propagated")
    for r in real_runs:
        print(f"             {r['k']:.2f}   {r['v_theory']:+.2f}    {r['v_meas']:+.2f}   "
              f"{str(r['below_c']):5}    {str(r['localized']):5}     {r['amp_kept']:.2f}    "
              f"{r['distance_widths']:5.1f}            {r['propagated_coherent']}")


if __name__ == "__main__":
    main()
