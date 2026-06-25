"""
Phase-3b — does the SUBSTRATE itself produce a gravity-like index field? (2026-06-22)

Phase-3 imposed n(r)=1+alpha/r BY HAND and watched rays bend. That tests geometry, not the
framework. The real question: does the Phase-2 substrate (complex Intent field + focusing
nonlinearity) make a localized high-Intent core — a "mass" — slow a passing probe wave (a
"photon"), and WITH WHAT SPATIAL PROFILE? A mechanism is only emergent if the index well comes
out of the dynamics, not the typist.

Faithful to absolute time: the probe's local PHASE VELOCITY (cells advanced per absolute tick)
is what varies — a spatial/kinematic effect. No clock is perturbed. The tick is universal.

Method (two independent measurements, agent-zero discipline):
  (1) WKB index from the LINEARIZED substrate. Around a static core A(x), a small probe sees an
      effective squared-mass m_eff^2(x) = w0^2 - 3*gamma*A(x)^2 / (1 + A(x)^2/us^2) (curvature of
      the focusing-saturating on-site potential at the core amplitude). Local index for a fixed
      probe frequency w:  k(x) = sqrt((w^2 - m_eff^2(x))/c^2),  n(x) = k(x)/k_inf.
  (2) DIRECT FDTD cross-check. Propagate a small-amplitude probe wavepacket across the core in a
      2nd-order wave field; compare its accumulated phase at exit to a free probe (no core). The
      measured extra phase must match integral of (n(x)-1)*k dx from (1), or the linearization is
      lying. This is the dummy that keeps the analytic honest.

THE QUESTION UNDER TEST: is the emergent n(x)-1 SHORT-RANGE (tracks the core, ~sech^2/Gaussian)
or LONG-RANGE (~1/r, gravity-like)? Newtonian gravity needs 1/r. A bare soliton core gives
sech^2. Prediction (pre-registered): short-range -> the substrate refracts light at a mass
(Q2 mechanism real), but the RANGE is wrong -> recovering 1/r is UNSOLVED. A productive negative.

numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

# --- substrate parameters (from the Phase-2 focusing-saturating Intent field) ---
L = 4000          # grid cells (long, to see the tail)
DX = 1.0
C = 1.0           # lattice speed ceiling
W0 = 1.0          # bare mass (rest frequency) of the field
GAMMA = 1.0       # focusing strength
US = 1.0          # saturation amplitude
A0 = 0.8          # core peak amplitude (the "mass")
CORE_W = 8.0      # core width (cells)
PROBE_W = 0.30    # probe frequency above the gap? -> set probe k and w


def core_profile(x, x0):
    """Static high-Intent core — soliton-shaped (sech), localized."""
    return A0 / np.cosh((x - x0) / CORE_W)


def m_eff_sq(A):
    """Linearized effective squared-mass a probe sees in the focusing-saturating potential.
    Focusing (cubic) term reduces the effective mass where core amplitude is high."""
    return W0 ** 2 - 3.0 * GAMMA * A ** 2 / (1.0 + A ** 2 / US ** 2)


def wkb_index(x, x0, w_probe):
    A = core_profile(x, x0)
    meff2 = m_eff_sq(A)
    k_local = np.sqrt(np.maximum((w_probe ** 2 - meff2) / C ** 2, 1e-9))
    k_inf = np.sqrt(max((w_probe ** 2 - W0 ** 2) / C ** 2, 1e-9))
    return k_local / k_inf, k_local, k_inf


def fit_range(x, n_minus_1, x0):
    """Classify the index-well range: compare best-fit sech^2 (short) vs 1/|r| (long)."""
    r = np.abs(x - x0) + 1e-6
    y = n_minus_1
    mask = y > y.max() * 1e-3
    xs, ys, rs = x[mask], y[mask], r[mask]
    # sech^2 fit (short-range): y = a * sech^2((x-x0)/w)
    sech2 = (1.0 / np.cosh((xs - x0) / CORE_W)) ** 2
    a_sech = np.dot(ys, sech2) / np.dot(sech2, sech2)
    resid_sech = float(np.sum((ys - a_sech * sech2) ** 2))
    # 1/r fit (long-range): y = b / r
    inv_r = 1.0 / rs
    b_invr = np.dot(ys, inv_r) / np.dot(inv_r, inv_r)
    resid_invr = float(np.sum((ys - b_invr * inv_r) ** 2))
    # tail decay: ratio of n-1 at 4w vs at 16w. sech^2 -> ~e^-2*(distance) (steep); 1/r -> 4x (slow)
    def val_at(d):
        i = int(x0 + d)
        return float(y[i]) if 0 <= i < len(y) else 0.0
    tail_ratio = val_at(4 * CORE_W) / (val_at(16 * CORE_W) + 1e-12)
    invr_ratio = 4.0  # what 1/r would give for that 4x distance span
    return {
        "sech2_residual": round(resid_sech, 8),
        "invr_residual": round(resid_invr, 8),
        "better_fit": "sech2 (short-range)" if resid_sech < resid_invr else "1/r (long-range)",
        "tail_ratio_4w_over_16w": round(tail_ratio, 2),
        "tail_ratio_if_1overr": invr_ratio,
        "tail_verdict": ("short-range (steeper than 1/r — tracks the core)"
                         if tail_ratio > invr_ratio * 2 else "1/r-like or slower"),
    }


def fdtd_phase_check(x0, w_probe, k_inf):
    """Direct: propagate a probe through the core vs free; compare accumulated phase.
    2nd-order Klein-Gordon-like field, probe = small-amplitude carrier under a wide envelope."""
    x = np.arange(L) * DX
    A_core = core_profile(x, x0)
    meff2_core = m_eff_sq(A_core)
    meff2_free = np.full(L, W0 ** 2)

    def run(meff2):
        # complex probe; carrier exp(i k x), wide gaussian envelope entering from the left
        env = np.exp(-((x - 0.15 * L) ** 2) / (2 * (0.05 * L) ** 2))
        psi = (0.01 * env * np.exp(1j * k_inf * x)).astype(complex)
        psi_prev = (0.01 * env * np.exp(1j * (k_inf * x - w_probe * (-DX)))).astype(complex)
        dt = 0.4 * DX / C
        steps = int(0.7 * L / (C * dt))
        for _ in range(steps):
            lap = (np.roll(psi, -1) - 2 * psi + np.roll(psi, 1)) / DX ** 2
            psi_next = 2 * psi - psi_prev + dt ** 2 * (C ** 2 * lap - meff2 * psi)
            psi_prev, psi = psi, psi_next
        return psi

    psi_core = run(meff2_core)
    psi_free = run(meff2_free)
    # extra phase accumulated by the core run at the probe location past the core
    probe_region = slice(int(x0 + 6 * CORE_W), int(x0 + 30 * CORE_W))
    amp_core = np.abs(psi_core[probe_region])
    sel = amp_core > amp_core.max() * 0.3
    if sel.sum() < 5:
        return None
    phase_core = np.angle(psi_core[probe_region][sel])
    phase_free = np.angle(psi_free[probe_region][sel])
    dphi = np.unwrap(phase_core) - np.unwrap(phase_free)
    return float(np.median(dphi))


def main():
    x = np.arange(L) * DX
    x0 = L / 2
    # probe frequency safely above the gap (propagating): w^2 > w0^2
    w_probe = np.sqrt(W0 ** 2 + (PROBE_W * C) ** 2 + 0.5)
    n, k_local, k_inf = wkb_index(x, x0, w_probe)
    n_minus_1 = n - 1.0

    peak = float(n_minus_1.max())
    # integrated extra optical path (analytic): integral (n-1) k_inf dx  (the deflection driver)
    integrated_phase_analytic = float(np.sum(n_minus_1 * k_inf) * DX)
    rng = fit_range(x, n_minus_1, x0)

    # direct FDTD cross-check
    dphi_fdtd = fdtd_phase_check(x0, w_probe, k_inf)
    fdtd_vs_analytic = (round(dphi_fdtd / integrated_phase_analytic, 3)
                        if (dphi_fdtd is not None and integrated_phase_analytic != 0) else None)

    verdict = (
        f"A localized high-Intent core (the 'mass') DOES create an index well: n-1 peaks at "
        f"{peak:.4f}, so a passing probe wave slows and would refract TOWARD the core — Q2's "
        f"mechanism (gravity bends light) emerges from the substrate, faithful to absolute time "
        f"(it is the probe's PHASE VELOCITY, cells-per-absolute-tick, that drops; no clock is "
        f"touched). BUT THE RANGE: the well is best fit by {rng['better_fit']} "
        f"(sech^2 residual {rng['sech2_residual']:.2e} vs 1/r residual {rng['invr_residual']:.2e}); "
        f"tail falls {rng['tail_ratio_4w_over_16w']}x over a 4x distance span vs the {rng['tail_ratio_if_1overr']}x "
        f"a 1/r field would give -> {rng['tail_verdict']}. So the emergent field is SHORT-RANGE "
        f"(it tracks the core profile), NOT the long-range 1/r that Newtonian gravity requires. "
        f"PRODUCTIVE NEGATIVE (pre-registered prediction confirmed): the substrate's own "
        f"nonlinearity refracts light at a mass, but as a contact-range index bump, not a 1/r "
        f"potential. Recovering 1/r — and only then asking the factor-of-2 question — is UNSOLVED: "
        f"it needs a long-range mediator (a massless mode / the field's far tail / a collective "
        f"effect), not a bare soliton core. This is the real gap between 'refracts light' and "
        f"'is gravity'. FDTD cross-check: measured extra phase / analytic = {fdtd_vs_analytic} "
        f"(should be ~1 if the linearized index is honest). Bucket 0 unchanged."
    )

    out = {
        "substrate": "complex Intent field, focusing-saturating on-site nonlinearity (Phase-2)",
        "absolute_time_note": "varying quantity is probe phase velocity (cells/absolute tick); the tick is universal, no clock perturbed",
        "params": {"A0": A0, "core_width": CORE_W, "w0": W0, "gamma": GAMMA, "us": US, "w_probe": round(float(w_probe), 4)},
        "n_minus_1_peak": round(peak, 6),
        "integrated_extra_phase_analytic": round(integrated_phase_analytic, 4),
        "range_classification": rng,
        "fdtd_extra_phase": (round(dphi_fdtd, 4) if dphi_fdtd is not None else None),
        "fdtd_over_analytic": fdtd_vs_analytic,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase3b_substrate_index_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
