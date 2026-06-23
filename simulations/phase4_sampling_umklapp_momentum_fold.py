"""
Phase-4 — the universe-as-sampler: does a discrete tick rate produce Umklapp (momentum mod G)?
(2026-06-22)

THE FRAME (dp 2026-06-22). Not cell-properties yet — the stage before: a GRID OF CELLS
ADVANCING STATE ON GLOBAL CLOCK TICKS. The question: if the substrate is a discrete-time state
machine sampled at a base rate, are quantization / Planck-scale / high-energy effects what you'd
EXPECT to fall out? This rests on a fractal-leverage bet: a Planck-scale discrete-time grid, if it
is a state machine, should behave *meaningfully similarly* to a digital signal processor — the
similarity holds within the discrete-sampling MRH, which is the MRH under test. Nyquist and
Umklapp are INVITED analogies (MRH-appropriate), not smuggled. The one thing NOT imported: any
"cutoff co-dilates with patterns" move — the sampling cutoff is the grid's property, pattern-
unaware, by the model's current statement.

THE TEST. A sampler has a band limit (spatial: |k| <= pi/a, the first Brillouin zone = the
Nyquist wavenumber). A LINEAR signal below it samples faithfully. A NONLINEAR interaction MIXES:
patterns at k1, k2 generate sum content at k1+k2. If k1+k2 exceeds the band limit, the sampler
cannot represent it — it ALIASES. On a grid, aliasing of crystal-momentum IS Umklapp: the daughter
appears at (k1+k2 - G), G = 2pi/a (a reciprocal-lattice vector). DIRECTIONAL momentum needs a
COMPLEX field e^{ikx} (a real cos has a sign-symmetric spectrum and cannot show direction); the
same Phase-1.6 lesson. With a complex field the fold is unambiguous: a forward daughter (k1+k2>pi)
reappears moving BACKWARD at k1+k2-2pi.

WHAT THIS PREDICTS FOR REALITY (if the fractal bet holds). Continuum QFT vacuum has continuous
translation symmetry -> EXACT momentum conservation (Noether). A discrete substrate has only
discrete translation symmetry -> momentum conserved mod G ~ Planck momentum -> "vacuum Umklapp":
extreme-momentum interactions show apparent momentum non-conservation by G / anomalous back-
scatter. A distinct LIV channel from Phase-2's dispersion-LIV, and (dp) reachable at lower
per-quantum energy because mixing STACKS momentum toward the cutoff.

THE CONTROL (agent-zero). Refine the grid (smaller a -> larger Nyquist G). The SAME physical
interaction must then NOT fold. If refinement kills the fold, it is a genuine SAMPLING-CUTOFF
effect, not coded-in physics.

numpy only. Headless. Complex field, chi^2 (quadratic) nonlinearity = minimal sum-frequency
generator; folding is nonlinearity-agnostic.
"""
import json
import os
import numpy as np


def wrap(k):
    """Fold a wavenumber into the first Brillouin zone (-pi, pi] — the sampler's band limit."""
    return (k + np.pi) % (2 * np.pi) - np.pi


def peak_k(field, d=1.0):
    """Wavenumber of the dominant spectral bin of a complex field (signed = directional)."""
    spec = np.abs(np.fft.fft(field))
    spec[0] = 0.0
    kf = 2 * np.pi * np.fft.fftfreq(len(field), d=d)
    return kf[int(np.argmax(spec))]


def static_complex_shg(N=512):
    """Pump = single complex tone e^{ikx}; chi^2 product (ψ^2) = e^{2ikx}. Sweep k across the
    half-zone pi/2 and show the daughter's SIGNED wavenumber folds to wrap(2k) — reversing to
    negative k (backward) once 2k>pi. Directional Umklapp, unambiguous (complex field)."""
    x = np.arange(N)
    rows = []
    for kf in np.linspace(0.15, 0.95, 17):
        bin_k = int(round(kf * np.pi * N / (2 * np.pi)))
        k = 2 * np.pi * bin_k / N
        psi = np.exp(1j * k * x)
        daughter = peak_k(psi * psi)                # chi^2: ψ^2 = e^{2ikx}, sampled -> aliased bin
        rows.append({
            "k_over_pi": round(k / np.pi, 3),
            "naive_2k_over_pi": round(2 * k / np.pi, 3),
            "measured_daughter_over_pi": round(daughter / np.pi, 3),
            "wrap2k_pred_over_pi": round(wrap(2 * k) / np.pi, 3),
            "folds": bool(2 * k > np.pi),
            "reverses_direction": bool(daughter < -1e-6 and k > 0),  # forward pump -> backward daughter
        })
    return rows


def two_pump_collision(N=512, a=1.0, refine=1):
    """Two complex pumps at PHYSICAL wavenumbers q1,q2; chi^2 makes the sum-frequency daughter at
    q1+q2. Grid spacing a/refine sets the band limit pi/(a/refine). q1+q2 chosen ABOVE the coarse
    band limit so it folds; refine raises the band limit (control -> no fold)."""
    a_eff = a / refine
    M = N * refine
    x = np.arange(M) * a_eff
    q1, q2 = 1.7, 1.8                 # sum = 3.5 > pi (~3.14) -> folds on the coarse grid
    band = np.pi / a_eff
    psi = np.exp(1j * q1 * x) + np.exp(1j * q2 * x)
    prod = psi * psi                  # cross term 2 e^{i(q1+q2)x} dominates the non-pump content
    spec = np.abs(np.fft.fft(prod)); spec[0] = 0.0
    kf = 2 * np.pi * np.fft.fftfreq(M, d=a_eff)
    # exclude the self-terms 2q1, 2q2 to isolate the sum-frequency daughter q1+q2
    mask = (np.abs(kf - 2 * q1) > 0.1) & (np.abs(kf - 2 * q2) > 0.1)
    idx = np.where(mask)[0]
    daughter = kf[idx[int(np.argmax(spec[idx]))]]
    folded_pred = wrap((q1 + q2) * a_eff) / a_eff
    return {
        "refine": refine, "a_eff": round(a_eff, 4), "band_limit": round(band, 4),
        "q1": q1, "q2": q2, "sum": round(q1 + q2, 4),
        "sum_exceeds_band": bool(q1 + q2 > band),
        "naive_daughter": round(q1 + q2, 4),
        "folded_pred": round(folded_pred, 4),
        "measured_daughter": round(float(daughter), 4),
        "is_umklapp_fold": bool(abs(daughter - folded_pred) < 0.05 and (q1 + q2) > band),
    }


def main():
    static = static_complex_shg()
    coarse = two_pump_collision(refine=1)
    fine = two_pump_collision(refine=4)

    fold_exact = all(abs(r["measured_daughter_over_pi"] - r["wrap2k_pred_over_pi"]) < 0.02 for r in static)
    reversals = [r for r in static if r["folds"]]
    all_folds_reverse = all(r["reverses_direction"] for r in reversals)
    control_kills_fold = coarse["is_umklapp_fold"] and not fine["is_umklapp_fold"]

    verdict = (
        f"DISCRETE SAMPLER -> UMKLAPP, demonstrated, directional, and controlled. (1) Static "
        f"complex chi^2 sweep: the daughter's SIGNED wavenumber = wrap(2k) at every k "
        f"(exact: {fold_exact}); once 2k crosses the half-zone pi/2 the daughter REVERSES to "
        f"negative k (forward pump -> backward daughter: {all_folds_reverse}). The sampler cannot "
        f"hold above-band content, so the mixing product folds — momentum conserved only MOD "
        f"G=2pi/a. (2) Two-pump collision: q1+q2={coarse['sum']} > coarse band {coarse['band_limit']} "
        f"-> daughter at {coarse['measured_daughter']} = folded {coarse['folded_pred']} (not naive "
        f"{coarse['naive_daughter']}); umklapp={coarse['is_umklapp_fold']}. (3) CONTROL "
        f"(agent-zero): refine 4x (band {coarse['band_limit']} "
        f"-> {fine['band_limit']}) and the SAME interaction does NOT fold (daughter "
        f"{fine['measured_daughter']} = unfolded {fine['naive_daughter']}); refine-kills-fold="
        f"{control_kills_fold}. So the fold is a genuine SAMPLING-CUTOFF effect, vanishing as the "
        f"Nyquist limit rises — not coded-in physics. "
        f"TESTABLE PREDICTION (what the analogy buys): IF the universe is a discrete-time grid "
        f"(the fractal Planck-DSP bet), vacuum momentum is conserved only mod G ~ Planck momentum "
        f"-> 'vacuum Umklapp': extreme-momentum-transfer interactions show apparent momentum "
        f"non-conservation by G / anomalous back-scatter — a distinct LIV channel from Phase-2's "
        f"dispersion-LIV, reachable at lower per-quantum energy because mixing STACKS momentum "
        f"toward the cutoff (dp). HONESTY: Umklapp is textbook condensed-matter (momentum mod G); "
        f"NOT novel. The novel falsifiable claim is its application to the VACUUM substrate "
        f"(momentum mod a Planck-scale G), currently UNOBSERVED (momentum conservation is exact to "
        f"high precision) -> a risky Bucket-1 bet, Planck-suppressed for single interactions, NOT a "
        f"confirmation. Bucket 0 unchanged. All of it is conditional on the fractal bet that a "
        f"Planck grid behaves like a sampler — the hypothesis under test, stated as such."
    )

    out = {
        "frame": "grid-of-cells advancing on global ticks, tested AS a discrete-time sampler (DSP); fractal bet stated, not assumed",
        "static_complex_shg_sweep": static,
        "two_pump_coarse": coarse,
        "two_pump_fine_control": fine,
        "fold_exact_wrap2k": bool(fold_exact),
        "folds_reverse_direction": bool(all_folds_reverse),
        "control_refine_kills_fold": bool(control_kills_fold),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase4_sampling_umklapp_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n  k/pi   naive 2k/pi   measured daughter/pi   wrap(2k)/pi   (FOLD = reverses)")
    for r in static:
        print(f"  {r['k_over_pi']:.3f}    {r['naive_2k_over_pi']:.3f}        "
              f"{r['measured_daughter_over_pi']:+.3f}                {r['wrap2k_pred_over_pi']:+.3f}     "
              f"{'FOLD<-' if r['folds'] else ''}")


if __name__ == "__main__":
    main()
