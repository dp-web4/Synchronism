"""
Phase-3 — 3D light propagation & gravitational deflection: the factor-of-2 discriminator
(2026-06-22, CBP-Claude, autonomous — dp: 1D/2D are napkin sketches with impenetrable-wall
artifacts; reality is 3D; degrees of freedom matter. What keeps photons linear in 3D, and why
does gravity bend them, weakly but nonzero?)

PROPOSAL under test:
  Q1 (linearity): in a HOMOGENEOUS, isotropic substrate a propagating mode goes straight by
     symmetry (translational invariance -> momentum conservation; eikonal ray in a uniform
     medium = straight line, Fermat). Straightness is the default, the absence of a gradient.
  Q2 (gravity bends it): mass = a localized high-Intent core that makes the substrate locally
     inhomogeneous — it lowers the local reconstruction rate c (an effective refractive index
     n>1). A gradient in c bends the ray toward the slower region (gradient-index optics).

THE SHARP TEST — the factor of 2 (the most famous discriminator in the history of gravity):
  n(r) = 1 + alpha * (GM/c^2) / r.
  Deflection (weak field): Δθ = 2*alpha*(GM/c^2)/b.
    alpha = 1  -> Δθ = 2GM/(c^2 b)  = Newton/Soldner = a SCALAR "variable-c only" substrate
    alpha = 2  -> Δθ = 4GM/(c^2 b)  = General Relativity (Eddington 1919, observed)
  A substrate that models mass as ONLY a reconstruction-rate (c) reduction gives alpha=1 —
  HALF the observed bending — and is FALSIFIED, exactly as scalar gravity was. To match GR it
  must reduce BOTH the tick-rate (time) AND the grid spacing (space), equally: that is the
  origin of the 2. This experiment confirms the alpha-scaling and the factor of 2 by direct 3D
  eikonal ray-tracing (genuinely 3D — no missing-dimension walls), and isolates the substrate-
  structure requirement.

Geometric units: c = 1, GM = 1, so the gravitational radius r_g = GM/c^2 = 1; weak field is
b >> 1. RK4 integration of the eikonal ray equation in 3D:
    dr/ds = t ;  d(n t)/ds = grad(n)  =>  dt/ds = (grad n - (grad n . t) t) / n.
numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

C = 1.0
GM = 1.0
RG = GM / (C * C)        # = 1 (gravitational radius)
B_LIST = [20.0, 50.0, 100.0, 200.0]
ALPHAS = [0.0, 1.0, 2.0]   # 0 = homogeneous (straight), 1 = scalar/Newton, 2 = GR
DS = 0.2
X0 = 4000.0              # start/stop |x| (>> b: asymptotically flat)


def n_and_grad(pos, alpha):
    r = np.linalg.norm(pos) + 1e-12
    n = 1.0 + alpha * RG / r
    gradn = -alpha * RG * pos / r ** 3
    return n, gradn


def deriv(state, alpha):
    pos = state[:3]; t = state[3:]
    n, gradn = n_and_grad(pos, alpha)
    dt = (gradn - np.dot(gradn, t) * t) / n
    return np.concatenate([t, dt])


def trace_ray(b, alpha):
    pos = np.array([-X0, b, 0.0])
    t = np.array([1.0, 0.0, 0.0])
    state = np.concatenate([pos, t])
    steps = 0
    max_steps = int(4 * X0 / DS)
    rmin = np.inf
    while state[0] < X0 and steps < max_steps:
        k1 = deriv(state, alpha)
        k2 = deriv(state + 0.5 * DS * k1, alpha)
        k3 = deriv(state + 0.5 * DS * k2, alpha)
        k4 = deriv(state + DS * k3, alpha)
        state = state + (DS / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        tnorm = np.linalg.norm(state[3:])
        state[3:] /= tnorm                       # keep tangent unit (eikonal hygiene)
        rmin = min(rmin, np.linalg.norm(state[:3]))
        steps += 1
    tx, ty = state[3], state[4]
    deflection = float(np.arctan2(-ty, tx))      # incoming along +x (ty=0); bend toward mass
    return abs(deflection), float(rmin)


def main():
    results = []
    for alpha in ALPHAS:
        for b in B_LIST:
            dtheta, rmin = trace_ray(b, alpha)
            analytic = 2.0 * alpha * RG / b
            results.append({
                "alpha": alpha, "b": b,
                "deflection_meas": round(dtheta, 6),
                "deflection_analytic_2aRG_over_b": round(analytic, 6),
                "ratio_meas_over_analytic": (round(dtheta / analytic, 3) if analytic > 0 else None),
                "dtheta_times_b": round(dtheta * b, 4),     # should ~ 2*alpha*RG = 2*alpha
                "closest_approach": round(rmin, 1),
            })

    # Q1: homogeneous substrate (alpha=0) -> straight line (deflection ~ 0)
    straight = [r for r in results if r["alpha"] == 0.0]
    homogeneous_is_straight = all(r["deflection_meas"] < 1e-4 for r in straight)
    # factor of 2: deflection(alpha=2)/deflection(alpha=1) at matched b
    ratios = []
    for b in B_LIST:
        d1 = next(r["deflection_meas"] for r in results if r["alpha"] == 1.0 and r["b"] == b)
        d2 = next(r["deflection_meas"] for r in results if r["alpha"] == 2.0 and r["b"] == b)
        ratios.append(d2 / d1)
    factor = float(np.mean(ratios))
    # alpha=1 matches Newton/Soldner, alpha=2 matches GR (within ray-trace accuracy)
    a1_is_newton = all(abs(r["dtheta_times_b"] - 2.0) < 0.15 for r in results if r["alpha"] == 1.0)
    a2_is_GR = all(abs(r["dtheta_times_b"] - 4.0) < 0.30 for r in results if r["alpha"] == 2.0)

    verdict = (
        f"Q1 confirmed: a HOMOGENEOUS substrate (alpha=0) propagates the ray dead straight "
        f"(deflection < 1e-4: {homogeneous_is_straight}) — linearity is the default, by "
        f"translational symmetry; no active mechanism needed. Q2 confirmed: a local c-gradient "
        f"(mass) bends the ray, weakly, toward the slower region (gradient-index). THE FACTOR "
        f"OF 2: deflection(alpha=2)/deflection(alpha=1) = {factor:.3f} across impact parameters. "
        f"A scalar 'variable-c only' substrate (alpha=1) reproduces Newton/Soldner 2GM/(c^2 b) "
        f"(matches: {a1_is_newton}) — HALF the observed bending, the value Eddington 1919 ruled "
        f"out. Matching GR's 4GM/(c^2 b) (alpha=2; matches: {a2_is_GR}) requires DOUBLING the "
        f"effect, which physically means the substrate must perturb BOTH the reconstruction rate "
        f"(time / c) AND the grid spacing (space), in equal measure. So the classic test hands "
        f"the substrate a sharp, falsifiable requirement: a mass cannot be modeled as only a "
        f"slowdown of Intent reconstruction (that gives half the bending and is dead on arrival) "
        f"— it must equally contract the grid. THE OPEN QUESTION (the real frontier): does the "
        f"framework's structure (c as reconstruction rate + the two-level time + the grid "
        f"itself) NATURALLY produce the factor-of-2 (space and time perturbed equally), or must "
        f"the 2 be put in by hand? If natural -> the substrate explains gravitational lensing. "
        f"If by hand -> it is fitting GR, not deriving it. Not yet answered; this is where to "
        f"keep proposing. NOT novel physics (the n=1+2GM/c^2r effective-medium picture of GR is "
        f"textbook); the value is locating the exact substrate requirement and the exact open "
        f"question. Bucket 0 unchanged."
    )

    out = {
        "units": "c=1, GM=1, r_g=1 (weak field b>>1)", "ds": DS, "x_range": X0,
        "alpha_meaning": {"0": "homogeneous (no mass)", "1": "scalar / variable-c-only (Newton/Soldner)",
                          "2": "GR-equivalent effective index"},
        "rays": results,
        "homogeneous_is_straight": bool(homogeneous_is_straight),
        "factor_of_2_measured": round(factor, 3),
        "alpha1_matches_newton": bool(a1_is_newton),
        "alpha2_matches_GR": bool(a2_is_GR),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase3_light_deflection_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n alpha   b      deflection   analytic(2aRG/b)  Δθ·b   closest_approach")
    for r in results:
        print(f"  {r['alpha']:.0f}   {r['b']:6.0f}   {r['deflection_meas']:.6f}    "
              f"{r['deflection_analytic_2aRG_over_b']:.6f}      {r['dtheta_times_b']:.3f}   {r['closest_approach']:.1f}")
    print(f"\n factor of 2  (alpha=2 / alpha=1) = {factor:.3f}")


if __name__ == "__main__":
    main()
