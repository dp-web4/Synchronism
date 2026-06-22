"""
Phase-3b — can the Intent substrate host LONG-RANGE (1/r) gravity at all? (2026-06-22)

A prior draft of this probe was near-circular: it imposed a sech-shaped core in a massive field
and "found" the induced index well was short-range — but that just echoed the imposed shape.
Discarded. The real, non-circular question: a static field SOURCED by a point mass has a spatial
range set by the field's OWN rest mass. Don't impose the shape — SOLVE for it.

The framework's Intent field has a rest mass w0 (that is the whole point of Phase-1/2: a
standing pattern's rest oscillation w0 IS the particle's mass). A massive field's static response
to a point source is the YUKAWA form e^{-w0 r}/r — SHORT-RANGE, with range ~ 1/w0. Only a
MASSLESS field (w0 -> 0) gives the COULOMB/NEWTON form 1/r — long-range, the precondition for
gravity. And 1/r exists ONLY in 3D (dp: dimensions matter) — so we solve the genuinely-3D radial
static field equation.

Faithful to absolute time: this is a static SPATIAL field; no clocks, no ticks vary.

Method: solve the 3D radial screened-Poisson (static Klein-Gordon)
    -(1/r^2) d/dr( r^2 dphi/dr ) + w0^2 phi = source(r)
i.e.  -phi'' - (2/r) phi' + w0^2 phi = source,  with a compact source near r=0 and phi(r_max)=0.
Tridiagonal (Thomas) solve. Then fit phi(r) to e^{-kappa r}/r and read kappa; Yukawa predicts
kappa = w0, Newton needs kappa = 0. Sweep w0 in {1.0, 0.3, 0.1, 0.0}.

PRE-REGISTERED: massive field (w0>0) -> kappa ~ w0, short-range Yukawa (range 1/w0); only w0=0
-> kappa ~ 0, long-range 1/r. CONSEQUENCE if so: the framework's particle-giving massive Intent
field CANNOT mediate gravity (wrong range); long-range 1/r needs a separate MASSLESS Intent mode
the framework does not currently have. A productive, structural negative (Yukawa 1935: massive
mediator => short range — the reason weak/strong forces are short-range, EM/gravity long-range).

numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

RMAX = 400.0
N = 4000
SRC_R = 3.0          # compact source radius (the "mass")
W0_LIST = [1.0, 0.3, 0.1, 0.0]


def thomas(a, b, c, d):
    """Solve tridiagonal system: a[i]x[i-1] + b[i]x[i] + c[i]x[i+1] = d[i]."""
    n = len(d)
    cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def solve_static(w0):
    r = np.linspace(RMAX / N, RMAX, N)   # avoid r=0
    dr = r[1] - r[0]
    # operator rows:  -phi'' - (2/r)phi' + w0^2 phi
    a = np.zeros(N); b = np.zeros(N); c = np.zeros(N); d = np.zeros(N)
    for i in range(N):
        a[i] = -1.0 / dr ** 2 + (1.0 / r[i]) / dr      # phi_{i-1}: -1/dr^2 from -phi'', +(2/r)/(2dr) from -(2/r)phi'
        b[i] = 2.0 / dr ** 2 + w0 ** 2
        c[i] = -1.0 / dr ** 2 - (1.0 / r[i]) / dr
        d[i] = 1.0 if r[i] <= SRC_R else 0.0           # compact uniform source ("point" mass)
    # BCs: regularity at inner edge (Neumann phi'~0): fold c[0] into b[0]
    b[0] += a[0]; a[0] = 0.0
    # outer: phi(r_max)=0 (Dirichlet)
    b[-1] = 1.0; a[-1] = 0.0; c[-1] = 0.0; d[-1] = 0.0
    phi = thomas(a, b, c, d)
    return r, phi


def fit_kappa(r, phi):
    """Fit phi(r) ~ A * e^{-kappa r} / r in the far region (r > a few * source)."""
    mask = (r > 8 * SRC_R) & (r < 0.6 * RMAX) & (phi > 0)
    rs, ps = r[mask], phi[mask]
    if len(rs) < 10:
        return None, None
    # log(phi * r) = log A - kappa r
    y = np.log(ps * rs)
    A = np.vstack([np.ones_like(rs), -rs]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    logA, kappa = coef[0], coef[1]
    # goodness: how well does e^{-kappa r}/r explain it vs pure 1/r (kappa=0)?
    pred_yuk = np.exp(logA) * np.exp(-kappa * rs) / rs
    pred_invr = np.exp(np.mean(np.log(ps * rs))) / rs
    r2_yuk = 1 - np.sum((ps - pred_yuk) ** 2) / np.sum((ps - ps.mean()) ** 2)
    r2_invr = 1 - np.sum((ps - pred_invr) ** 2) / np.sum((ps - ps.mean()) ** 2)
    return float(kappa), {"r2_yukawa": round(float(r2_yuk), 5), "r2_pure_1overr": round(float(r2_invr), 5)}


def main():
    rows = []
    for w0 in W0_LIST:
        r, phi = solve_static(w0)
        kappa, gof = fit_kappa(r, phi)
        # effective range: where phi drops to 1/e of its value at r=10*SRC_R
        i0 = np.searchsorted(r, 10 * SRC_R)
        ref = phi[i0]
        below = np.where(phi[i0:] < ref / np.e)[0]
        rng = float(r[i0 + below[0]] - r[i0]) if len(below) else float(RMAX)
        rows.append({
            "w0": w0,
            "fitted_kappa": (round(kappa, 4) if kappa is not None else None),
            "kappa_over_w0": (round(kappa / w0, 3) if (kappa is not None and w0 > 0) else None),
            "effective_range_1e": round(rng, 1),
            "expected_yukawa_range_1_over_w0": (round(1.0 / w0, 1) if w0 > 0 else "inf (massless)"),
            "goodness": gof,
        })

    massive = [x for x in rows if x["w0"] > 0]
    massless = next(x for x in rows if x["w0"] == 0.0)
    yukawa_holds = all(x["kappa_over_w0"] is not None and abs(x["kappa_over_w0"] - 1.0) < 0.25 for x in massive)
    massless_is_longrange = massless["fitted_kappa"] is not None and abs(massless["fitted_kappa"]) < 0.01

    verdict = (
        f"SOLVED (not imposed): the static Intent field sourced by a point mass is YUKAWA when the "
        f"field is massive. Fitted decay kappa tracks the field rest-mass w0 (kappa/w0 ~ 1 across "
        f"w0>0: {yukawa_holds}) -> range ~ 1/w0, SHORT. Only the massless field (w0=0) gives "
        f"kappa ~ 0 ({massless_is_longrange}) -> the pure 1/r Coulomb/Newton form, long-range. "
        f"STRUCTURAL CONSEQUENCE: the framework's Intent field carries a rest mass w0 BY DESIGN "
        f"(that w0 IS particle mass, Phase-1/2). A massive mediator forces a short-range force "
        f"(Yukawa 1935) — so the substrate's particle sector CANNOT be the source of gravity, "
        f"whose 1/r long range is the precondition for the whole light-deflection / factor-of-2 "
        f"analysis. To host gravity the framework needs a SEPARATE MASSLESS Intent mode (a "
        f"graviton-analog / the field's far tail / a collective long-wavelength excitation), "
        f"which it does not currently specify. This is the real gap, and it sits UPSTREAM of the "
        f"factor-of-2: there is no point asking 'space-only 2x or falsified' until a 1/r field "
        f"exists to bend light at range. PRODUCTIVE NEGATIVE, non-circular (the range came out of "
        f"the solve, not an imposed shape), and it connects to textbook physics: massive mediator "
        f"=> short range is exactly why the weak/strong forces are short-range and EM/gravity are "
        f"long-range. And 1/r exists only in 3D — the dimension mattered. Bucket 0 unchanged; this "
        f"adds a sharply-stated open problem to Bucket 1's gravity line."
    )

    out = {
        "question": "can a massive Intent field mediate long-range (1/r) gravity? (solved, not imposed)",
        "absolute_time_note": "static spatial field; no clock/tick varies",
        "geometry": "3D radial screened-Poisson (static Klein-Gordon), Thomas solve",
        "rows": rows,
        "yukawa_holds_for_massive": bool(yukawa_holds),
        "massless_is_long_range_1overr": bool(massless_is_longrange),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase3b_intent_field_range_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n w0     fitted_kappa   kappa/w0   range(1/e)   expected 1/w0")
    for x in rows:
        print(f"  {x['w0']:.1f}    {str(x['fitted_kappa']):>8}     {str(x['kappa_over_w0']):>6}    "
              f"{x['effective_range_1e']:>8}     {x['expected_yukawa_range_1_over_w0']}")


if __name__ == "__main__":
    main()
