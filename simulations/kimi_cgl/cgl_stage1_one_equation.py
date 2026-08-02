"""
K1 of the CGL arc (charter: explorations/2026-08-01-kimi-cgl-arc-plan.md) — ONE complex
Ginzburg-Landau parameter family vs the Phase-1 Stage-1 battery.

  dA/dt = A + (1+ib) lap(A) - (1+ic)|A|^2 A        (A complex, 1D periodic lattice, L=256)

The bet (pre-registered): the repo's two substrates — dissipative diffusion (Phase-1 arm A,
0% pass) and conservative focusing breather (arm D, 83% pass) — are regimes of ONE CGL
equation, i.e. some (b,c) region passes the Stage-1 bar (>10% of 24 random ICs produce
localized + genuinely oscillating + perturbation-surviving states) while the b,c -> 0 corner
reproduces diffusion-like failure. A phase boundary between the two is NOT a kill.

Pre-registered kill criteria (charter, registered 2026-08-01 before any run):
  K1a: NO (b,c) point passes the >10% bar -> kill.
  K1b: the breather behavior cannot be reproduced in kind by any (1+ic)|A|^2 A choice
       (best CGL pass rate far below the in-harness D-arm anchor) -> kill.

Harness mirrors simulations/phase1_stage1_localized_oscillation_substrate.py:
  - L=256 periodic lattice, 24 seeded random Gaussian-pulse ICs, perturbation kick at
    midpoint, same localization/oscillation classifier on |A|, same >10% criterion.
  - Arms A and D of the original battery are PORTED VERBATIM into this file and run as
    calibration anchors in the same classifier (expected: A ~0%, D high).

CLASSIFIER EXTENSION (documented choice, per the K1 brief): the phase-1 oscillation test is
direction reversals of the core AMPLITUDE. A uniformly phase-rotating CGL localized state
A = rho(x) exp(-i w t) has constant |A| yet is physically oscillating — scoring it
"non-oscillating" would repeat the arm-A metric flaw in reverse. So for the complex field we
ADD a parallel test: net unwrapped core-phase travel >= 4*pi over the post-perturbation
window (with amplitude retained) also counts as oscillating. Both sub-flags are logged per
run (osc_amp / osc_phase), and pass rates under the STRICT amplitude-only phase-1 criterion
are reported alongside (n_pass_amp_only) so the extension is auditable and the anchors
(phase-1 real fields: phase travel = 0) remain exactly comparable.

Integration: explicit RK4 on the lattice. CFL/stability reasoning: two stiff modes —
(i) the Laplacian's band edge |lambda|max = 4/dx^2 = 4 (dx=1), multiplied by (1+ib);
(ii) the local saturation term, whose radial linearization at amplitude r is
~-(1+ic)(3r^2-1), stiffness ~3|1+ic|r^2. RK4 is stable for |lambda*dt| <~ 2.78 on the
negative real axis (~2.83 on the imaginary axis). A b-only rule dt = 0.6/sqrt(1+b^2)
covers (i) but VIOLATES (ii) at c=7 and blew up in the smoke test, so dt must cover both:
  dt = 2.4 / (4*|1+ib| + 12*|1+ic|)
keeps each contribution <~1.2 < 2.78 with margin: the diffusive band edge 4|1+ib|dt <= 2.4
by construction, and the saturation term stays stable up to transient overshoot |A| ~ 2
(3*4*|1+ic|*dt = 12|1+ic|dt <= 2.4). The linear-growth time scale (e^t) is well resolved.
Runs that still go non-finite are flagged "blew_up" and logged as such (numerical event,
not a physical regime) rather than silently absorbed.

numpy + stdlib only; seeded (seed=1); headless. Writes simulations/results/kimi_cgl/
cgl_stage1_sweep.json. Full run is the default; --n-runs / --b-grid / --c-grid / --steps
allow a reduced run (state the reduction in the JSON's "run_note").
"""
import argparse
import json
import os
import time

import numpy as np

L = 256               # lattice sites (1D, periodic), mirrors phase-1
N_RUNS = 24           # random ICs per (b,c) point, mirrors phase-1
SEED = 1
T_PRE = 2000          # RK4 steps before the perturbation kick
T_POST = 2000         # steps after the kick
SAMPLE_EVERY = 50     # mirrors phase-1 sampling cadence
PASS_FRACTION_CRITERION = 0.10   # Stage-1 bar, unchanged
PHASE_TRAVEL_CRITERION = 4 * np.pi   # core-phase oscillation test (documented above)

# (b,c) coarse grid, 5x5. Log-ish spacing spanning the three regimes named in the brief:
#   b,c -> 0       : near-diffusion limit (arm-A-like behavior expected)
#   b,c ~ 1.5-3    : 1D CGL localized-state / Nozaki-Bekki hole territory
#   b,c large (7)  : NLS-like, conservative-dominant limit (arm-D-like behavior expected)
B_GRID_DEFAULT = [0.0, 0.5, 1.5, 3.0, 7.0]
C_GRID_DEFAULT = [0.0, 0.5, 1.5, 3.0, 7.0]


# ---------------------------------------------------------------- CGL core

def laplacian(A):
    """Periodic 5-point-free 3-point lattice Laplacian, dx=1."""
    return np.roll(A, -1) - 2.0 * A + np.roll(A, 1)


def make_rhs(b, c):
    def rhs(A):
        return A + (1.0 + 1j * b) * laplacian(A) - (1.0 + 1j * c) * np.abs(A) ** 2 * A
    return rhs


def rk4_step(A, dt, rhs):
    k1 = rhs(A)
    k2 = rhs(A + 0.5 * dt * k1)
    k3 = rhs(A + 0.5 * dt * k2)
    k4 = rhs(A + dt * k3)
    return A + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def seed_pulse_complex(rng):
    """Random localized IC, mirrors phase-1 seed_pulse: random width/amplitude/position,
    real Gaussian pulse + small complex noise (zero initial phase gradient)."""
    x = np.arange(L)
    pos = rng.integers(L // 4, 3 * L // 4)
    width = rng.uniform(3, 10)
    amp = rng.uniform(0.4, 0.95)
    A = amp * np.exp(-0.5 * ((x - pos) / width) ** 2)
    A = A + 0.01 * (rng.standard_normal(L) + 1j * rng.standard_normal(L))
    return A.astype(complex)


def eff_width(u):
    """Effective occupied sites: (sum u^2)^2 / sum u^4. Small = localized; ~L = spread."""
    s2 = np.sum(u * u)
    s4 = np.sum(u ** 4)
    if s4 < 1e-12:
        return float(L)
    return float(s2 * s2 / s4)


def run_cgl_point(rng, b, c, n_runs, t_pre, t_post):
    """One (b,c) point: n_runs seeded ICs, kick at midpoint, classify each run."""
    # CFL: covers both stiff modes — diffusive band edge 4|1+ib| and saturation stiffness
    # 3|1+ic||A|^2 (overshoot to |A|~2 allowed) — each contribution <= 2.4 < 2.78 (RK4).
    dt = 2.4 / (4.0 * np.hypot(1.0, b) + 12.0 * np.hypot(1.0, c))
    rhs = make_rhs(b, c)
    runs = [run_cgl_single(rng, rhs, dt, t_pre, t_post) for _ in range(n_runs)]
    return summarize_runs(b, c, dt, runs)


def run_cgl_single(rng, rhs, dt, t_pre, t_post):
    A = seed_pulse_complex(rng)
    amps, widths, phases = [], [], []
    core_site = None   # fixed at the post-window argmax site, for the phase-travel test
    blew_up = False
    for t in range(t_pre + t_post):
        if t == t_pre:
            A = A + 0.05 * (rng.standard_normal(L) + 1j * rng.standard_normal(L))  # kick
        A = rk4_step(A, dt, rhs)
        if not np.all(np.isfinite(A)):
            blew_up = True
            break
        if t % SAMPLE_EVERY == 0:
            a = np.abs(A)
            amps.append(float(np.max(a)))
            widths.append(eff_width(a))
            if t >= t_pre:
                if core_site is None:
                    core_site = int(np.argmax(a))
                phases.append(float(np.angle(A[core_site])))
    res = classify_complex(amps, widths, phases)
    res["blew_up"] = blew_up
    return res


def classify_complex(amps, widths, phases):
    """Phase-1 classifier on |A| (localized + amplitude-reversal oscillation), PLUS the
    documented core-phase-travel oscillation test for the complex field. Both sub-flags
    logged; pass = localized AND (amp-reversals OR phase travel), amplitude retained."""
    amps = np.array(amps)
    widths = np.array(widths)
    n = len(amps)
    if n < 8:   # blew up before meaningful sampling — numerical event, not a pass
        return {
            "pass": False, "pass_amp_only": False, "localized": False,
            "osc_amp": False, "osc_phase": False, "reversals": 0,
            "phase_travel_rad": 0.0, "final_width": float(L), "spread_ratio": 0.0,
            "amp_retained": 0.0,
        }
    post = slice(n // 2, None)
    apost = amps[post]
    wpost = widths[post]
    last = slice(-max(4, len(apost) // 3), None)
    early = slice(0, max(4, len(wpost) // 3))
    final_width = float(np.mean(wpost[last]))
    early_width = float(np.mean(wpost[early]))
    final_amp = float(np.mean(apost[last]))
    init_amp = float(np.mean(amps[:max(1, n // 10)]))
    localized = (final_width < 0.15 * L) and (final_width <= 1.15 * early_width)
    d = np.diff(apost)
    s = np.sign(d)
    s = s[s != 0]
    reversals = int(np.sum(s[1:] != s[:-1])) if len(s) > 1 else 0
    amp_ok = final_amp > 0.25 * init_amp
    osc_amp = (reversals >= 4) and amp_ok
    phase_travel = float(np.abs(np.unwrap(np.array(phases))[-1] - np.unwrap(np.array(phases))[0])) \
        if len(phases) > 1 else 0.0
    osc_phase = (phase_travel >= PHASE_TRAVEL_CRITERION) and amp_ok
    passed = bool(localized and (osc_amp or osc_phase))
    return {
        "pass": passed,
        "pass_amp_only": bool(localized and osc_amp),   # strict phase-1 criterion
        "localized": bool(localized),
        "osc_amp": bool(osc_amp),
        "osc_phase": bool(osc_phase),
        "reversals": reversals,
        "phase_travel_rad": round(phase_travel, 1),
        "final_width": round(final_width, 1),
        "spread_ratio": round(final_width / (early_width + 1e-9), 2),
        "amp_retained": round(final_amp / (init_amp + 1e-9), 3),
    }


def summarize_runs(b, c, dt, runs):
    n = len(runs)
    n_pass = sum(r["pass"] for r in runs)
    n_pass_amp_only = sum(r["pass_amp_only"] for r in runs)
    frac = n_pass / n
    mean_amp = float(np.mean([r["amp_retained"] for r in runs]))
    mean_width = float(np.mean([r["final_width"] for r in runs]))
    frac_localized = sum(r["localized"] for r in runs) / n
    frac_osc = sum(r["osc_amp"] or r["osc_phase"] for r in runs) / n
    n_blew_up = sum(r.get("blew_up", False) for r in runs)
    return {
        "b": b, "c": c, "dt": round(dt, 4),
        "pass_fraction": round(frac, 3), "n_pass": n_pass, "n_runs": n,
        "n_pass_amp_only": n_pass_amp_only,
        "n_blew_up": n_blew_up,
        "stage1_pass": bool(frac > PASS_FRACTION_CRITERION),
        "mean_amp_retained": round(mean_amp, 3),
        "mean_final_width": round(mean_width, 1),
        "frac_localized": round(frac_localized, 3),
        "frac_oscillating": round(frac_osc, 3),
        "regime": regime_tag(frac, mean_amp, mean_width, frac_osc, n_blew_up / n),
        "example": runs[0],
    }


def regime_tag(pass_frac, mean_amp_retained, mean_width, frac_osc, frac_blew_up=0.0):
    """Simple observable-based tag: decay / breather-like / turbulence / diffusion-like.
    >50% numerical blow-up overrides everything: that is an integration event, not physics."""
    if frac_blew_up > 0.5:
        return "numerically-unstable"
    if mean_amp_retained < 0.05:
        return "decay"
    if pass_frac > PASS_FRACTION_CRITERION:
        return "breather-like"
    if mean_width > 0.5 * L:
        return "turbulence" if frac_osc > 0.3 else "diffusion-like"
    return "diffusion-like"


# ------------------------------------------------- phase-1 anchors (ported verbatim)

DT_ANCHOR = 0.05
C2 = 1.0
OMEGA0_2 = 1.0
U_MAX = 1.0
N_SAT = 2
GAMMA = 1.2
U_S = 0.6
T_PRE_ANCHOR = 4000
T_POST_ANCHOR = 4000


def monotonic_R(absu):
    return np.clip(1.0 - (np.minimum(absu, U_MAX) / U_MAX) ** N_SAT, 0.0, 1.0)


def lap_coupled(u, Rfield):
    Rl = 0.5 * (Rfield + np.roll(Rfield, -1))
    right = Rl * (np.roll(u, -1) - u)
    left = np.roll(Rl, 1) * (np.roll(u, 1) - u)
    return right + left


def seed_pulse_real(rng):
    x = np.arange(L)
    pos = rng.integers(L // 4, 3 * L // 4)
    width = rng.uniform(3, 10)
    amp = rng.uniform(0.4, 0.95)
    u = amp * np.exp(-0.5 * ((x - pos) / width) ** 2)
    u += 0.01 * rng.standard_normal(L)
    return u


def core_amp(u):
    return float(np.max(np.abs(u)))


def classify_real(amps, widths, energy_drift):
    """Phase-1 classifier, verbatim (real fields => no phase channel)."""
    amps = np.array(amps)
    widths = np.array(widths)
    n = len(amps)
    post = slice(n // 2, None)
    apost = amps[post]
    wpost = widths[post]
    last = slice(-max(4, len(apost) // 3), None)
    early = slice(0, max(4, len(wpost) // 3))
    final_width = float(np.mean(wpost[last]))
    early_width = float(np.mean(wpost[early]))
    final_amp = float(np.mean(apost[last]))
    init_amp = float(np.mean(amps[:max(1, n // 10)]))
    localized = (final_width < 0.15 * L) and (final_width <= 1.15 * early_width)
    d = np.diff(apost)
    s = np.sign(d)
    s = s[s != 0]
    reversals = int(np.sum(s[1:] != s[:-1])) if len(s) > 1 else 0
    oscillating = (reversals >= 4) and (final_amp > 0.25 * init_amp)
    passed = bool(localized and oscillating)
    return {
        "pass": passed, "localized": bool(localized), "oscillating": bool(oscillating),
        "reversals": reversals,
        "final_width": round(final_width, 1),
        "spread_ratio": round(final_width / (early_width + 1e-9), 2),
        "amp_retained": round(final_amp / (init_amp + 1e-9), 3),
        "energy_drift": None if energy_drift is None else round(energy_drift, 4),
    }


def run_arm_A(rng):
    """Anchor A: first-order diffusion + monotonic R (ported verbatim from phase-1)."""
    u = seed_pulse_real(rng)
    D = 0.2
    amps, widths = [], []
    for t in range(T_PRE_ANCHOR + T_POST_ANCHOR):
        if t == T_PRE_ANCHOR:
            u = u + 0.05 * rng.standard_normal(L)
        R = monotonic_R(np.abs(u))
        u = u + DT_ANCHOR * D * lap_coupled(u, R)
        if t % 50 == 0:
            amps.append(core_amp(u))
            widths.append(eff_width(u))
    return classify_real(amps, widths, energy_drift=None)


def run_arm_D(rng):
    """Anchor D: 2nd-order wave + focusing-saturating on-site (ported verbatim)."""
    u = seed_pulse_real(rng)
    v = np.zeros(L)

    def accel(u):
        coupling = C2 * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
        onsite = -OMEGA0_2 * u + GAMMA * (u ** 3) / (1.0 + (u / U_S) ** 2)
        return coupling + onsite

    def energy(u, v):
        kin = 0.5 * np.sum(v * v)
        grad = 0.5 * C2 * np.sum((np.roll(u, -1) - u) ** 2)
        pot = 0.5 * OMEGA0_2 * np.sum(u * u)
        return float(kin + grad + pot)

    a = accel(u)
    e0 = energy(u, v)
    e_max_dev = 0.0
    amps, widths = [], []
    for t in range(T_PRE_ANCHOR + T_POST_ANCHOR):
        if t == T_PRE_ANCHOR:
            v = v + 0.05 * rng.standard_normal(L)
            e0 = energy(u, v)
        u = u + DT_ANCHOR * v + 0.5 * DT_ANCHOR * DT_ANCHOR * a
        a_new = accel(u)
        v = v + 0.5 * DT_ANCHOR * (a + a_new)
        a = a_new
        if t % 50 == 0:
            amps.append(core_amp(u))
            widths.append(eff_width(u))
            e_max_dev = max(e_max_dev, abs(energy(u, v) - e0) / (abs(e0) + 1e-9))
    return classify_real(amps, widths, energy_drift=e_max_dev)


def run_anchor(name, fn, rng, n_runs):
    runs = [fn(rng) for _ in range(n_runs)]
    n_pass = sum(r["pass"] for r in runs)
    frac = n_pass / n_runs
    drifts = [r["energy_drift"] for r in runs if r["energy_drift"] is not None]
    return {
        "arm": name,
        "pass_fraction": round(frac, 3), "n_pass": n_pass, "n_runs": n_runs,
        "stage1_pass": bool(frac > PASS_FRACTION_CRITERION),
        "mean_final_width": round(float(np.mean([r["final_width"] for r in runs])), 1),
        "mean_amp_retained": round(float(np.mean([r["amp_retained"] for r in runs])), 3),
        "max_energy_drift": (round(max(drifts), 4) if drifts else None),
    }


# ------------------------------------------------------------------ sweep + path

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-runs", type=int, default=N_RUNS)
    ap.add_argument("--b-grid", type=float, nargs="+", default=B_GRID_DEFAULT)
    ap.add_argument("--c-grid", type=float, nargs="+", default=C_GRID_DEFAULT)
    ap.add_argument("--steps", type=int, default=T_PRE,
                    help="steps pre AND post kick (total = 2x)")
    ap.add_argument("--path-steps", type=int, default=8)
    ap.add_argument("--note", type=str, default="full run (defaults)")
    args = ap.parse_args()

    t0 = time.time()
    rng = np.random.default_rng(SEED)

    # 1) calibration anchors in THIS harness (ported phase-1 arms A and D)
    anchors = [
        run_anchor("A: 1st-order diffusion + monotonic R", run_arm_A, rng, args.n_runs),
        run_anchor("D: 2nd-order wave + focusing-saturating on-site", run_arm_D, rng, args.n_runs),
    ]

    # 2) (b,c) coarse grid
    grid = []
    for b in args.b_grid:
        for c in args.c_grid:
            grid.append(run_cgl_point(rng, b, c, args.n_runs, args.steps, args.steps))
            print(f"[grid] b={b:4.1f} c={c:4.1f} pass={grid[-1]['pass_fraction']:.3f} "
                  f"{grid[-1]['regime']}", flush=True)

    # 3) path check: best breather-like point -> most diffusion-like point, straight line
    best = max(grid, key=lambda g: g["pass_fraction"])
    diff_candidates = [g for g in grid if g["regime"] == "diffusion-like"]
    if diff_candidates:
        most_diff = min(diff_candidates, key=lambda g: (g["pass_fraction"], -g["mean_final_width"]))
    else:
        most_diff = min(grid, key=lambda g: g["pass_fraction"])
    path = []
    n_path = max(2, args.path_steps)
    for i in range(n_path):
        f = i / (n_path - 1)
        b = (1 - f) * best["b"] + f * most_diff["b"]
        c = (1 - f) * best["c"] + f * most_diff["c"]
        pt = run_cgl_point(rng, b, c, args.n_runs, args.steps, args.steps)
        pt["path_fraction_along"] = round(f, 3)
        path.append(pt)
        print(f"[path] f={f:.2f} b={b:5.2f} c={c:5.2f} pass={pt['pass_fraction']:.3f} "
              f"{pt['regime']}", flush=True)

    # 4) pre-registered kill-criteria evaluation (coarse data only; verdict at the arc wake)
    best_cgl = max(g["pass_fraction"] for g in grid)
    anchor_a = anchors[0]["pass_fraction"]
    anchor_d = anchors[1]["pass_fraction"]
    k1a_kill = best_cgl <= PASS_FRACTION_CRITERION
    k1b_kill = best_cgl < 0.5 * anchor_d   # "far below the D-arm anchor" -> flagged, coarse
    out = {
        "stage": "K1 — one equation, two regimes (CGL (b,c) sweep vs Stage-1 battery)",
        "charter": "explorations/2026-08-01-kimi-cgl-arc-plan.md",
        "equation": "dA/dt = A + (1+ib) lap(A) - (1+ic)|A|^2 A  (1D periodic, L=256)",
        "lattice": L, "n_runs": args.n_runs, "seed": SEED,
        "t_pre_steps": args.steps, "t_post_steps": args.steps,
        "sample_every": SAMPLE_EVERY,
        "dt_rule": "dt = 2.4/(4*|1+ib| + 12*|1+ic|)  (RK4: covers diffusive band edge 4|1+ib| and saturation stiffness 3|1+ic||A|^2 to |A|~2; each term <= 2.4 < 2.78)",
        "pass_fraction_criterion": PASS_FRACTION_CRITERION,
        "phase_travel_criterion_rad": round(PHASE_TRAVEL_CRITERION, 3),
        "classifier_note": (
            "phase-1 classifier on |A| (localized + >=4 core-amplitude direction reversals), "
            "PLUS documented core-phase-travel (>=4pi unwrapped, post-kick, amplitude retained) "
            "as an oscillation channel for the complex field; n_pass_amp_only gives the strict "
            "phase-1 criterion. Anchors are real fields: phase channel identically zero."
        ),
        "b_grid": list(args.b_grid), "c_grid": list(args.c_grid),
        "grid_rationale": (
            "log-ish 5x5: (0,0) near-diffusion limit; (1.5,1.5)-(3,3) known 1D CGL "
            "localized-state territory; (7,7) NLS-like conservative-dominant limit."
        ),
        "run_note": args.note,
        "anchors": anchors,
        "grid": grid,
        "path": {
            "from_best_breather": {"b": best["b"], "c": best["c"],
                                   "pass_fraction": best["pass_fraction"]},
            "to_most_diffusion": {"b": most_diff["b"], "c": most_diff["c"],
                                  "pass_fraction": most_diff["pass_fraction"],
                                  "regime": most_diff["regime"]},
            "points": path,
        },
        "kill_criteria_coarse_eval": {
            "K1a_no_point_above_bar": bool(k1a_kill),
            "K1b_best_far_below_D_anchor": bool(k1b_kill),
            "best_cgl_pass_fraction": best_cgl,
            "anchor_A_pass_fraction": anchor_a,
            "anchor_D_pass_fraction": anchor_d,
            "note": "coarse-grid flag only; the verdict belongs to the arc wake's analysis.",
        },
        "runtime_seconds": round(time.time() - t0, 1),
    }
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "kimi_cgl")
    os.makedirs(outdir, exist_ok=True)
    path_out = os.path.join(outdir, "cgl_stage1_sweep.json")
    with open(path_out, "w") as f:
        json.dump(out, f, indent=2)

    print("\n=== anchors (this harness) ===")
    for a in anchors:
        print(f"  {a['arm'][:52]:52} pass={a['pass_fraction']:.3f}")
    print("\n=== (b,c) grid ===")
    print("  b \\ c | " + "  ".join(f"{c:>7.1f}" for c in args.c_grid))
    for b in args.b_grid:
        row = [g for g in grid if g["b"] == b]
        print(f"  {b:5.1f} | " + "  ".join(f"{g['pass_fraction']:7.3f}" for g in row))
    print(f"\nK1a kill (coarse): {k1a_kill}   K1b kill (coarse): {k1b_kill}")
    print(f"best CGL pass: {best_cgl:.3f}  anchor A: {anchor_a:.3f}  anchor D: {anchor_d:.3f}")
    print(f"runtime: {out['runtime_seconds']}s")
    print(f"wrote {path_out}")


if __name__ == "__main__":
    main()
