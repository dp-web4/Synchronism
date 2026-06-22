"""
Phase-1 / CA-challenge Stage 1 — can a local discrete-grid rule produce a STABLE,
LOCALIZED, OSCILLATING pattern that survives perturbation? And critically: does the
framework's OWN monotonic-saturation mechanism do it, or only a focusing ingredient that
breaks Foundation 3?

Grounded in the documented nulls (do NOT reinvent):
- S617/S666: the first-order transfer rule ∂I/∂t = ∇·[D·R(I)·∇I] is 1-DOF scalar diffusion
  — dissipative, cannot oscillate.
- S19/S20/S665: a second-order wave with monotonic saturation R(I)=1-(|I|/I_max)^n in the
  coupling is DEFOCUSING (local wave speed c·√R decreases at high I → high-I center lags,
  low-I edges run ahead → the pulse DISPERSES (width 13→81) and pulses REPEL). "Any
  monotonically-decreasing R(I) is defocusing; self-confinement impossible from R(I) alone."
- Open obligation (saturation-reframe-corrections-2026-05-28): the additional ingredient
  beyond independent J that could escape the dispersal is a FOCUSING nonlinearity — which
  "breaks the saturation-as-pattern-stability axiom (Foundation 3)."

THE TEST. Four local rules on a 1D lattice (periodic), each seeded with random localized
pulses, evolved, perturbed mid-run, evolved again. A run PASSES Stage 1 if its final state
is (a) localized (effective width ≪ L), (b) oscillating (core amplitude keeps varying, not
decayed), and (c) survived the perturbation. CA-challenge criterion: a rule family passes
if > 10% of random ICs pass.

  ARM A  first-order diffusion + monotonic R   — the literal Foundation-1/3 substrate (S617 null)
  ARM B  second-order wave   + monotonic-R coupling — Foundation-3-faithful (S19 defocusing null)
  ARM C  second-order wave   + linear              — control (linear waves disperse)
  ARM D  second-order wave   + focusing-saturating on-site — the flagged escape ingredient

EXPECTED (honest pre-registration): A, B, C fail (no localized oscillation); D passes.
If so, the result is a genuine CONSTRAINT, not a win: stable particle-like patterns are
achievable on the discrete grid, but NOT from monotonic saturation (the framework's own
Foundation 3) — they require a focusing nonlinearity, which costs that axiom. That moves the
CA challenge from "can it oscillate at all" (no, for the framework's rule) to "what must the
substrate give up to make particles" (Foundation 3). Discrete breathers in a focusing KG
lattice are textbook — so ARM D passing is NOT novel physics; it is a precise localization
of what the framework's substrate lacks. Reported in the four-bucket honesty frame.

Velocity-Verlet (symplectic) for the wave arms so any decay is physical, not numerical;
energy drift is reported as the honesty check. numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

L = 256              # lattice sites (1D, periodic)
DT = 0.05
C2 = 1.0             # (wave speed)^2 base
OMEGA0_2 = 1.0       # on-site linear restoring (mass term)
U_MAX = 1.0          # saturation scale
N_SAT = 2            # saturation exponent for monotonic R = 1-(|u|/U_MAX)^n
GAMMA = 1.2          # focusing strength (arm D)
U_S = 0.6            # focusing saturation scale (arm D)
T_PRE = 4000         # steps before perturbation
T_POST = 4000        # steps after perturbation
N_RUNS = 24          # random ICs per arm
SEED = 1
PASS_FRACTION_CRITERION = 0.10   # CA-challenge Stage-1 threshold


def monotonic_R(absu):
    """Framework Foundation-3 saturation: resistance falls as |u| -> U_MAX (defocusing)."""
    return np.clip(1.0 - (np.minimum(absu, U_MAX) / U_MAX) ** N_SAT, 0.0, 1.0)


def lap_coupled(u, Rfield):
    """Discrete Laplacian with amplitude-dependent coupling R (periodic).
    Flux between neighbours uses the mean R of the pair (conservative)."""
    Rl = 0.5 * (Rfield + np.roll(Rfield, -1))   # coupling on bond j..j+1
    right = Rl * (np.roll(u, -1) - u)
    left = np.roll(Rl, 1) * (np.roll(u, 1) - u)
    return right + left


def seed_pulse(rng):
    """Random localized initial pulse: random width, amplitude, position, + small noise."""
    x = np.arange(L)
    pos = rng.integers(L // 4, 3 * L // 4)
    width = rng.uniform(3, 10)
    amp = rng.uniform(0.4, 0.95)
    u = amp * np.exp(-0.5 * ((x - pos) / width) ** 2)
    u += 0.01 * rng.standard_normal(L)
    return u


def eff_width(u):
    """Effective number of occupied sites: (Σu²)² / Σu⁴ . Small = localized; ~L = spread."""
    s2 = np.sum(u * u)
    s4 = np.sum(u ** 4)
    if s4 < 1e-12:
        return float(L)
    return float(s2 * s2 / s4)


def core_amp(u):
    return float(np.max(np.abs(u)))


def run_first_order_diffusion(rng):
    """ARM A: I_t = D * d/dx[R(I) dI/dx]. Dissipative; expected no oscillation."""
    u = seed_pulse(rng)
    D = 0.2
    amps, widths = [], []
    for t in range(T_PRE + T_POST):
        if t == T_PRE:
            u = u + 0.05 * rng.standard_normal(L)  # perturbation kick
        R = monotonic_R(np.abs(u))
        u = u + DT * D * lap_coupled(u, R)
        if t % 50 == 0:
            amps.append(core_amp(u)); widths.append(eff_width(u))
    return classify(amps, widths, energy_drift=None)


def run_wave(rng, mode):
    """ARM B/C/D: second-order wave via velocity-Verlet.
    mode = 'monoR'  -> monotonic-R coupling, linear on-site   (Foundation-3 faithful)
         = 'linear' -> constant coupling, linear on-site       (control)
         = 'focus'  -> constant coupling, focusing-saturating on-site (escape ingredient)"""
    u = seed_pulse(rng)
    v = np.zeros(L)

    def accel(u):
        if mode == 'monoR':
            R = monotonic_R(np.abs(u))
            coupling = C2 * lap_coupled(u, R)
            onsite = -OMEGA0_2 * u
        elif mode == 'linear':
            coupling = C2 * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
            onsite = -OMEGA0_2 * u
        elif mode == 'focus':
            coupling = C2 * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
            # focusing (soft) on-site, saturated to avoid blow-up:
            # F = -w0^2 u + gamma u^3/(1+(u/u_s)^2)  -> self-attraction balances dispersion
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
    for t in range(T_PRE + T_POST):
        if t == T_PRE:
            v = v + 0.05 * rng.standard_normal(L)  # perturbation kick (energy injection)
            e0 = energy(u, v)                       # reset reference post-kick
        # velocity Verlet
        u = u + DT * v + 0.5 * DT * DT * a
        a_new = accel(u)
        v = v + 0.5 * DT * (a + a_new)
        a = a_new
        if t % 50 == 0:
            amps.append(core_amp(u)); widths.append(eff_width(u))
            e_max_dev = max(e_max_dev, abs(energy(u, v) - e0) / (abs(e0) + 1e-9))
    return classify(amps, widths, energy_drift=e_max_dev)


def classify(amps, widths, energy_drift):
    """Honest classifier. NOTE: an earlier version mis-scored Arm A (diffusion) as a pass —
    monotonic decay has nonzero amplitude-variance (looked like 'oscillation') and a
    still-spreading pulse can sit just under the width threshold (looked like 'localized').
    The diffusion arm is the built-in dummy: it provably cannot oscillate (S617), so it
    caught the metric flaw. Fixed here: oscillation = repeated DIRECTION REVERSALS of the
    core amplitude (peaks AND troughs), not mere variance; localized = small AND not still
    spreading across the post-perturbation window."""
    amps = np.array(amps); widths = np.array(widths)
    n = len(amps)
    post = slice(n // 2, None)                      # post-perturbation half
    apost = amps[post]; wpost = widths[post]
    last = slice(-max(4, len(apost) // 3), None)
    early = slice(0, max(4, len(wpost) // 3))
    final_width = float(np.mean(wpost[last]))
    early_width = float(np.mean(wpost[early]))
    final_amp = float(np.mean(apost[last]))
    init_amp = float(np.mean(amps[:max(1, n // 10)]))
    # localized: well below box size AND not still spreading (kills diffusion false-positive)
    localized = (final_width < 0.15 * L) and (final_width <= 1.15 * early_width)
    # oscillating: core amplitude reverses direction repeatedly (genuine peaks+troughs).
    # Monotonic decay -> ~0 reversals. Breather -> many.
    d = np.diff(apost)
    s = np.sign(d); s = s[s != 0]
    reversals = int(np.sum(s[1:] != s[:-1])) if len(s) > 1 else 0
    oscillating = (reversals >= 4) and (final_amp > 0.25 * init_amp)
    passed = bool(localized and oscillating)
    return {
        "pass": passed, "localized": bool(localized), "oscillating": bool(oscillating),
        "reversals": reversals,
        "final_width": round(final_width, 1), "spread_ratio": round(final_width / (early_width + 1e-9), 2),
        "amp_retained": round(final_amp / (init_amp + 1e-9), 3),
        "energy_drift": None if energy_drift is None else round(energy_drift, 4),
    }


def run_arm(name, fn, rng):
    runs = [fn(rng) for _ in range(N_RUNS)]
    n_pass = sum(r["pass"] for r in runs)
    frac = n_pass / N_RUNS
    drifts = [r["energy_drift"] for r in runs if r["energy_drift"] is not None]
    return {
        "arm": name,
        "pass_fraction": round(frac, 3),
        "n_pass": n_pass, "n_runs": N_RUNS,
        "stage1_pass": bool(frac > PASS_FRACTION_CRITERION),
        "mean_final_width": round(float(np.mean([r["final_width"] for r in runs])), 1),
        "mean_amp_retained": round(float(np.mean([r["amp_retained"] for r in runs])), 3),
        "max_energy_drift": (round(max(drifts), 4) if drifts else None),
        "example": runs[0],
    }


def main():
    rng = np.random.default_rng(SEED)
    arms = [
        ("A: 1st-order diffusion + monotonic R (Foundation-1/3 substrate; S617 null)",
         run_first_order_diffusion),
        ("B: 2nd-order wave + monotonic-R coupling (Foundation-3 faithful; S19 defocusing)",
         lambda r: run_wave(r, 'monoR')),
        ("C: 2nd-order wave + linear (control)",
         lambda r: run_wave(r, 'linear')),
        ("D: 2nd-order wave + focusing-saturating on-site (the escape ingredient)",
         lambda r: run_wave(r, 'focus')),
    ]
    results = [run_arm(name, fn, rng) for name, fn in arms]

    foundation3_passes = results[1]["stage1_pass"]    # monotonic saturation builds particles?
    focusing_passes = results[3]["stage1_pass"]       # focusing ingredient builds particles?

    if focusing_passes and not foundation3_passes:
        verdict = (
            "CA-challenge Stage 1 is ACHIEVABLE on the discrete grid — but NOT from the "
            "framework's own monotonic saturation (Foundation 3), which disperses (reproduces "
            "S19/S665). Stable localized oscillating patterns require a FOCUSING nonlinearity, "
            "which breaks Foundation 3's saturation-as-pattern-stability axiom. This is a "
            "CONSTRAINT, not a confirmation: it precisely names what the substrate must give up "
            "to make particles. (Discrete breathers in a focusing lattice are textbook — Arm D "
            "passing is not novel physics; it localizes what the framework's rule lacks.) "
            "Next: can Arm-D patterns INTERACT / gain mass-like and charge-like behavior "
            "(Stage 2+)? Open."
        )
    elif foundation3_passes:
        verdict = (
            "Unexpected: the monotonic-saturation rule (Foundation 3) produced localized "
            "oscillation > criterion. Contradicts S19/S665 — AUDIT before any claim "
            "(check for numerical artifact / boundary trapping / insufficient runtime)."
        )
    else:
        verdict = (
            "No arm passed Stage 1 above criterion. Even the focusing ingredient did not "
            "self-confine under these parameters — sweep (gamma, u_s, omega0, runtime) before "
            "concluding the rule family fails; record as a tighter null."
        )

    out = {
        "lattice": L, "dt": DT, "t_pre": T_PRE, "t_post": T_POST, "n_runs": N_RUNS,
        "seed": SEED, "pass_fraction_criterion": PASS_FRACTION_CRITERION,
        "arms": results,
        "foundation3_monotonic_saturation_passes": bool(foundation3_passes),
        "focusing_ingredient_passes": bool(focusing_passes),
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results",
                        "phase1_stage1_localized_oscillation_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n arm                                              pass_frac  final_width  amp_kept  E-drift")
    for r in results:
        print(f"  {r['arm'][:46]:46}  {r['pass_fraction']:.2f}       "
              f"{r['mean_final_width']:6.1f}     {r['mean_amp_retained']:.2f}     {r['max_energy_drift']}")


if __name__ == "__main__":
    main()
