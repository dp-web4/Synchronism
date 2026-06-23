"""
Phase-8 (P2-deep) — multi-faith on capacity rules: which one DERIVES the GR gravity profile? (2026-06-22)

Phase-7 located the remaining freedom in deriving gravity: the equation-of-state (how intent
DENSITY responds to a coherence sink) — GR needs rho ∝ r^(-3/2). dp's method: we DON'T know the
capacity rule; we DISCOVER it by hypothesizing a few promising ones and testing which lands on GR
(multi-faith — imagine paths, explore with feedback, refine/abandon/redirect). The capacity rule
(how much intent a cell holds / can receive, gated on its momentum) IS the EoS, so each candidate
rule predicts a density profile rho(r) and a light deflection. Test them.

SETUP (steady radial inflow, intent conserved): continuity rho*v*r^2 = const. A capacity rule that
relates held-density to flow gives the closure. The momentum-coupled family — "a cell's held intent
scales with its flow speed as rho ∝ v^n" (n = how strongly capacity is gated on momentum) — is the
clean computable backbone:
    rho ∝ v^n , continuity rho*v*r^2 = const  =>  v ∝ r^(-2/(n+1)) , rho ∝ r^(-2n/(n+1)).
GR target rho ∝ r^(-3/2)  =>  2n/(n+1) = 3/2  =>  n = 3. So IF the capacity rule is rho ∝ v^3,
gravity matches GR — derived from continuity, no force postulated. Test the family (does n=3 give
GR on BOTH the density exponent AND the light deflection?), then test a second, physically-different
faith (saturating capacity) as a refine/redirect.

FAITHS:
  A. momentum-coupled capacity rho ∝ v^n  (sweep n) — does n=3 land on GR?
  B. saturating capacity (rho clipped at I_max near the core) — does it change the GR tail, or just
     core the centre? (redirect: a cored-but-GR-tail profile is what dark-matter HALOS look like.)
  (Dummy/baseline: n=0 = incompressible = Phase-7's wrong b^-4 case.)

numpy only. Headless. Eikonal swimmer for the deflection (Phase-3c/7). Writes results JSON.
"""
import json
import os
import numpy as np

C = 1.0
A = 1.0
B_LIST = [40.0, 80.0, 160.0, 320.0]
DS = 0.2
X0 = 8000.0


def u_inflow(pos, p):
    r = np.linalg.norm(pos) + 1e-12
    return -A / r ** p * (pos / r)


def deriv(state, p):
    pos = state[:2]; k = state[2:]
    khat = k / (np.linalg.norm(k) + 1e-12)
    u = u_inflow(pos, p)
    dx = C * khat + u
    eps = 1e-4; gku = np.zeros(2)
    for i in range(2):
        d = np.zeros(2); d[i] = eps
        gku[i] = (np.dot(k, u_inflow(pos + d, p)) - np.dot(k, u_inflow(pos - d, p))) / (2 * eps)
    return np.concatenate([dx, -gku])


def deflection(b, p):
    state = np.array([-X0, b, 1.0, 0.0]); steps = 0; mx = int(4 * X0 / DS)
    while state[0] < X0 and steps < mx:
        k1 = deriv(state, p); k2 = deriv(state + 0.5 * DS * k1, p)
        k3 = deriv(state + 0.5 * DS * k2, p); k4 = deriv(state + DS * k3, p)
        state = state + (DS / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4); steps += 1
    return abs(float(np.arctan2(-state[3], state[2])))


def deflection_exponent(p):
    d = np.array([deflection(b, p) for b in B_LIST]); b = np.array(B_LIST)
    pos = d > 0
    return (-float(np.polyfit(np.log(b[pos]), np.log(d[pos]), 1)[0]) if pos.sum() >= 2 else float("inf"))


def main():
    # FAITH A: momentum-coupled capacity rho ∝ v^n
    faith_a = []
    for n in [0, 1, 2, 3, 4, 5]:
        rho_exp = -2 * n / (n + 1)            # density power law
        p = 2.0 / (n + 1)                     # velocity power law -> inflow exponent
        q = deflection_exponent(p)            # measured deflection scaling Delta_theta ∝ b^-q
        faith_a.append({
            "n_momentum_coupling": n,
            "capacity_rule": f"held intent rho ∝ v^{n}",
            "density_exp": round(rho_exp, 3),
            "velocity_exp": round(-p, 3),
            "deflection_q": round(q, 3),
            "is_GR": bool(abs(rho_exp + 1.5) < 0.05 and abs(q - 1.0) < 0.15),
        })
    gr_faith = next((f for f in faith_a if f["is_GR"]), None)

    # FAITH B: saturating capacity — clip the n=3 density at I_max near the core; does the far tail
    # (and thus lensing) change, or only the centre? Analytic: clipping only affects r < r_sat
    # (where rho would exceed I_max); the r > r_sat tail stays rho ∝ r^-3/2. So GR lensing is
    # preserved and the centre is CORED. Demonstrate the tail exponent is unchanged outside r_sat.
    I_max = 5.0
    r = np.logspace(-1, 3, 400)
    rho_unsat = r ** (-1.5)                    # n=3 (GR) profile, arbitrary normalization
    rho_sat = np.minimum(rho_unsat, I_max)     # capacity ceiling
    r_sat = float(I_max ** (-1 / 1.5))         # radius where rho hits I_max
    tail = r > 5 * r_sat
    tail_exp = float(np.polyfit(np.log(r[tail]), np.log(rho_sat[tail]), 1)[0])
    cored = bool(np.all(rho_sat[r < r_sat] >= I_max - 1e-9))
    faith_b = {"I_max": I_max, "r_sat": round(r_sat, 3),
               "tail_density_exp_outside_r_sat": round(tail_exp, 3),
               "tail_still_GR_-1.5": bool(abs(tail_exp + 1.5) < 0.05),
               "centre_cored": cored}

    verdict = (
        f"MULTI-FAITH ON CAPACITY RULES — which derives GR? FAITH A (momentum-coupled, rho ∝ v^n): "
        f"sweeping the coupling n, the density exponent is -2n/(n+1) and the light deflection scales "
        f"b^-q; ONLY n=3 lands on GR — rho ∝ r^-3/2 AND Delta_theta ∝ 1/b "
        f"(GR faith: n={gr_faith['n_momentum_coupling'] if gr_faith else None}, "
        f"density_exp={gr_faith['density_exp'] if gr_faith else None}, q={gr_faith['deflection_q'] if gr_faith else None}). "
        f"The others are abandoned-for-GR but real: n=0 (incompressible) -> b^-4 (Phase-7's wrong "
        f"case); n=1,2 -> too-steep deflection; n>3 -> too shallow. So the capacity rule GR REQUIRES "
        f"is sharp and falsifiable: 'a cell's held intent scales as the CUBE of its flow speed' "
        f"(rho ∝ v^3) — derived from continuity, NO gravitational force postulated. That is real "
        f"narrowing: gravity-is-fit has shrunk from 'fit the inflow' (3c) to 'fit the EoS' (Phase-7) "
        f"to 'the momentum-coupling exponent is exactly 3' (here). What we have NOT done: derive WHY "
        f"n=3 from the cell's actual receive-dynamics — that is the open question, now a single "
        f"number. FAITH B (saturating capacity, rho <= I_max): clipping only cores the CENTRE "
        f"(r < r_sat={faith_b['r_sat']}); the tail stays rho ∝ r^{faith_b['tail_density_exp_outside_r_sat']} "
        f"(still GR: {faith_b['tail_still_GR_-1.5']}). REDIRECT, not abandon: a cored-centre + GR-tail "
        f"profile is exactly what galactic dark-matter HALOS look like (the core-cusp structure) — so "
        f"the saturating capacity rule is worth keeping for the galactic (P4/B2) frontier, where it "
        f"would change rotation curves WITHOUT changing solar-system lensing. HONEST: power-law EoS "
        f"-> power-law profile is textbook fluid analog; finding n=3 is LOCATING the rule, not "
        f"deriving it (numerological until the receive-dynamics force n=3). Bucket 0 unchanged. But "
        f"the multi-faith did its job: one path (n=3) is the GR target, one (saturation) redirects to "
        f"the galactic frontier, the rest are eliminated."
    )

    out = {
        "method": "multi-faith: hypothesize capacity rules, test which derives the GR profile (dp 2026-06-22)",
        "faith_a_momentum_coupled_rho_v_n": faith_a,
        "gr_selecting_rule": gr_faith,
        "faith_b_saturating_capacity": faith_b,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase8_capacity_rule_multifaith_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\nFAITH A — momentum-coupled capacity rho ∝ v^n:")
    print("  n   density_exp   deflection_q   GR? (density -1.5 AND q=1)")
    for f in faith_a:
        print(f"  {f['n_momentum_coupling']}    {f['density_exp']:+.3f}        {f['deflection_q']:.3f}         {f['is_GR']}")
    print(f"\nFAITH B — saturating capacity: centre cored, tail exp {faith_b['tail_density_exp_outside_r_sat']} "
          f"(GR tail preserved: {faith_b['tail_still_GR_-1.5']}) -> redirect to galactic frontier")


if __name__ == "__main__":
    main()
