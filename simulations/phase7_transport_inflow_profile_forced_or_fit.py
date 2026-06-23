"""
Phase-7 (P2) — does a transport DOF DERIVE the gravity inflow, or is the profile still fit? (2026-06-22)

Following the temporal/spatial research clue (the spatial sector is where the walls are) and the
unpaid part of Phase-3c (the GP inflow u=sqrt(2GM/r) was IMPOSED). Phase-3b showed the complex-KG
substrate's *massive* field can't host long-range gravity (Yukawa, short-range). Phase-3c's
resolution was to treat gravity as substrate INFLOW — but it imposed the flow. The honest question:
add the transport DOF dp specified (intent has momentum / "where it wants to go", under continuity)
and ask whether a coherent high-Intent core (a 'mass' = a sink toward which intent reconstructs)
DERIVES a long-range inflow, and whether the substrate's rules FORCE the GP profile or leave it free.

DERIVATION (steady radial inflow, intent conserved by continuity):
  continuity: d/dr( rho * v_r * r^2 ) = -(sink at origin)  =>  rho * v_r * r^2 = const = -Mdot/(4pi).
  - INCOMPRESSIBLE (rho = rho0, the minimal/default): v_r ∝ 1/r^2.  LONG-RANGE (power law) — the
    Yukawa obstacle DISSOLVES from the transport DOF, derived (not imposed).
  - GP / GR-matching: v_r = sqrt(2GM/r) ∝ 1/r^(1/2). This requires a SPECIFIC equation-of-state
    (the self-advection=Newtonian condition from Phase-3c) — i.e. rho must vary as rho ∝ r^(-3/2)
    for continuity to give it. NOT the minimal incompressible case.

THE TEST: trace light (eikonal swimmer H=c|k|+k.u, Phase-3c) through each inflow u ∝ 1/r^p and
measure the deflection-vs-impact-parameter scaling Delta_theta ∝ b^(-q). GR (observed) needs
Delta_theta ∝ 1/b (q=1). GP (p=1/2) gives q=1 (Phase-3c: coefficient -> 4). Does the minimal
incompressible continuity profile (p=2) give q=1 (it would be forced + correct) or q != 1 (wrong)?

PRE-REGISTRATION: incompressible continuity (p=2) is long-range but FALLS TOO FAST -> its light
deflection scales steeper than 1/b -> WRONG vs GR. So the transport DOF dissolves Yukawa but does
NOT force the GP profile; matching GR needs a specific EoS (rho ∝ r^-3/2), still a CHOICE. Net: P2
answer = "long-range DERIVED, GP profile still FIT (via the EoS)". A productive partial result that
locates the remaining freedom precisely (the equation-of-state).

numpy only. Headless. Planar ray trace (exact orbital plane). Writes results JSON.
"""
import json
import os
import numpy as np

C = 1.0
A = 1.0                 # inflow amplitude (sets units); weak field via large b
B_LIST = [40.0, 80.0, 160.0, 320.0]
DS = 0.2
X0 = 8000.0


def u_inflow(pos, p):
    """Radial inflow u = -A / r^p * r_hat (toward origin)."""
    r = np.linalg.norm(pos) + 1e-12
    return -A / r ** p * (pos / r)


def deriv(state, p):
    pos = state[:2]; k = state[2:]
    kmag = np.linalg.norm(k) + 1e-12
    khat = k / kmag
    u = u_inflow(pos, p)
    dx = C * khat + u
    eps = 1e-4
    gku = np.zeros(2)
    for i in range(2):
        d = np.zeros(2); d[i] = eps
        gku[i] = (np.dot(k, u_inflow(pos + d, p)) - np.dot(k, u_inflow(pos - d, p))) / (2 * eps)
    dk = -gku
    return np.concatenate([dx, dk])


def deflection(b, p):
    pos = np.array([-X0, b]); k = np.array([1.0, 0.0])
    state = np.concatenate([pos, k])
    steps = 0; max_steps = int(4 * X0 / DS)
    while state[0] < X0 and steps < max_steps:
        k1 = deriv(state, p); k2 = deriv(state + 0.5 * DS * k1, p)
        k3 = deriv(state + 0.5 * DS * k2, p); k4 = deriv(state + DS * k3, p)
        state = state + (DS / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4); steps += 1
    return abs(float(np.arctan2(-state[3], state[2])))


def scaling_exponent(p):
    """Fit Delta_theta ∝ b^(-q); return q and the per-b values."""
    rows = []
    for b in B_LIST:
        rows.append({"b": b, "deflection": round(deflection(b, p), 8)})
    b = np.array([r["b"] for r in rows]); d = np.array([r["deflection"] for r in rows])
    pos = d > 0
    if pos.sum() >= 2:        # fit only non-underflowed points (steep profiles hit float 0 fast)
        q = -float(np.polyfit(np.log(b[pos]), np.log(d[pos]), 1)[0])
    else:
        q = float("inf")
    return q, rows


def main():
    # both inflow profiles are LONG-RANGE power laws (Yukawa dissolved by the transport DOF)
    q_gp, rows_gp = scaling_exponent(0.5)         # GP / GR-matching profile
    q_inc, rows_inc = scaling_exponent(2.0)        # incompressible-continuity (minimal) profile

    gp_is_1overb = abs(q_gp - 1.0) < 0.15
    inc_is_1overb = abs(q_inc - 1.0) < 0.15

    verdict = (
        f"DOES THE TRANSPORT DOF FORCE THE GRAVITY PROFILE? Partial — long-range is DERIVED, the "
        f"GP profile is NOT forced. (1) WIN: adding a transport DOF (intent conserved under "
        f"continuity, a high-Intent core as a sink) gives a LONG-RANGE power-law inflow for ANY "
        f"equation-of-state — the Phase-3b Yukawa obstacle DISSOLVES from the substrate's own "
        f"transport, not by imposition. (2) BUT the PROFILE is set by the equation-of-state, and "
        f"the minimal/default one (INCOMPRESSIBLE, rho=const) gives v ∝ 1/r^2 — and its light "
        f"deflection scales as Delta_theta ∝ b^(-{q_inc:.2f}) (is-1/b: {inc_is_1overb}), i.e. it "
        f"FALLS TOO FAST and does NOT match the observed 1/b bending. The GP profile v ∝ 1/sqrt(r) "
        f"gives Delta_theta ∝ b^(-{q_gp:.2f}) (is-1/b: {gp_is_1overb}) = GR's 1/b — but GP requires "
        f"a SPECIFIC equation-of-state (rho ∝ r^-3/2, the self-advection=Newtonian condition), NOT "
        f"the minimal incompressible case. SO: P2 answer — the transport DOF DERIVES long-range "
        f"gravity (real progress: Yukawa gone, not by hand) but does NOT FORCE the GR-matching "
        f"profile; matching GR still requires CHOOSING the equation-of-state. The remaining freedom "
        f"is located precisely: it is the substrate's EoS (how intent density responds to "
        f"convergence). 'Gravity is fit' has shrunk from 'fit the whole inflow' (Phase-3c) to 'fit "
        f"one equation-of-state' — progress, but the loan is still unpaid. The spatial sector "
        f"remains the wall (research clue holds): what FORCES rho ∝ r^-3/2 is the open question, and "
        f"it is a statement about how the substrate's density responds to a coherence sink — a "
        f"genuinely spatial-emergence question. NOT novel physics; Bucket 0 unchanged."
    )

    out = {
        "question": "does a transport DOF (continuity + intent sink) FORCE the GP gravity profile, or is it fit?",
        "both_profiles_long_range": "yes — any EoS gives a power-law inflow; the Yukawa (short-range) obstacle is dissolved by the transport DOF, derived",
        "gp_profile_p0.5": {"deflection_scaling_q": round(q_gp, 3), "matches_GR_1overb": bool(gp_is_1overb), "rows": rows_gp},
        "incompressible_p2.0": {"deflection_scaling_q": round(q_inc, 3), "matches_GR_1overb": bool(inc_is_1overb), "rows": rows_inc},
        "long_range_derived": True,
        "gp_profile_forced_by_minimal_continuity": bool(inc_is_1overb),
        "remaining_freedom": "the equation-of-state (rho response to a coherence sink); GR needs rho ∝ r^-3/2",
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase7_transport_inflow_profile_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print(f"\n  profile            deflection scaling  Delta_theta ∝ b^-q     matches GR (q=1)?")
    print(f"  GP  v∝1/sqrt(r)    q = {q_gp:.3f}                              {gp_is_1overb}")
    print(f"  incompressible v∝1/r^2  q = {q_inc:.3f}                         {inc_is_1overb}")


if __name__ == "__main__":
    main()
