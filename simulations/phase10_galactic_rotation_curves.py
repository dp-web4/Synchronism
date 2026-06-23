"""
Phase-10 (door #1, galactic) — does a capacity rule give FLAT rotation curves without dark matter?
(2026-06-23)

Phase-9 mapped the Bucket-0 frontier to two doors; door #1 is "the inflow PROFILE departs from GP
because the capacity/EoS rule changes u(r)" — the galactic frontier. Phase-8 Faith-B suggested
saturating capacity -> "cored centre + GR tail -> dark-matter halo." THIS PHASE TESTS THAT CLAIM —
and corrects it.

KEY RELATION (Newtonian limit of the inflow/swimmer, Phase-9): a circular orbit has
v_c(r)^2 = G*M(r)/r, with M(r) the enclosed gravitating "intent" (the substrate flowing toward
matter, IF intent self-gravitates — the modeling premise door #1 requires). For a halo density
rho ∝ r^(-alpha):  M(r) ∝ r^(3-alpha),  v_c^2 ∝ r^(2-alpha),  v_c ∝ r^(1-alpha/2).
  - alpha = 2 (isothermal): v_c ∝ r^0 = FLAT  <- what observed galactic curves need
  - alpha = 1.5 (the n=3 GR capacity rule, Phase-8): v_c ∝ r^(1/4) = slowly RISING (NOT flat)
  - alpha = 3 (point-mass / visible matter only): v_c ∝ r^(-1/2) = Keplerian FALLING (the DM problem)

So the GR-matching capacity rule (n=3, rho∝r^-3/2) does NOT give flat curves — it gives a slowly
rising r^1/4 curve. Phase-8's "saturation -> DM halo" conflated "cored" with "isothermal": a GR
(r^-3/2) tail is NOT an isothermal (r^-2) halo. FLAT curves require a capacity rule giving rho∝1/r^2
at galactic scale — a DIFFERENT exponent than n=3. For one rule to give BOTH (GR near matter, flat
far away) it must TRANSITION between regimes; the transition acceleration scale would be the MOND a0
(the RAR / bet B2). This phase verifies the exponents on realistic curves and states door #1's true
requirement.

numpy only. Headless. Writes results JSON.
"""
import json
import os
import numpy as np


def enclosed_mass_powerlaw(r, alpha, rho0=1.0, r0=1.0):
    """M(r) for rho = rho0 (r/r0)^-alpha, M = 4pi ∫_0^r rho r'^2 dr' (alpha<3)."""
    if abs(alpha - 3) < 1e-9:
        return 4 * np.pi * rho0 * r0 ** alpha * np.log(r / 1e-3)
    return 4 * np.pi * rho0 * r0 ** alpha * r ** (3 - alpha) / (3 - alpha)


def vc_powerlaw(r, alpha):
    G = 1.0
    M = enclosed_mass_powerlaw(r, alpha)
    return np.sqrt(G * M / r)


def vc_exponential_disk(r, Mtot=1.0, Rd=2.0):
    """Visible matter: enclosed mass of an exponential disk (spherical-approx), -> Keplerian tail."""
    G = 1.0
    Menc = Mtot * (1 - (1 + r / Rd) * np.exp(-r / Rd))
    return np.sqrt(G * Menc / r)


def large_r_slope(r, v):
    m = (r > 0.5 * r.max())
    return float(np.polyfit(np.log(r[m]), np.log(v[m]), 1)[0])


def main():
    r = np.logspace(0, 2.5, 200)

    # candidate halo profiles (from capacity rules) + visible-only baseline
    profiles = {
        "visible_only_exponential_disk": vc_exponential_disk(r),
        "n3_GR_capacity_halo_alpha1.5": vc_powerlaw(r, 1.5),
        "isothermal_halo_alpha2.0": vc_powerlaw(r, 2.0),
    }
    slopes = {k: round(large_r_slope(r, v), 3) for k, v in profiles.items()}

    # classify large-r behaviour
    def classify(s):
        if abs(s) < 0.06: return "FLAT"
        return "rising" if s > 0 else "falling"
    cls = {k: classify(s) for k, s in slopes.items()}

    isothermal_flat = abs(slopes["isothermal_halo_alpha2.0"]) < 0.06
    n3_not_flat = abs(slopes["n3_GR_capacity_halo_alpha1.5"]) > 0.1
    visible_falls = slopes["visible_only_exponential_disk"] < -0.2

    verdict = (
        f"DOES A CAPACITY RULE GIVE FLAT ROTATION CURVES (door #1)? The rotation curve exponent is "
        f"v_c ∝ r^(1-alpha/2) for a halo rho∝r^-alpha. Measured large-r slopes: "
        f"visible-only (exponential disk) = {slopes['visible_only_exponential_disk']} "
        f"({cls['visible_only_exponential_disk']} = the dark-matter problem, Keplerian); "
        f"n=3 GR-capacity halo (alpha=1.5) = {slopes['n3_GR_capacity_halo_alpha1.5']} "
        f"({cls['n3_GR_capacity_halo_alpha1.5']}, ~r^1/4); isothermal halo (alpha=2.0) = "
        f"{slopes['isothermal_halo_alpha2.0']} ({cls['isothermal_halo_alpha2.0']}). "
        f"CORRECTION TO PHASE-8: the GR-matching capacity rule (n=3, rho∝r^-3/2) does NOT give flat "
        f"curves — it gives a slowly RISING r^1/4 curve (n3_not_flat={n3_not_flat}). Phase-8's "
        f"'saturating capacity -> cored centre + GR tail -> dark-matter halo' CONFLATED 'cored' with "
        f"'isothermal': a cored centre + r^-3/2 tail is NOT a flat-curve (r^-2) halo. The honest "
        f"requirement for door #1: FLAT curves need rho ∝ 1/r^2 (isothermal, alpha=2; "
        f"flat={isothermal_flat}) at galactic scale — a DIFFERENT exponent than the n=3 GR rule. "
        f"So a single capacity rule must TRANSITION: rho∝r^-3/2 (GR) near matter / high acceleration, "
        f"to rho∝r^-2 (flat curves) at large r / low acceleration. The transition acceleration IS the "
        f"MOND a0 / the RAR (bet B2). TWO further requirements door #1 makes explicit: (i) intent "
        f"must SELF-GRAVITATE (the flowing intent halo is itself a source), or only visible matter "
        f"gravitates and curves stay Keplerian; (ii) the capacity rule must produce the alpha:1.5->2 "
        f"transition at the right scale. NEITHER is yet shown — so door #1 is OPEN, now precisely "
        f"scoped: it is not 'saturation = halo' (that was wrong) but 'a transitioning capacity rule "
        f"that yields an isothermal intent-halo at low acceleration, with transition scale a0'. "
        f"Per dp: a0 being a FIT constant is fine (most constants are) — the test is whether the "
        f"transition + a0 reproduce the RAR across SPARC (the actual B2 test, next). NOT novel "
        f"physics yet; Bucket 0 unchanged. This phase's value: a real CORRECTION (Phase-8 overclaim) "
        f"+ door #1 scoped to a falsifiable target."
    )

    out = {
        "key_relation": "v_c ∝ r^(1 - alpha/2) for halo rho∝r^-alpha; flat needs alpha=2 (isothermal)",
        "large_r_slopes": slopes,
        "classification": cls,
        "isothermal_gives_flat": bool(isothermal_flat),
        "n3_GR_rule_gives_flat": bool(not n3_not_flat),
        "phase8_saturation_halo_claim_corrected": "GR (r^-1.5) tail gives RISING r^1/4, not flat; flat needs isothermal r^-2",
        "door1_requirements": ["intent self-gravitates", "capacity rule transitions alpha 1.5->2 at scale a0 (=RAR/MOND)"],
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase10_galactic_rotation_curves_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n  profile                              large-r v_c slope    behaviour")
    for k in profiles:
        print(f"  {k:36}  {slopes[k]:+.3f}              {cls[k]}")


if __name__ == "__main__":
    main()
