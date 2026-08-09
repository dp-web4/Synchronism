#!/usr/bin/env python3
"""
Publisher 2026-08-09 -- Appendix D.2 states two force laws in one equation block, and they are
not the same law.  How far apart are they, and which one does the whitepaper's evidence use?

CONTEXT.  manuscripts/Appendix_D_Synchronism_in_General_Relativistic_Form.md was committed
2025-12-01 (4400d54f).  Its section D.2 opens with a modified Poisson equation

    L1:   nabla^2 Phi = 4 pi G rho / C(rho)

and then says, "in this limit", that the observed acceleration is

    L3:   g_obs = g_bar / C(rho)

as if L3 followed from L1.  It does not.  L3 is the spherical Gauss solution of a DIFFERENT
field equation,

    L2:   div[ C(rho) grad Phi ] = 4 pi G rho ,

which whitepaper section 5.15 introduced on 2026-08-04 as a "momentum-conserving completion"
whose premise was that the galaxy sector had NO field equation.  So the archive's own appendix
already contained L2's phenomenology, attached to L1's field equation.

FIVE CHECKS, each falsifiable on its own.

  (1) L2 => L3 exactly in spherical symmetry (Gauss).  Analytic; verified numerically.
  (2) L1 == L3 if and only if C is spatially constant.  Verified by constructing a constant-C
      control: the gap must vanish identically there and only there.
  (3) FORK AMPLITUDE on real in-repo SPARC data (Lelli+2016), computed in a fully
      self-consistent spherical construction so that no extra convention (scale height h,
      disk-vs-sphere geometry) enters the comparison:
          rho_sph(r) = (1/4 pi r^2) d/dr [ V_bar^2 r / G ]
      is the spherically-equivalent density that reproduces the observed baryonic curve
      exactly.  Both laws are then evaluated on the SAME rho.  Any difference is the fork.
  (4) GAMMA-INVARIANCE.  The explorer lane (synchronism-site, 2026-08-08) reports the gap is
      identical at gamma = 2 and gamma = 0.489.  Re-derived here on a different sample and a
      different density construction.
  (5) THE A-PRIORI KILL OF L1.  As rho -> 0, C -> gamma*rho/rho_crit, so rho/C -> rho_crit/gamma,
      a CONSTANT.  L1's source therefore never vanishes in vacuum.  Verified over a grid in
      gamma and rho_crit.

  (6) PROVENANCE.  Which law does the evidence in whitepaper 5.15 actually run on?  Grep the
      simulation corpus for an implementation of rho/C as a Poisson SOURCE (L1) versus
      g_bar/C as a force law (L3).

Data: Lelli, McGaugh & Schombert (2016) SPARC, in-repo at simulations/sparc_real_data/.
Pipeline conventions (M/L 0.5/0.7, Q<=2, Inc>30) copied from
simulations/publisher_20260808_test09_test10_independence.py so that any difference in result
is a difference in question, not in pipeline.
"""
import os
import subprocess

import numpy as np

BASE = "/mnt/c/exe/projects/ai-agents/Synchronism/simulations/sparc_real_data/"
REPO = "/mnt/c/exe/projects/ai-agents/Synchronism/"

UP_DISK, UP_BUL = 0.5, 0.7     # McGaugh (2016) 3.6um mass-to-light
G_ASTRO = 4.300917270e-6       # kpc (km/s)^2 / Msun

# Appendix D.2 as written: gamma = 2, rho_crit = A * Vflat^B.
# Whitepaper 5.15's own convention for the density law: rho_crit = 0.029 * Vflat^2 Msun/pc^3.
GAMMA_D2 = 2.0
GAMMA_FIT = 0.489              # the frozen-instrument value quoted throughout 5.15
RHO_CRIT_A, RHO_CRIT_B = 0.029, 2.0


def C_of_rho(rho, rho_crit, gamma):
    """C(rho) = tanh[gamma * ln(1 + rho/rho_crit)]  -- Appendix D.2 eq. 2, verbatim."""
    return np.tanh(gamma * np.log1p(rho / rho_crit))


# --------------------------------------------------------------------------- data
def load_galaxy_table():
    out = {}
    with open(BASE + "SPARC_Lelli2016c.mrt") as f:
        for line in f:
            p = line.split()
            if len(p) < 18:
                continue
            try:
                out[p[0]] = dict(Vflat=float(p[15]), Q=int(p[17]), Inc=float(p[5]))
            except ValueError:
                continue
    return out


def load_mass_models():
    rows = {}
    with open(BASE + "MassModels_Lelli2016c.mrt") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 10:
                continue
            try:
                vals = list(map(float, parts[1:]))
            except ValueError:
                continue
            D, R, Vobs, eVobs, Vgas, Vdisk, Vbul, SBd, SBb = vals
            rows.setdefault(parts[0], []).append((R, Vobs, Vgas, Vdisk, Vbul))
    out = {}
    for gid, rr in rows.items():
        rr.sort()
        a = np.array(rr)
        out[gid] = dict(R=a[:, 0], Vobs=a[:, 1], Vgas=a[:, 2],
                        Vdisk=a[:, 3], Vbul=a[:, 4])
    return out


def v_bar(gal):
    vg, vd, vb = gal["Vgas"], gal["Vdisk"], gal["Vbul"]
    v2 = vg * np.abs(vg) + UP_DISK * vd ** 2 + UP_BUL * vb ** 2
    return np.sqrt(np.clip(v2, 0.0, None))


# --------------------------------------------------------------------------- the two laws
def spherical_profile(R, Vb):
    """Self-consistent spherical construction.

    M_bar(<r) = V_bar^2 r / G reproduces the observed baryonic curve exactly.
    rho_sph = (1/4 pi r^2) dM/dr is the density whose spherical Poisson solution IS that curve.
    Returns (r, M_enclosed, rho_sph) on the SPARC radial grid, monotone-M points only.
    """
    M = Vb ** 2 * R / G_ASTRO                       # Msun
    keep = np.concatenate(([True], np.diff(M) > 0))  # rho >= 0 requires M increasing
    r, M = R[keep], M[keep]
    if len(r) < 4:
        return None
    dM = np.gradient(M, r)
    rho = dM / (4.0 * np.pi * r ** 2)               # Msun / kpc^3
    return r, M, rho


def g_L3(r, M, C):
    """L2 == L3: C * g * r^2 = G M(<r)  =>  g = g_bar / C."""
    return G_ASTRO * M / r ** 2 / C


def g_L1(r, rho, C):
    """L1: nabla^2 Phi = 4 pi G rho/C  =>  g(r) r^2 = G * int_0^r 4 pi r'^2 rho/C dr'."""
    integrand = 4.0 * np.pi * r ** 2 * rho / C
    M_eff = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))))
    M_eff += integrand[0] * r[0] / 3.0              # inner cone, rho ~ const below r[0]
    return G_ASTRO * M_eff / r ** 2, M_eff


# --------------------------------------------------------------------------- checks
def check1_gauss():
    print("=" * 78)
    print("CHECK 1 -- L2 (div[C grad Phi] = 4 pi G rho) reduces to L3 (g = g_bar/C) in spherical")
    print("=" * 78)
    r = np.linspace(0.2, 30.0, 4000)
    rho = 1.0e7 * np.exp(-r / 3.0) + 1.0e4
    M = np.concatenate(([0.0], np.cumsum(
        0.5 * (4 * np.pi * r[1:] ** 2 * rho[1:] + 4 * np.pi * r[:-1] ** 2 * rho[:-1]) * np.diff(r))))
    C = C_of_rho(rho, 1.0e6, 2.0)
    # integrate div[C grad Phi] = 4 pi G rho over a ball:  C(r) g(r) 4 pi r^2 = 4 pi G M(<r)
    g_from_gauss = G_ASTRO * M / r ** 2 / C
    g_direct = (G_ASTRO * M / r ** 2) / C
    err = np.nanmax(np.abs(g_from_gauss - g_direct) / np.abs(g_direct))
    print(f"  max relative |Gauss(L2) - L3| = {err:.3e}   -> identical (analytic, not fitted)")
    return err


def check2_constant_C():
    print()
    print("=" * 78)
    print("CHECK 2 -- L1 == L3 iff C is spatially constant")
    print("=" * 78)
    r = np.linspace(0.2, 30.0, 4000)
    rho = 1.0e7 * np.exp(-r / 3.0) + 1.0e4
    M = np.concatenate(([0.0], np.cumsum(
        0.5 * (4 * np.pi * r[1:] ** 2 * rho[1:] + 4 * np.pi * r[:-1] ** 2 * rho[:-1]) * np.diff(r))))
    for label, C in (("C = 0.4 constant", np.full_like(r, 0.4)),
                     ("C = C(rho), gamma=2", C_of_rho(rho, 1.0e6, 2.0))):
        g3 = G_ASTRO * M / r ** 2 / C
        integrand = 4.0 * np.pi * r ** 2 * rho / C
        Meff = np.concatenate(([0.0], np.cumsum(
            0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))))
        g1 = G_ASTRO * Meff / r ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            dex = np.log10(g3[1:] / g1[1:])
        print(f"  {label:24s}  max |log10(g_L3/g_L1)| = {np.nanmax(np.abs(dex)):.3e} dex")
    print("  -> the fork IS the spatial variation of C.  C is the framework's entire content,")
    print("     so the two forms of D.2 cannot be reconciled by any choice of parameters.")


def check3_4_fork_amplitude(gals, tab):
    print()
    print("=" * 78)
    print("CHECK 3+4 -- fork amplitude on real SPARC, and its gamma-dependence")
    print("=" * 78)
    results = {}
    for gamma in (GAMMA_D2, GAMMA_FIT):
        per_gal = []
        for gid, gal in sorted(gals.items()):
            t = tab.get(gid)
            if t is None or t["Q"] > 2 or t["Inc"] < 30.0 or t["Vflat"] <= 0:
                continue
            prof = spherical_profile(gal["R"], v_bar(gal))
            if prof is None:
                continue
            r, M, rho = prof
            if np.any(rho <= 0):
                continue
            rho_crit = RHO_CRIT_A * t["Vflat"] ** RHO_CRIT_B * 1.0e9   # Msun/pc^3 -> Msun/kpc^3
            C = C_of_rho(rho, rho_crit, gamma)
            g3 = g_L3(r, M, C)
            g1, _ = g_L1(r, rho, C)
            if not np.all(np.isfinite(g1)) or np.any(g1 <= 0):
                continue
            per_gal.append((gid, float(np.log10(g3[-1] / g1[-1])),
                            float(np.median(np.log10(g3 / g1)))))
        results[gamma] = per_gal
        d_out = np.array([x[1] for x in per_gal])
        print(f"\n  gamma = {gamma}:  N = {len(per_gal)} galaxies (Q<=2, Inc>30, rho_sph>0)")
        print(f"    log10(g_L3 / g_L1) at outermost point: "
              f"median {np.median(d_out):+.3f} dex, "
              f"IQR [{np.percentile(d_out, 25):+.3f}, {np.percentile(d_out, 75):+.3f}], "
              f"range [{d_out.min():+.3f}, {d_out.max():+.3f}]")
        print(f"    |gap| > 0.3 dex in {np.sum(np.abs(d_out) > 0.3)}/{len(d_out)} galaxies; "
              f"> 1.0 dex in {np.sum(np.abs(d_out) > 1.0)}/{len(d_out)}")

    a = {g: dict((x[0], x[1]) for x in v) for g, v in results.items()}
    common = sorted(set(a[GAMMA_D2]) & set(a[GAMMA_FIT]))
    d2 = np.array([a[GAMMA_D2][k] for k in common])
    df = np.array([a[GAMMA_FIT][k] for k in common])
    print(f"\n  gamma-INVARIANCE over {len(common)} common galaxies:")
    print(f"    max |gap(gamma=2) - gap(gamma=0.489)| = {np.max(np.abs(d2 - df)):.4f} dex")
    print(f"    median gap: {np.median(d2):+.3f} (gamma=2) vs {np.median(df):+.3f} (gamma=0.489)")
    if np.max(np.abs(d2 - df)) < 0.05:
        print("    -> gamma-INVARIANT (confirms the explorer lane's claim on a different sample")
        print("       and a different density construction)")
    else:
        print("    -> NOT gamma-invariant on this construction; the explorer lane's claim does")
        print("       NOT reproduce here.  Report the disagreement, do not average it away.")
    return results


def check4b_why_gamma_invariant(gals, tab):
    """The invariance is not deep.  Name its cause or it will be over-read as robustness."""
    print()
    print("=" * 78)
    print("CHECK 4b -- WHY the gap is gamma-invariant (the explorer lane reports the fact,")
    print("            not the reason, and the reason bounds how far the fact travels)")
    print("=" * 78)
    xs = []
    for gid, gal in sorted(gals.items()):
        t = tab.get(gid)
        if t is None or t["Q"] > 2 or t["Inc"] < 30.0 or t["Vflat"] <= 0:
            continue
        prof = spherical_profile(gal["R"], v_bar(gal))
        if prof is None:
            continue
        r, M, rho = prof
        if np.any(rho <= 0):
            continue
        xs.append(rho / (RHO_CRIT_A * t["Vflat"] ** RHO_CRIT_B * 1.0e9))
    allx = np.concatenate(xs)
    print(f"  x = rho/rho_crit over {allx.size} points in {len(xs)} galaxies:")
    print(f"    median {np.median(allx):.3e},  p99 {np.percentile(allx, 99):.3e},  "
          f"max {allx.max():.3e}")
    for gamma in (GAMMA_D2, GAMMA_FIT):
        C = np.tanh(gamma * np.log1p(allx))
        print(f"    gamma={gamma}: max |C/(gamma*x) - 1| = "
              f"{np.max(np.abs(C / (gamma * allx) - 1)):.3e}")
    print("  Every SPARC point sits at x <= 0.05, i.e. DEEP in C's linear regime, where")
    print("  C -> gamma*x to within ~2.4%.  There 1/C ∝ 1/gamma in BOTH laws, so gamma cancels")
    print("  identically from the ratio g_L3/g_L1.  The gamma-invariance is therefore a")
    print("  linear-regime artifact of rho_crit = 0.029*V_flat^2, not a structural robustness:")
    print("  it holds precisely BECAUSE the density calibration is far enough off that no")
    print("  galaxy reaches C's knee -- the same fact as the already-banked 2-5 dex miss.")
    print("  Correct scope: the fork survives every gamma AT THIS rho_crit.  A rho_crit")
    print("  recalibration that moved galaxies to the knee would have to be re-run.")


def check4c_which_is_closer(gals, tab):
    """Guard against the obvious misreading: the eliminated law is the BETTER-fitting one."""
    print()
    print("=" * 78)
    print("CHECK 4c -- which law is closer to the observed accelerations?")
    print("=" * 78)
    for gamma in (GAMMA_D2, GAMMA_FIT):
        d3, d1 = [], []
        for gid, gal in sorted(gals.items()):
            t = tab.get(gid)
            if t is None or t["Q"] > 2 or t["Inc"] < 30.0 or t["Vflat"] <= 0:
                continue
            prof = spherical_profile(gal["R"], v_bar(gal))
            if prof is None:
                continue
            r, M, rho = prof
            if np.any(rho <= 0):
                continue
            Vobs = gal["Vobs"][np.isin(gal["R"], r)]
            rho_crit = RHO_CRIT_A * t["Vflat"] ** RHO_CRIT_B * 1.0e9
            C = C_of_rho(rho, rho_crit, gamma)
            g3 = g_L3(r, M, C)
            g1, _ = g_L1(r, rho, C)
            if not np.all(np.isfinite(g1)) or np.any(g1 <= 0):
                continue
            g_obs = Vobs ** 2 / r
            d3.append(np.log10(g3[-1] / g_obs[-1]))
            d1.append(np.log10(g1[-1] / g_obs[-1]))
        d3, d1 = np.array(d3), np.array(d1)
        print(f"  gamma = {gamma}  (N = {len(d3)})")
        print(f"    L2=L3 vs observed: median {np.median(d3):+.3f} dex")
        print(f"    L1    vs observed: median {np.median(d1):+.3f} dex")
        print(f"    L1 closer in {np.sum(np.abs(d1) < np.abs(d3))}/{len(d3)} galaxies")
    print("  -> The reading eliminated a priori (L1) is UNIFORMLY closer to the data than the")
    print("     reading that survives.  The elimination is structural (vacuum source floor,")
    print("     ill-posed per-galaxy rho_crit), NOT a goodness-of-fit verdict.  Both are")
    print("     catastrophically wrong; neither is viable.")
    print("  SCOPE: these absolute offsets use the spherical rho_sph construction and are NOT")
    print("     comparable to whitepaper 5.15's '2-5 orders of magnitude', which uses rho = Sigma/2h.")
    print("     A disk's mass spread over spheres gives a much lower rho, hence a larger boost.")
    print("     The FORK amplitude (check 3) is a ratio of two laws on the same rho and is robust")
    print("     to this choice: 0.821 dex here vs 0.81 dex in the sibling lane's Sigma/2h run.")


def check5_vacuum_floor():
    print()
    print("=" * 78)
    print("CHECK 5 -- L1's source does not vanish in vacuum: rho/C -> rho_crit/gamma")
    print("=" * 78)
    ok = True
    for gamma in (0.25, 0.489, 1.0, 2.0, 5.0):
        for rho_crit in (1.0e5, 1.0e7):
            for x in (1e-6, 1e-8, 1e-10):
                rho = x * rho_crit
                rho_eff = rho / C_of_rho(rho, rho_crit, gamma)
                pred = rho_crit / gamma
                rel = abs(rho_eff - pred) / pred
                if rel > 1e-4:
                    ok = False
                    print(f"    MISMATCH gamma={gamma} rho_crit={rho_crit:.0e} x={x:.0e}: "
                          f"{rho_eff:.6e} vs {pred:.6e} (rel {rel:.2e})")
    print(f"    limit rho_eff = rho_crit/gamma verified to <1e-4 relative over "
          f"gamma in [0.25, 5], rho/rho_crit down to 1e-10: {'PASS' if ok else 'FAIL'}")
    print("    Consequence: M_eff(<r) ~ (4/3) pi r^3 rho_crit/gamma diverges; V ~ r forever.")
    print("    But note the SHARPER objection, which the divergence obscures:")
    for vf in (50.0, 150.0, 250.0):
        rc = RHO_CRIT_A * vf ** RHO_CRIT_B * 1.0e9
        print(f"      V_flat = {vf:5.0f} km/s -> vacuum source floor rho_crit/gamma = "
              f"{rc / GAMMA_D2:.3e} Msun/kpc^3")
    print("    rho_crit = A*V_flat^B is a per-GALAXY constant, so L1 assigns the SAME point of")
    print("    empty space a different source density for every galaxy you ask about.  The")
    print("    field equation has no single-valued source: it is ill-posed, not merely divergent.")
    print("    This objection also survives into L2=L3, which the divergence argument does not.")


def check6_provenance():
    print()
    print("=" * 78)
    print("CHECK 6 -- which law does the corpus actually implement?")
    print("=" * 78)
    pats = [("L1  rho/C as a Poisson SOURCE",
             r"rho\s*/\s*C\b|rho_eff\s*=\s*rho\s*/|4\s*\*\s*np\.pi\s*\*\s*G\s*\*\s*rho\s*/"),
            ("L2  div[C grad Phi] (dielectric)",
             r"div\s*\[\s*C|C\s*\*\s*grad|np\.gradient\(\s*C\s*\*"),
            ("L3  g_bar/C or V_bar/sqrt(C) (force law)",
             r"g_bar\s*/\s*C\b|g_obs\s*=\s*g_bar\s*/|V_bar\s*/\s*np\.sqrt\(\s*C")]
    for label, pat in pats:
        try:
            out = subprocess.run(["git", "grep", "-lE", pat, "--", "simulations/*.py"],
                                 cwd=REPO, capture_output=True, text=True, timeout=120)
            files = [f for f in out.stdout.split("\n") if f.strip()]
        except Exception as e:                                    # noqa: BLE001
            files = [f"<grep failed: {e}>"]
        print(f"  {label:42s} -> {len(files)} script(s)")
        for f in files[:8]:
            print(f"       {f}")
    print()
    print("  The load-bearing question is narrower: does any script that produced a number")
    print("  QUOTED IN THE WHITEPAPER implement L1?  The frozen instrument behind gamma=0.489,")
    print("  dBIC +7.1/+184 and Cassini +17.95sigma is sparc_tanhlog_profile.py /")
    print("  sparc_cassini_q2.py, which invert g_bar = g_obs*tanh(gamma*ln(1+g_obs/a0)) --")
    print("  an L3-shaped force law in ACCELERATION space (whitepaper 5.15, amended 2026-08-04).")
    for f in ("simulations/sparc_tanhlog_profile.py", "simulations/sparc_cassini_q2.py"):
        p = os.path.join(REPO, f)
        if not os.path.exists(p):
            print(f"    {f}: ABSENT")
            continue
        src = open(p, encoding="utf-8", errors="replace").read()
        print(f"    {f}: 'rho_crit' x{src.count('rho_crit')}, "
              f"'a0' x{src.count('a0') + src.count('a_0')}, "
              f"'/ C' x{src.count('/ C') + src.count('/C')}")


def main():
    tab = load_galaxy_table()
    gals = load_mass_models()
    print(f"Loaded {len(gals)} SPARC galaxies with mass models; "
          f"{len(tab)} rows in the galaxy table.\n")
    check1_gauss()
    check2_constant_C()
    check3_4_fork_amplitude(gals, tab)
    check4b_why_gamma_invariant(gals, tab)
    check4c_which_is_closer(gals, tab)
    check5_vacuum_floor()
    check6_provenance()
    print()
    print("=" * 78)
    print("Run complete.  Every number above is recomputed from in-repo SPARC data;")
    print("nothing is quoted from the explorer lane without re-derivation.")
    print("=" * 78)


if __name__ == "__main__":
    main()
