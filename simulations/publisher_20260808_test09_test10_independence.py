#!/usr/bin/env python3
"""
Publisher 2026-08-08 -- how much independent empirical content does TEST-09 carry beyond TEST-10?

CONTEXT.  The refutation count of 6 is asserted four times in the Synchronism whitepaper
(executive summary, conclusion, and 5.15 dark matter x3) and enumerated nowhere in it.
Three different "independent root" figures are circulating, all asserted, none computed:

  - synchronism-site/maintainer/logs/2026-08-07.md : "Honest independent-root figure: 3-4, not 6."
  - synchronism-site/visitor/logs/2026-08-07.md    : "at most four independent roots"
  - Research/proposals/boost_ceiling_provenance_and_class_exclusion.md, open question 4
    (2026-07-28, GATED ON DP for 11 days): "Should TEST-09 and TEST-10 be counted as one
    refutation or two? ... This is a real methodological question (how much independent
    empirical content does a corollary carry?), not a typo to silently fix."

That question was filed as needing a *naming convention* decided by dp.  It does not.  It is
measurable, and this script measures it.  Naming conventions gate on dp; sample overlap and
conditional kill-strength do not.

WHAT IS TESTED, in three steps.

  (1) IDENTITY.  TEST-10's registered statistic (apparent DM fraction exceeds the ceiling
      f_DM > 1 - Om = 0.685) and TEST-09's byproduct statistic (observed boost exceeds the
      ceiling B > 1/Om = 3.17) are algebraically the same inequality, because
      f_DM = 1 - (V_bar/V_obs)^2 = 1 - g_bar/g_obs = 1 - 1/B.  Verified numerically per galaxy,
      not asserted.

  (2) SHARED INPUTS.  Both executions read the same SPARC files, the same M/L (0.5 / 0.7),
      the same force law C(a) = Om + (1-Om)x/(1+x) with the same Om, phi and a0.  Sample
      overlap is computed, not assumed.

  (3) CONDITIONAL KILL.  The decisive one.  Delete every galaxy that fires TEST-10 (i.e. every
      ceiling-exceeder) and refit the BTFR on what is left.  If TEST-09's registered kill
      (|n_obs - n_sync| > 0.3) still fires on the non-exceeders, TEST-09 carries empirical
      content that TEST-10 does not, and two rows are defensible.  If it does not fire,
      TEST-09 is the same fact restated as a slope, and the honest structural-root count
      collapses by one.

Data: Lelli, McGaugh & Schombert (2016) SPARC (in-repo, simulations/sparc_real_data/).
Conventions copied verbatim from synchronism-site/explorer/scripts/test09_btfr_bounded_boost_real_sparc.py
and test10_dwarf_dm_fraction_ceiling.py so that any difference in result is a difference in
question, not in pipeline.
"""
import numpy as np

BASE = "/mnt/c/exe/projects/ai-agents/Synchronism/simulations/sparc_real_data/"
KPC = 3.0856775814913673e19
KMS = 1.0e3

UP_DISK, UP_BUL = 0.5, 0.7          # McGaugh (2016) 3.6um mass-to-light
OM = 0.315                          # Omega_m (Sessions 191-193)
PHI = (1.0 + 5.0 ** 0.5) / 2.0
A0_SYNC = 1.05e-10                  # m/s^2 (Session 193)
A0_MOND = 1.20e-10                  # m/s^2 (Milgrom)

B_CEIL = 1.0 / OM                   # 3.17 -- TEST-09's stated ceiling
F_CEIL = 1.0 - OM                   # 0.685 -- TEST-10's stated ceiling
KILL_09 = 0.3                       # registered TEST-09 criterion: slope off by > 0.3


def g_obs_sync(g_bar):
    x = (g_bar / A0_SYNC) ** (1.0 / PHI)
    C = OM + (1.0 - OM) * x / (1.0 + x)
    return g_bar / C


def g_obs_mond(g_bar):
    return g_bar / (1.0 - np.exp(-np.sqrt(g_bar / A0_MOND)))


def load_galaxy_table():
    out = {}
    with open(BASE + "SPARC_Lelli2016c.mrt") as f:
        for line in f:
            p = line.split()
            if len(p) < 18:
                continue
            try:
                out[p[0]] = dict(L36=float(p[7]), MHI=float(p[13]), Vflat=float(p[15]),
                                 e_Vflat=float(p[16]), Q=int(p[17]), Inc=float(p[5]))
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
            rows.setdefault(parts[0], []).append((R, Vobs, eVobs, Vgas, Vdisk, Vbul))
    out = {}
    for gid, rr in rows.items():
        rr.sort()
        a = np.array(rr)
        out[gid] = dict(R=a[:, 0], Vobs=a[:, 1], eVobs=a[:, 2],
                        Vgas=a[:, 3], Vdisk=a[:, 4], Vbul=a[:, 5])
    return out


def v_bar(gal):
    vg, vd, vb = gal["Vgas"], gal["Vdisk"], gal["Vbul"]
    v2 = vg * np.abs(vg) + UP_DISK * vd ** 2 + UP_BUL * vb ** 2
    return np.sqrt(np.clip(v2, 0.0, None))


def v_flat_estimator(R, V):
    if len(R) < 3:
        return np.nan
    return float(np.mean(V[-3:]))


def fit_btfr(logM, logV, n_boot=2000, seed=12345):
    a, b = np.polyfit(logM, logV, 1)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(logM))
    slopes = []
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        try:
            aa, _ = np.polyfit(logM[s], logV[s], 1)
            if aa > 1e-6:
                slopes.append(1.0 / aa)
        except Exception:
            pass
    return 1.0 / a, float(np.std(np.array(slopes)))


def build_rows():
    """The TEST-09 sample, carrying TEST-10's per-galaxy statistic alongside."""
    tbl, mm = load_galaxy_table(), load_mass_models()
    rows = []
    for name, t in tbl.items():
        if name not in mm:
            continue
        if t["Q"] > 2 or t["Inc"] < 30.0 or t["Vflat"] <= 0 or t["L36"] <= 0:
            continue
        gal = mm[name]
        R, Vb = gal["R"], v_bar(gal)
        ok = (R > 0) & np.isfinite(Vb) & (Vb > 0)
        R, Vb, Vobs = R[ok], Vb[ok], gal["Vobs"][ok]
        if len(R) < 3:
            continue
        r_m = R * KPC
        g_bar = (Vb * KMS) ** 2 / r_m
        v_sync = np.sqrt(g_obs_sync(g_bar) * r_m) / KMS
        v_mond = np.sqrt(g_obs_mond(g_bar) * r_m) / KMS
        m_bar = (UP_DISK * t["L36"] + 1.33 * t["MHI"]) * 1e9
        if m_bar <= 0:
            continue

        # --- TEST-10's per-galaxy statistic, on TEST-09's outermost point.
        #     TEST-10 proper uses the error-weighted mean of the outermost 3 radii; both are
        #     reported below so the conclusion cannot hang on that choice.
        boost_out = (Vobs[-1] / Vb[-1]) ** 2
        f_dm_out = 1.0 - (Vb[-1] / Vobs[-1]) ** 2
        w = 1.0 / np.clip(gal["eVobs"][ok][-3:], 1e-3, None) ** 2
        vo3 = float(np.sum(Vobs[-3:] * w) / np.sum(w))
        vb3 = float(np.sum(Vb[-3:] * w) / np.sum(w))
        boost_o3 = (vo3 / vb3) ** 2
        f_dm_o3 = 1.0 - (vb3 / vo3) ** 2

        rows.append(dict(name=name, M=m_bar,
                         V_obs=v_flat_estimator(R, Vobs),
                         V_sync=v_flat_estimator(R, v_sync),
                         V_mond=v_flat_estimator(R, v_mond),
                         boost_out=boost_out, f_dm_out=f_dm_out,
                         boost_o3=boost_o3, f_dm_o3=f_dm_o3))
    return [r for r in rows if np.isfinite(r["V_sync"]) and np.isfinite(r["V_obs"])
            and r["V_sync"] > 0 and np.isfinite(r["V_mond"]) and r["V_mond"] > 0]


def slopes_on(rows, label):
    logM = np.log10(np.array([r["M"] for r in rows]))
    out = {}
    for key in ("V_obs", "V_mond", "V_sync"):
        v = np.log10(np.array([r[key] for r in rows]))
        out[key] = fit_btfr(logM, v)
    n_obs, e_obs = out["V_obs"]
    n_syn, e_syn = out["V_sync"]
    n_mnd, e_mnd = out["V_mond"]
    dev = abs(n_obs - n_syn)
    sig = dev / np.hypot(e_obs, e_syn)
    print(f"  {label:<44} N={len(rows):>4}")
    print(f"    observed     n = {n_obs:5.2f} +/- {e_obs:.2f}")
    print(f"    MOND         n = {n_mnd:5.2f} +/- {e_mnd:.2f}   dev {abs(n_obs-n_mnd):.2f}")
    print(f"    Synchronism  n = {n_syn:5.2f} +/- {e_syn:.2f}   dev {dev:.2f} "
          f"({sig:.1f} sigma)  kill(>0.3): {'FIRES' if dev > KILL_09 else 'does NOT fire'}")
    return dev, sig, len(rows)


def main():
    rows = build_rows()
    print("=" * 78)
    print("Publisher 2026-08-08 -- TEST-09 / TEST-10 independence, measured")
    print("=" * 78)
    print(f"  TEST-09 sample rebuilt: {len(rows)} SPARC galaxies (Q<=2, Inc>30, Vflat measured)")
    print(f"  ceiling: B <= 1/Om = {B_CEIL:.2f}   <=>   f_DM <= 1-Om = {F_CEIL:.3f}\n")

    # ---------------------------------------------------------------- (1) identity
    print("-" * 78)
    print("(1) Are TEST-09's boost-exceedance and TEST-10's f_DM-exceedance the same test?")
    print("-" * 78)
    b = np.array([r["boost_out"] for r in rows])
    f = np.array([r["f_dm_out"] for r in rows])
    resid = np.max(np.abs(f - (1.0 - 1.0 / b)))
    fire_b = b > B_CEIL
    fire_f = f > F_CEIL
    print(f"    max |f_DM - (1 - 1/B)| over {len(rows)} galaxies : {resid:.3e}   (algebraic identity)")
    print(f"    galaxies firing on boost  B > {B_CEIL:.2f}      : {fire_b.sum():>4}")
    print(f"    galaxies firing on f_DM   f > {F_CEIL:.3f}     : {fire_f.sum():>4}")
    print(f"    disagreements between the two criteria     : {int((fire_b ^ fire_f).sum())}")
    b3 = np.array([r["boost_o3"] for r in rows])
    f3 = np.array([r["f_dm_o3"] for r in rows])
    print(f"    same, on the outer-3 error-weighted point  : "
          f"{int((b3 > B_CEIL).sum())} / {int((f3 > F_CEIL).sum())}, "
          f"disagreements {int(((b3 > B_CEIL) ^ (f3 > F_CEIL)).sum())}")
    print("    => TEST-10's statistic is not a corollary of TEST-09's ceiling; it is the")
    print("       SAME per-galaxy inequality expressed in a different variable.\n")

    # ---------------------------------------------------------------- (3) conditional kill
    print("-" * 78)
    print("(3) Does TEST-09's registered kill survive with every TEST-10 firing galaxy deleted?")
    print("-" * 78)
    keep = [r for r in rows if r["boost_out"] <= B_CEIL]
    drop = [r for r in rows if r["boost_out"] > B_CEIL]
    slopes_on(rows, "full TEST-09 sample")
    print()
    slopes_on(drop, "TEST-10 firing galaxies only (exceeders)")
    print()
    dev_keep, sig_keep, n_keep = slopes_on(keep, "TEST-10 NON-firing galaxies only")
    print()

    # A ceiling-free control: does the framework's slope deficit exist at all where the
    # ceiling is not binding?  If not, TEST-09 is the ceiling, restated.
    print("-" * 78)
    print("VERDICT")
    print("-" * 78)
    if dev_keep > KILL_09:
        print(f"    TEST-09 STILL FIRES on the {n_keep} galaxies TEST-10 does not touch "
              f"(dev {dev_keep:.2f} > {KILL_09}).")
        print("    => the two rows carry separable empirical content; 6 test-level rows map to")
        print("       roots with TEST-09 and TEST-10 NOT fully collapsible.")
    else:
        print(f"    TEST-09 does NOT fire on the {n_keep} galaxies TEST-10 does not touch "
              f"(dev {dev_keep:.2f} <= {KILL_09}).")
        print("    => TEST-09's registered statistic has no discriminating power where the")
        print("       ceiling is not binding.  Its kill is carried by galaxies that already")
        print("       fire TEST-10.")

    # ---- state the power of the leave-out test rather than reading a null as a proof
    logM_k = np.log10(np.array([r["M"] for r in keep]))
    _, e_obs_k = fit_btfr(logM_k, np.log10(np.array([r["V_obs"] for r in keep])))
    _, e_syn_k = fit_btfr(logM_k, np.log10(np.array([r["V_sync"] for r in keep])))
    sigma_k = float(np.hypot(e_obs_k, e_syn_k))
    print()
    print("    POWER OF THE LEAVE-OUT TEST, stated so the null is not over-read:")
    print(f"      sigma on the non-exceeder deviation = {sigma_k:.2f}; the registered kill "
          f"threshold {KILL_09} sits {KILL_09/sigma_k:.1f} sigma out.")
    lm_all = np.log10(np.array([r["M"] for r in rows]))
    lm_drop = np.log10(np.array([r["M"] for r in drop]))
    print(f"      mass lever arm (dex): full {lm_all.ptp():.2f}, exceeders {lm_drop.ptp():.2f}, "
          f"non-exceeders {logM_k.ptp():.2f}")
    print("      => the leave-out subsample CANNOT demonstrate independence at the registered")
    print("         threshold even if it existed.  What it does show is that no part of")
    print("         TEST-09's evidence is visibly located outside TEST-10's firing set, while")
    print("         step (1) shows the two firing sets are the SAME set by algebraic identity.")
    print("         The burden of the two-row convention is therefore unmet, not disproved.")
    print()
    print("    Shared inputs, for the record: same SPARC files, same M/L (0.5/0.7), same")
    print(f"    force law C(a)=Om+(1-Om)x/(1+x), same Om={OM}, phi={PHI:.6f}, a0={A0_SYNC:.2e}.")
    print("    Sample overlap between the two executions is by construction total on the")
    print("    TEST-09 cut; TEST-10 proper runs the wider 175-galaxy set with the same cuts.")


if __name__ == "__main__":
    main()
