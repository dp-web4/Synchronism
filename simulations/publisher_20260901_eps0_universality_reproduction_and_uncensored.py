#!/usr/bin/env python3
"""
Publisher 2026-09-01 -- REPRODUCE the site explorer's 2026-08-30 universality test (synchronism-site c66718e,
explorer/findings/scripts/universality_eps0_vs_a0.py) in-repo, then EXTEND it by the first step its own seeded
topic (explorer/topics/eps0-mass-relation-the-last-escape.md) asks for, plus the null it says not to skip.

The site's claim (finding 'the-ceiling-is-a-measurement-and-the-measurement-excludes-it', sec. 7):
  parameter-matched (ONE free constant per galaxy each), under L2 on 153 SPARC discs, Upsilon_d = 0.5 fixed:
    class eps0 : global 126.53 at eps0=0.220 | per-galaxy chi2/N 39.70, wins 21% | scatter 1.197 dex (16-84), 42% at grid edge
    MOND  a0   : global  52.22 at a0=1.202e-10 | per-galaxy 10.30            | scatter 0.619 dex,         2% at grid edge
    Spearman vs log M_bar: rho_s(eps0) = +0.758 (p 7e-30)  vs  rho_s(a0) = +0.073 (p 0.37)
    Spearman vs log rho_mid: +0.162 (p 0.046)  vs  +0.033
  and the seed's step 1: "refit eps0 per galaxy on an UNCENSORED grid (extend to [0.005, 0.98]) so the 42% edge
  pile-up resolves; the current +0.758 is a lower bound BECAUSE of the censoring".
  and its trap: "declare the null by permutation, not from a Spearman table on correlated, censored data".

Solver: the site's committed L2 finite-volume code (l2_sparc_core.py + l2_field_equation_on_sparc.py), imported from
the sibling checkout exactly as publisher_20260829_* did.  Data: in-repo simulations/sparc_real_data/*.mrt (real SPARC).
Nuisances: Upsilon_disk = 0.5 FIXED both sides (as the site's script); distance and inclination NOT marginalised for
either side -- absolute chi2 is mis-scaled and every claim is a RATIO or a RANK between models on identical points.
"""
import os, sys, time, json
import numpy as np
from multiprocessing import get_context
from scipy.stats import spearmanr, wilcoxon

SITE = "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/findings/scripts"
sys.path.insert(0, SITE)
import l2_sparc_core as K                 # noqa: E402
import l2_field_equation_on_sparc as D    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "publisher_20260901_eps0_universality_results.json")
RHO_C_BEST = 3.5e-06
GAMMA = 0.489
EPS_SITE = np.array([0.02, 0.035, 0.05, 0.073, 0.089, 0.12, 0.16, 0.22, 0.315, 0.42, 0.56, 0.661, 0.80])
A0_SITE = np.logspace(np.log10(0.15e-10), np.log10(12e-10), 41)
# extension: uncensored grids (the seed's [0.005, 0.98]); a0 widened symmetrically so the control is not censored either
EPS_WIDE = np.unique(np.concatenate([np.logspace(np.log10(0.005), np.log10(0.98), 28), EPS_SITE]))
A0_WIDE = np.logspace(np.log10(0.02e-10), np.log10(60e-10), 61)

t0 = time.time()
gals = K.load_sparc()
G = [D.Gal(d) for d in gals.values()]
NN = np.array([int(g.ok.sum()) for g in G])
print(f"built {len(G)} galaxies in {time.time()-t0:.0f}s;  N = {NN.sum()} points", flush=True)


def eps_task(e0):
    Cf = K.C_framework(GAMMA, RHO_C_BEST, e0)
    return [g.score(g.solve(Cf, e0)[0] * g.vbar2)[0] for g in G]


def a0_task(a0):
    old = K.A0_KPC
    K.A0_KPC = a0 / 3.24078e-14
    try:
        return [g.score(K.mond_simple(g.gbar_sparc) * g.d["R"])[0] for g in G]
    finally:
        K.A0_KPC = old


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78, flush=True)


ncpu = max(1, min(8, os.cpu_count() or 1))
with get_context("fork").Pool(ncpu) as pool:
    M_eps_w = np.array(pool.map(eps_task, EPS_WIDE))          # (n_eps_wide, n_gal)
M_a0_w = np.array([a0_task(a) for a in A0_WIDE])              # algebraic, cheap
print(f"scored {len(EPS_WIDE)} eps0 x {len(G)} L2 solves + {len(A0_WIDE)} a0 values in {time.time()-t0:.0f}s on {ncpu} cpus", flush=True)

# the site's grids are subsets of the wide ones (EPS_SITE is included exactly; A0_SITE re-scored cheaply)
idx_site = np.array([int(np.argmin(np.abs(EPS_WIDE - e))) for e in EPS_SITE])
assert np.allclose(EPS_WIDE[idx_site], EPS_SITE)
M_eps_s = M_eps_w[idx_site]
M_a0_s = np.array([a0_task(a) for a in A0_SITE])

obs = {
    "log Vflat": np.array([np.log10(max(g.d["props"]["Vflat"], 10.0)) for g in G]),
    "log Rdisk": np.array([np.log10(max(g.d["props"]["Rdisk"], 0.05)) for g in G]),
    "log Mbar ": np.array([np.log10(max(g.Mtot, 1e6)) for g in G]),
    "log rho_m": np.array([np.log10(np.median(g.rho_mid[g.ok] / K.KPC3) + 1e-12) for g in G]),
    "log Breq ": np.array([np.log10(np.max(g.d["Vobs"][g.ok]**2 / g.vbar2[g.ok])) for g in G]),
}
res = {"n_gal": len(G), "N": int(NN.sum())}


def analyse(tag, EPS, M_eps, A0S, M_a0, site_ref=None):
    hdr(f"{tag}: eps0 grid [{EPS.min():.3g}, {EPS.max():.3g}] x{len(EPS)}   a0 grid [{A0S.min():.2e}, {A0S.max():.2e}] x{len(A0S)}")
    ie, ia = int(M_eps.sum(1).argmin()), int(M_a0.sum(1).argmin())
    g_e, g_a = M_eps[ie].sum() / NN.sum(), M_a0[ia].sum() / NN.sum()
    print(f"   GLOBAL  class: eps0 = {EPS[ie]:.3f} (ceiling {1/EPS[ie]:.2f})  chi2/N = {g_e:8.2f}")
    print(f"   GLOBAL  MOND : a0   = {A0S[ia]:.3e}                chi2/N = {g_a:8.2f}   ({A0S[ia]/1.2e-10:.3f}x literature)")
    best_e, arg_e = M_eps.min(0), M_eps.argmin(0)
    best_a, arg_a = M_a0.min(0), M_a0.argmin(0)
    e0g, a0g = EPS[arg_e], A0S[arg_a]
    pg_e, pg_a = best_e.sum() / NN.sum(), best_a.sum() / NN.sum()
    wins = float(np.mean(best_e < best_a))
    pw = wilcoxon(best_e / NN, best_a / NN).pvalue
    print(f"   PER-GAL class: chi2/N = {pg_e:8.2f}  median red = {np.median(best_e/NN):6.2f}")
    print(f"   PER-GAL MOND : chi2/N = {pg_a:8.2f}  median red = {np.median(best_a/NN):6.2f}   class wins {100*wins:.0f}%   Wilcoxon p = {pw:.2e}")
    le, la = np.log10(e0g), np.log10(a0g)
    edge_e = float(np.mean((arg_e == 0) | (arg_e == len(EPS) - 1)))
    edge_a = float(np.mean((arg_a == 0) | (arg_a == len(A0S) - 1)))
    lo_e = float(np.mean(arg_e == 0)); hi_e = float(np.mean(arg_e == len(EPS) - 1))
    sp_e = float(np.percentile(le, 84) - np.percentile(le, 16)); sp_a = float(np.percentile(la, 84) - np.percentile(la, 16))
    print(f"   SCATTER eps0: median {np.median(e0g):.4f}  16-84 {sp_e:.3f} dex  sd {np.std(le):.3f}  at-edge {100*edge_e:.0f}% (low {100*lo_e:.0f}% / high {100*hi_e:.0f}%)")
    print(f"   SCATTER a0  : median {np.median(a0g):.3e}  16-84 {sp_a:.3f} dex  sd {np.std(la):.3f}  at-edge {100*edge_a:.0f}%")
    d_e = M_eps[ie].sum() - best_e.sum(); d_a = M_a0[ia].sum() - best_a.sum()
    print(f"   dchi2 bought by non-universality: class {d_e:.0f}  MOND {d_a:.0f}  ratio {d_e/d_a:.2f}x")
    print(f"\n   {'observable':<12s} {'rho_s(eps0)':>12s} {'p':>10s}  |  {'rho_s(a0)':>10s} {'p':>10s}")
    corr = {}
    for k, v in obs.items():
        r1, p1 = spearmanr(v, le); r2, p2 = spearmanr(v, la)
        corr[k.strip()] = dict(eps0=[float(r1), float(p1)], a0=[float(r2), float(p2)])
        print(f"   {k:<12s} {r1:+12.3f} {p1:10.2e}  |  {r2:+10.3f} {p2:10.2e}")
    out = dict(global_eps0=float(EPS[ie]), global_chi2N_class=float(g_e), global_a0=float(A0S[ia]), global_chi2N_mond=float(g_a),
               pergal_chi2N_class=float(pg_e), pergal_chi2N_mond=float(pg_a), class_wins=wins, wilcoxon_p=float(pw),
               scatter_eps0_dex=sp_e, scatter_a0_dex=sp_a, edge_eps0=edge_e, edge_eps0_low=lo_e, edge_eps0_high=hi_e, edge_a0=edge_a,
               dchi2_ratio=float(d_e / d_a), spearman=corr,
               per_galaxy=dict(gid=[g.gid for g in G], eps0=e0g.tolist(), a0=a0g.tolist(), chi2_e=best_e.tolist(), chi2_a=best_a.tolist(),
                               n=NN.tolist(), logMbar=obs["log Mbar "].tolist()))
    if site_ref:
        print("\n   vs the site's printed values:")
        checks = [("global chi2/N class", g_e, site_ref["g_e"]), ("global chi2/N MOND", g_a, site_ref["g_a"]),
                  ("per-gal chi2/N class", pg_e, site_ref["pg_e"]), ("per-gal chi2/N MOND", pg_a, site_ref["pg_a"]),
                  ("class wins", 100 * wins, site_ref["wins"]), ("eps0 16-84 dex", sp_e, site_ref["sp_e"]),
                  ("a0 16-84 dex", sp_a, site_ref["sp_a"]), ("eps0 at-edge %", 100 * edge_e, site_ref["edge_e"]),
                  ("rho_s(eps0, logMbar)", corr["log Mbar"]["eps0"][0], site_ref["rs_e"]),
                  ("rho_s(a0, logMbar)", corr["log Mbar"]["a0"][0], site_ref["rs_a"]),
                  ("rho_s(eps0, log rho_m)", corr["log rho_m"]["eps0"][0], site_ref["rs_rho"])]
        allok = True
        for name, mine, theirs in checks:
            ok = abs(mine - theirs) <= max(0.02 * abs(theirs), 0.006)
            allok &= ok
            print(f"      {name:<24s} here {mine:9.3f}  site {theirs:9.3f}  {'ok' if ok else 'DIFF'}")
        print(f"   ALL WITHIN 2% (or 0.006 abs): {allok}")
        out["site_reproduced"] = bool(allok)
    return out, le, la


SITE_REF = dict(g_e=126.53, g_a=52.22, pg_e=39.70, pg_a=10.30, wins=21, sp_e=1.197, sp_a=0.619, edge_e=42, rs_e=0.758, rs_a=0.073, rs_rho=0.162)
res["A_site_grid"], le_s, la_s = analyse("A. REPRODUCTION on the site's own grids", EPS_SITE, M_eps_s, A0_SITE, M_a0_s, SITE_REF)
res["B_uncensored"], le_w, la_w = analyse("B. EXTENSION -- uncensored grids (the seed's step 1)", EPS_WIDE, M_eps_w, A0_WIDE, M_a0_w)

# --------------------------------------------------------------- C. permutation null (the seed's trap)
hdr("C. PERMUTATION NULL for rho_s(log eps0, log M_bar) -- galaxy labels shuffled (the seed: 'declare the null by permutation')")
rng = np.random.default_rng(20260901)
lm = obs["log Mbar "]
NPERM = 20000
perm = {}
for tag, le in (("site grid", le_s), ("uncensored", le_w)):
    obs_r = spearmanr(lm, le)[0]
    null = np.array([spearmanr(lm, rng.permutation(le))[0] for _ in range(NPERM)])
    p_perm = float((np.sum(np.abs(null) >= abs(obs_r)) + 1) / (NPERM + 1))
    q999 = float(np.percentile(np.abs(null), 99.9))
    print(f"   {tag:<11s} observed rho_s = {obs_r:+.3f}   null |rho_s| 99.9th pct = {q999:.3f}   max |null| = {np.abs(null).max():.3f}   p_perm <= {p_perm:.1e}  ({NPERM} perms)")
    perm[tag] = dict(observed=float(obs_r), null_q999=q999, null_max=float(np.abs(null).max()), p_perm_upper=p_perm)
    # the same null for the a0 control
    la = la_s if tag == "site grid" else la_w
    obs_a = spearmanr(lm, la)[0]
    null_a = np.array([spearmanr(lm, rng.permutation(la))[0] for _ in range(2000)])
    print(f"   {'':<11s} control a0: observed rho_s = {obs_a:+.3f}   fraction of null >= |obs| = {np.mean(np.abs(null_a) >= abs(obs_a)):.2f}")
    perm[tag]["a0_control"] = dict(observed=float(obs_a), frac_null_ge=float(np.mean(np.abs(null_a) >= abs(obs_a))))
res["C_permutation"] = perm

# --------------------------------------------------------------- D. the seed's step 2+3 AS NUMBERS, not a verdict
hdr("D. log eps0 = k*log M_bar + b, and the same for a0 -- scatter about the fit (seed steps 2+3; NO decision rule pre-registered here)")
fits = {}
for tag, le, la, arg_e_edge in (("site grid", le_s, la_s, None), ("uncensored", le_w, la_w, None)):
    for name, y in (("eps0", le), ("a0", la)):
        A = np.vstack([lm, np.ones_like(lm)]).T
        k, b = np.linalg.lstsq(A, y, rcond=None)[0]
        resid = y - (k * lm + b)
        rms = float(np.sqrt(np.mean(resid**2))); mad = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
        sd_raw = float(np.std(y))
        print(f"   {tag:<11s} {name:<5s} slope k = {k:+.3f} dex/dex   intercept {b:+.2f}   scatter about fit: rms {rms:.3f} dex, robust {mad:.3f} dex   (raw sd {sd_raw:.3f}; fraction removed {(1-rms/sd_raw)*100:.0f}%)")
        fits[f"{tag}/{name}"] = dict(k=float(k), b=float(b), rms=rms, mad=mad, sd_raw=sd_raw)
print("   (MOND's own relation-scatter benchmark, quoted by the seed for comparison only: the RAR's 0.11 dex.)")
res["D_fits"] = fits
res["elapsed_s"] = time.time() - t0
with open(OUT_JSON, "w") as f:
    json.dump(res, f, indent=1)
print(f"\n[results -> {os.path.basename(OUT_JSON)}]   total {time.time()-t0:.0f}s")
