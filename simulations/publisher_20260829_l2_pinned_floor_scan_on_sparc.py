#!/usr/bin/env python3
"""
Publisher 2026-08-29 -- REPRODUCE the site explorer's 08-28 "L2 solved on real SPARC" numbers in-repo,
then EXTEND them with the one number the RG literature still owed the sector: how much the fit
moves when the floor is PINNED at the derived Omega_m = 0.315 versus left on a fine grid.

Depends on the site lane's solver, committed in synchronism-site 4b9c681 (explorer: session 2026-08-28):
  explorer/findings/scripts/l2_sparc_core.py, l2_field_equation_on_sparc.py
which in turn read the in-repo SPARC files simulations/sparc_real_data/{MassModels_Lelli2016c,SPARC_Lelli2016c}.mrt.
The site's OUTPUT files for this run are uncommitted in the shared checkout (modified 2026-08-28 09:07 PDT after the
08:20 commit; no finding write-up exists) -- which is why the numbers are re-derived here rather than cited.

Field equation (L2):  div[ C(rho) grad Phi ] = 4 pi G rho, axisymmetric finite volume, Dirichlet outer BC,
midplane g_R interpolated to the SPARC radii; boost B = g_L2/g_N transferred onto SPARC's V_bar^2.
Forms:  framework  C = floor + (1-floor) tanh(gamma ln(rho/rho_c + 1))          (gamma = 0.489)
        RG         eps = eps0 + (1-eps0) 0.5 [tanh(Q ln(rho/rho_c)) + 1]       (Q = 0.47)
Scoring: chi2 on V_obs with SPARC errors; Upsilon_disk profiled over {0.3,0.5,0.7} with a 0.1-dex lognormal prior
(site convention).  chi2/N below is PER POINT (the site's header says "N_eff = galaxies", its code divides by points).
"""
import os, sys, time, json
import numpy as np
from multiprocessing import get_context

SITE = "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/findings/scripts"
sys.path.insert(0, SITE)
import l2_sparc_core as K                 # noqa: E402
import l2_field_equation_on_sparc as M    # noqa: E402

OM = 0.315
t0 = time.time()
gal = K.load_sparc()
gids = sorted(gal)
G = []
GU = {u: [] for u in M.UPS_GRID}
for gid in gids:
    G.append(M.Gal(gal[gid]))
    for u in M.UPS_GRID:
        GU[u].append(M.Gal(gal[gid], up_disk=u))
NPTS = sum(int(g.ok.sum()) for g in G)
print(f"built {len(G)} galaxies x {1+len(M.UPS_GRID)} Upsilon values in {time.time()-t0:.0f}s; "
      f"{NPTS} scored points (transfer-guard dropped {sum(g.n_dropped for g in G)})", flush=True)


def make_C(kind, ex, rc, fl):
    return K.C_framework(ex, rc, fl) if kind == "fw" else K.C_refracted(fl, ex, rc)


def score_point(args):
    kind, ex, rc, fl = args
    Cf = make_C(kind, ex, rc, fl)
    per = []
    for i in range(len(G)):
        sc = []
        for u in M.UPS_GRID:
            gu = GU[u][i]
            sc.append(gu.score(gu.solve(Cf, fl)[0] * gu.vbar2))
        per.append(M.profiled(sc))
    s = M.summarise(f"{kind} ex={ex} rc={rc:.3g} fl={fl}", per)
    return args, s


def score_ref(f):
    per = []
    for i in range(len(G)):
        per.append(M.profiled([GU[u][i].score(f(GU[u][i])) for u in M.UPS_GRID]))
    return M.summarise("", per)


refs = {
    "Newton (C=1)": lambda gg: gg.vbar2,
    "MOND simple mu (a0=1.2e-10)": lambda gg: K.mond_simple(gg.gbar_sparc) * gg.d["R"],
    "MOND RAR nu (McGaugh+16)": lambda gg: K.mond_rar(gg.gbar_sparc) * gg.d["R"],
}
REF = {k: score_ref(f) for k, f in refs.items()}
mond = np.array(REF["MOND simple mu (a0=1.2e-10)"]["per_gal_chi2"])

# ---- 1. REPRODUCTION rows (site's printed values, Upsilon profiled column) ----
SITE_ROWS = {
    ("fw", 0.489, 0.161, OM):     ("L2 Jeans gamma=.489 rho_c=0.161 floor=Om",         160.27),
    ("fw", 2.0,   0.161, OM):     ("L2 Jeans gamma=2 rho_c=0.161 floor=Om",            107.89),
    ("RG", 1.79,  4.3e-3, 0.661): ("L2 RG DMS-unique e0=.661 Q=1.79 rc=4.3e-3",        239.92),
    ("RG", 0.47,  8.3e-3, 0.089): ("L2 RG E0 (site label 'Cesare+20') e0=.089 q=.47", 911.46),
    ("fw", 0.489, 10**-3.5, OM):  ("grid framework floor=0.315 rho_c=3.16e-4",         69.67),
    ("RG", 0.47,  10**(-3.5+4/6), OM): ("grid RG floor=0.315 rho_c=1.47e-3",           68.44),
    ("RG", 0.47,  10**(-3.5+3*4/6), 0.661): ("grid RG floor=0.661 rho_c=3.16e-2",      236.15),
}
SITE_REFS = {"MOND simple mu (a0=1.2e-10)": 21.25, "MOND RAR nu (McGaugh+16)": 21.55, "Newton (C=1)": 465.43}

# ---- 2. EXTENSION: fine floor scan at the site's best rho_c and its two grid neighbours ----
FLOORS = [0.15, 0.20, 0.25, 0.315, 0.40, 0.50, 0.661]
RCS = [10**-3.5, 10**(-3.5+4/6), 10**(-3.5+2*4/6)]          # 3.16e-4, 1.47e-3, 6.81e-3 Msun/pc^3
SCAN = [(kind, ex, rc, fl) for kind, ex in (("RG", 0.47), ("fw", 0.489)) for rc in RCS for fl in FLOORS]

points = list(SITE_ROWS.keys()) + [p for p in SCAN if p not in SITE_ROWS]
print(f"scoring {len(points)} (form, exponent, rho_c, floor) points on {os.cpu_count()} cpus", flush=True)
t1 = time.time()
try:
    with get_context("fork").Pool(min(6, os.cpu_count() or 1)) as pool:
        out = pool.map(score_point, points)
except Exception as e:  # noqa
    print("pool failed, running serially:", e, flush=True)
    out = [score_point(p) for p in points]
RES = {args: s for args, s in out}
print(f"scored in {time.time()-t1:.0f}s\n", flush=True)

print("=" * 78 + "\n1. REPRODUCTION of the site's 08-28 numbers (chi2/N per point, Upsilon profiled)\n" + "=" * 78)
ok_all = True
for k, v in SITE_REFS.items():
    mine = REF[k]["chi2_per_pt"]; flag = abs(mine - v) / v < 0.02
    ok_all &= flag
    print(f"   {k:<52s} site {v:8.2f}   here {mine:8.2f}   {'PASS' if flag else 'FAIL'}")
for args, (label, v) in SITE_ROWS.items():
    mine = RES[args]["chi2_per_pt"]; flag = abs(mine - v) / v < 0.02
    ok_all &= flag
    print(f"   {label:<52s} site {v:8.2f}   here {mine:8.2f}   {'PASS' if flag else 'FAIL'}")
print(f"   ALL WITHIN 2%: {ok_all}")

print("\n" + "=" * 78 + "\n2. EXTENSION: floor scan under L2 (chi2/N per point; * = best in column)\n" + "=" * 78)
best_overall = None
for kind, ex in (("RG", 0.47), ("fw", 0.489)):
    print(f"   --- {kind} form (exponent {ex}) ---")
    print("   floor \\ rho_c " + " ".join(f"{r:9.3g}" for r in RCS))
    tab = np.array([[RES[(kind, ex, rc, fl)]["chi2_per_pt"] for rc in RCS] for fl in FLOORS])
    for j, fl in enumerate(FLOORS):
        row = " ".join(f"{tab[j, i]:9.2f}" + ("*" if j == np.argmin(tab[:, i]) else " ") for i in range(len(RCS)))
        print(f"   {fl:5.3f}         {row}")
    jm, im = np.unravel_index(np.argmin(tab), tab.shape)
    print(f"   scan minimum: floor={FLOORS[jm]}, rho_c={RCS[im]:.3g}, chi2/N={tab[jm, im]:.2f}")
    for i, rc in enumerate(RCS):
        col = tab[:, i]; jb = int(np.argmin(col)); j315 = FLOORS.index(0.315)
        print(f"   rho_c={rc:.3g}: pinned floor 0.315 vs column best (floor {FLOORS[jb]}): "
              f"Delta chi2 = {(col[j315]-col[jb])*NPTS:+.0f} over {NPTS} points; "
              f"vs floor 0.661: {(col[j315]-col[FLOORS.index(0.661)])*NPTS:+.0f}; "
              f"vs floor 0.15: {(col[j315]-col[FLOORS.index(0.15)])*NPTS:+.0f}")
    if best_overall is None or tab[jm, im] < best_overall[1]:      # (bug fixed after first run: the rc loop had shadowed i)
        best_overall = ((kind, ex, RCS[im], FLOORS[jm]), tab[jm, im])

args, val = best_overall
c = np.array(RES[args]["per_gal_chi2"])
print(f"\n   BEST scanned point overall: {args}  chi2/N = {val:.2f}  "
      f"(MOND simple {REF['MOND simple mu (a0=1.2e-10)']['chi2_per_pt']:.2f}, "
      f"RAR nu {REF['MOND RAR nu (McGaugh+16)']['chi2_per_pt']:.2f}, Newton {REF['Newton (C=1)']['chi2_per_pt']:.2f})")
print(f"   vs MOND simple at galaxy level: wins {np.mean(c < mond)*100:.0f}% of {len(c)}; "
      f"median log10(chi2/chi2_MOND) = {np.median(np.log10(np.maximum(c,1e-3)/np.maximum(mond,1e-3))):+.2f}; "
      f"sum dchi2 = {c.sum()-mond.sum():+.0f}")
from scipy.stats import wilcoxon
def vs_mond(label, cc):
    st = wilcoxon(np.log(np.maximum(cc, 1e-3)), np.log(np.maximum(mond, 1e-3)))
    print(f"   {label}: wins {np.mean(cc < mond)*100:.0f}% of {len(cc)} galaxies vs MOND simple; "
          f"median log10(chi2/chi2_MOND) = {np.median(np.log10(np.maximum(cc,1e-3)/np.maximum(mond,1e-3))):+.2f}; "
          f"Wilcoxon p = {st.pvalue:.2e}")
vs_mond("best scanned point", c)
vs_mond("pinned floor 0.315 (RG form, rho_c=1.47e-3)", np.array(RES[("RG", 0.47, RCS[1], OM)]["per_gal_chi2"]))
vs_mond("pinned floor 0.315 (fw form, rho_c=3.16e-4)", np.array(RES[("fw", 0.489, RCS[0], OM)]["per_gal_chi2"]))

# ---- 3. the floor-marginal: minimum over rho_c for each floor, both forms ----
print("\n" + "=" * 78 + "\n3. FLOOR-MARGINAL (min over the 3 rho_c values), both forms; Delta chi2 relative to 0.315\n" + "=" * 78)
for kind, ex in (("RG", 0.47), ("fw", 0.489)):
    marg = {fl: min(RES[(kind, ex, rc, fl)]["chi2_per_pt"] for rc in RCS) for fl in FLOORS}
    base = marg[0.315]
    print(f"   {kind}: " + "  ".join(f"{fl}:{marg[fl]:.1f}({(marg[fl]-base)*NPTS:+.0f})" for fl in FLOORS))

json.dump({"npts": NPTS, "refs": {k: v["chi2_per_pt"] for k, v in REF.items()},
           "points": {str(k): v["chi2_per_pt"] for k, v in RES.items()}},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "publisher_20260829_l2_pinned_floor_scan_results.json"), "w"), indent=1)
print("\nexit 0")
