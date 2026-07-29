#!/usr/bin/env python3
"""
TEST-10 headline-statistic robustness diagnostic (Publisher, 2026-07-29).

NOT the registered ceiling-definition sweep. `Research/proposals/boost_ceiling_provenance_
and_class_exclusion.md` registers a sweep that recomputes BOTH TEST-09's BTFR slope and
TEST-10's exceedance fraction per candidate ceiling definition, under a pre-fixed verdict
rule. That sweep remains UNRUN. This script answers a narrower question that does not need
it, and does not discharge it.

The question: the 2026-07-28 triage caveated "69% exceed the 68.5% ceiling" as
convention-dependent (a boost is a dynamical/baryonic ratio, whose cosmic value is
Om/Ob = 6.40, not 1/Om = 3.17) and prescribed citing f_DM,max = 0.927 (=> B >= 13.7) as
the definition-free replacement. f_DM,max is a SINGLE-GALAXY extreme-value statistic.
Trading a convention-dependence for an outlier-dependence is not obviously an improvement,
and a referee who dislikes the first will dislike the second.

So: compute the full exceedance curve B_req = 1/(1 - f_DM) over the same sample, and ask
whether there is a statistic that is BOTH definition-free (survives every candidate
ceiling) AND robust (backed by many galaxies, not one).

Reuses the executed script's loader and cuts verbatim by import, so the sample is identical
to the 2026-07-15 run (reproduced 2026-07-29: 106/153 = 69%, f_DM,max = 0.927).
"""
import sys
import numpy as np

sys.path.insert(0, "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/scripts")
import test10_dwarf_dm_fraction_ceiling as t10  # noqa: E402

OM, OB = 0.315, 0.0493
SEED = 20260729


def build_rows():
    """Identical pipeline to t10.main(), returning per-galaxy rows instead of printing."""
    meta, mm = t10.load_meta(), t10.load_mass_models()
    rows = []
    for name, arr in mm.items():
        m = meta.get(name)
        if m is None or m["q"] > 2 or m["inc"] < 30.0:
            continue
        arr = arr[np.argsort(arr[:, 0])]
        outer = arr[-3:] if len(arr) >= 3 else arr
        r, vobs, evobs, vgas, vdisk, vbul = outer.T
        ok = (vobs > 0) & (evobs > 0)
        if not ok.any():
            continue
        r, vobs, evobs, vgas, vdisk, vbul = (a[ok] for a in (r, vobs, evobs, vgas, vdisk, vbul))
        vbar2 = vgas * np.abs(vgas) + t10.UP_DISK * vdisk ** 2 + t10.UP_BUL * vbul ** 2
        if (vbar2 <= 0).any():
            continue
        w = 1.0 / evobs ** 2
        vobs_o = np.sqrt(np.average(vobs ** 2, weights=w))
        vbar2_o = np.average(vbar2, weights=w)
        f_obs = 1.0 - vbar2_o / vobs_o ** 2
        m_bar = (t10.UP_DISK * m["l36"] + 1.33 * m["mhi"]) * 1e9
        rows.append((name, m_bar, f_obs))
    return rows


def main():
    rows = build_rows()
    f = np.array([r[2] for r in rows])
    n = len(f)
    # required boost to produce the observed apparent DM fraction
    b_req = 1.0 / (1.0 - f)

    print(f"Sample: N = {n} SPARC galaxies (Q<=2, i>30), outer points. "
          f"Sanity: {(f > 1 - OM).sum()}/{n} exceed 1-Om = {1-OM:.3f}\n")

    print("Candidate ceilings and how many galaxies exceed each")
    print(f"{'definition':<28} {'B_max':>7} {'f_DM cap':>9} {'exceed':>10} {'%':>6}")
    cands = [
        ("1/Om  (site's choice)", 1.0 / OM),
        ("(Om-Ob)/Ob", (OM - OB) / OB),
        ("Om/Ob  (baryon budget)", OM / OB),
        ("-- 2x the most permissive", 2.0 * OM / OB),
        ("f_DM,max  (proposal's cite)", 1.0 / (1.0 - f.max())),
    ]
    for label, b in cands:
        k = int((b_req > b).sum())
        print(f"{label:<28} {b:>7.2f} {1-1/b:>9.3f} {k:>7}/{n} {100*k/n:>5.0f}%")

    print("\nRequired-boost quantiles (bootstrap 95% CI, 10k resamples)")
    rng = np.random.default_rng(SEED)
    boot = rng.choice(b_req, size=(10000, n), replace=True)
    for q in (50, 75, 90, 95, 99, 100):
        pt = np.percentile(b_req, q)
        lo, hi = np.percentile(np.percentile(boot, q, axis=1), [2.5, 97.5])
        tag = "  <- single galaxy" if q == 100 else ""
        print(f"  p{q:<3} B_req = {pt:>6.2f}   CI [{lo:.2f}, {hi:.2f}]{tag}")

    print("\nExcluded-ceiling statement at each support level")
    print("  'a bounded-boost class with B_max <= X is excluded by K SPARC galaxies'")
    for k in (1, 5, 10, 20, 40, 76):
        if k <= n:
            x = np.sort(b_req)[-k]
            print(f"    K = {k:>3} galaxies  ->  B_max <= {x:>6.2f} excluded")

    # how much of the claim rests on the single largest point
    b2 = np.sort(b_req)[-2]
    print(f"\nDrop-the-max sensitivity: B_req,max = {b_req.max():.2f} -> 2nd = {b2:.2f} "
          f"({100*(b_req.max()-b2)/b_req.max():.0f}% of the headline bound sits on one galaxy)")


if __name__ == "__main__":
    main()
