#!/usr/bin/env python3
"""
Which coupling is orientation-inverted: the whitepaper's or the site's?

Publisher, 2026-08-01.

WHY. `synchronism-site/explorer/findings/two-coherence-orientations-chemistry-uses-the-
flipped-one.md` (2026-07-29) concludes that "the governing equation makes coherence *rise*
with density, and every sector that touches data needs it to *fall*", and proposes a global
orientation flip C -> C(rho_crit/rho). That headline is stated about *the framework*.

But the framework has two different couplings of C to gravity in two different places:

  (A) whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md:36
          G_eff = G / C(rho)        =>  v^2 = v_bar^2 / C
      enhancement goes as 1/C, i.e. LOW C means MORE missing gravity.

  (B) synchronism-site /galaxy-plotter (code, quoted in the finding's section 1)
          v_syn^2 = v_bar^2 + (V_flat * C)^2
      enhancement goes as C^2, i.e. HIGH C means MORE missing gravity.

These demand OPPOSITE orientations from the same rotation curve. The finding computes the
required non-baryonic term against density (mean r = -0.824) and reads that as a refutation
of the governing equation -- but that reading only follows under coupling (B). Under (A) the
same anti-correlation is what a density-increasing C predicts.

So the test is: invert each coupling for the C it REQUIRES at each radius, and correlate that
required C against a local density proxy. Whichever coupling demands a C that FALLS with
density is the inverted one.

Density proxy: SPARC's own SBdisk + SBbul (L/pc^2, 3.6um) -- a directly measured local
surface brightness. No profile assumption, no estimator choice. Radius is reported as a
second, independent proxy (density falls outward in every disk) precisely because the
07-29 finding's own section 4 shows estimator choice has silently carried headline verdicts
on this program three times.

Sample and cuts identical to the executed TEST-10 script (Q<=2, inc>30), loader imported
from it so the sample cannot drift.
"""
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/scripts")
import test10_dwarf_dm_fraction_ceiling as t10  # noqa: E402

BASE = t10.BASE
UP_DISK, UP_BUL = t10.UP_DISK, t10.UP_BUL


def load_with_sb():
    """Table 2 incl. surface brightness: R, Vobs, eVobs, Vgas, Vdisk, Vbul, SBdisk, SBbul."""
    out = {}
    with open(BASE + "MassModels_Lelli2016c.mrt") as f:
        for line in f:
            p = line.split()
            if len(p) < 10:
                continue
            try:
                name = p[0]
                vals = [float(x) for x in p[2:10]]
            except ValueError:
                continue
            out.setdefault(name, []).append(vals)
    return {k: np.array(v) for k, v in out.items()}


def spearman(x, y):
    """Rank correlation without scipy."""
    def rank(a):
        order = np.argsort(a, kind="mergesort")
        r = np.empty(len(a), float)
        r[order] = np.arange(len(a), dtype=float)
        # average ties
        _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
        for i in np.where(cnt > 1)[0]:
            m = inv == i
            r[m] = r[m].mean()
        return r
    rx, ry = rank(np.asarray(x, float)), rank(np.asarray(y, float))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d > 0 else np.nan


def main():
    meta = t10.load_meta()
    mm = load_with_sb()

    rows = []
    for name, arr in sorted(mm.items()):
        m = meta.get(name)
        if m is None or m["q"] > 2 or m["inc"] < 30.0:
            continue
        arr = arr[np.argsort(arr[:, 0])]
        r, vobs, evobs, vgas, vdisk, vbul, sbd, sbb = arr.T

        vbar2 = (vgas * np.abs(vgas)
                 + UP_DISK * vdisk * np.abs(vdisk)
                 + UP_BUL * vbul * np.abs(vbul))

        sigma = sbd + sbb                      # local surface density proxy [L/pc^2]
        ok = (vobs > 0) & (vbar2 > 0) & (sigma > 0) & (vobs ** 2 > vbar2)
        if ok.sum() < 5:
            continue

        r, vobs, vbar2, sigma = r[ok], vobs[ok], vbar2[ok], sigma[ok]

        # V_flat proxy for coupling (B): outer-3 mean of Vobs (the site plotter's own scale).
        v_flat = vobs[-3:].mean()

        # Required C, coupling (A): v^2 = v_bar^2 / C   =>  C = v_bar^2 / v_obs^2
        C_A = vbar2 / vobs ** 2
        # Required C, coupling (B): v^2 = v_bar^2 + (V_flat C)^2 => C = sqrt(v_obs^2-v_bar^2)/V_flat
        C_B = np.sqrt(vobs ** 2 - vbar2) / v_flat

        rows.append(dict(
            name=name, n=int(ok.sum()),
            sA=spearman(np.log(sigma), C_A), sB=spearman(np.log(sigma), C_B),
            rA=spearman(r, C_A), rB=spearman(r, C_B),
        ))

    sA = np.array([x["sA"] for x in rows])
    sB = np.array([x["sB"] for x in rows])
    rA = np.array([x["rA"] for x in rows])
    rB = np.array([x["rB"] for x in rows])

    print(f"galaxies passing cuts with >=5 usable radii: {len(rows)}")
    print()
    print("Spearman of the REQUIRED coherence against the local density proxy")
    print("(positive = required C rises with density = agrees with C(rho) as written)")
    print()
    print(f"{'coupling':<52} {'median':>8} {'mean':>8} {'frac>0':>8}")
    print("-" * 80)
    print(f"{'(A) whitepaper  G_eff = G/C   [C = vbar^2/vobs^2]':<52} "
          f"{np.median(sA):>8.3f} {sA.mean():>8.3f} {np.mean(sA > 0):>8.1%}")
    print(f"{'(B) site plotter  vsyn^2 = vbar^2 + (Vflat C)^2':<52} "
          f"{np.median(sB):>8.3f} {sB.mean():>8.3f} {np.mean(sB > 0):>8.1%}")
    print()
    print("Same, against radius (independent proxy; density falls outward)")
    print(f"{'(A) whitepaper':<52} {np.median(rA):>8.3f} {rA.mean():>8.3f} "
          f"{np.mean(rA > 0):>8.1%}")
    print(f"{'(B) site plotter':<52} {np.median(rB):>8.3f} {rB.mean():>8.3f} "
          f"{np.mean(rB > 0):>8.1%}")
    print()
    print("(radius correlations have the OPPOSITE expected sign to the density ones,")
    print(" since radius is an inverse density proxy)")
    print()

    # the five galaxies the 07-29 finding tabulated
    print("the five galaxies tabulated in the 07-29 finding:")
    print(f"{'galaxy':<14}{'n':>4}{'  (A) vs logSigma':>20}{'  (B) vs logSigma':>20}")
    for want in ["DDO154", "NGC2403", "NGC3198", "UGC128", "NGC7331"]:
        for x in rows:
            if x["name"].replace(" ", "").upper() == want:
                print(f"{x['name']:<14}{x['n']:>4}{x['sA']:>20.3f}{x['sB']:>20.3f}")


if __name__ == "__main__":
    main()
