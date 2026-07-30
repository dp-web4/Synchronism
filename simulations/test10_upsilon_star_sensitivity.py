#!/usr/bin/env python3
"""
TEST-10 class-exclusion headline: stellar mass-to-light (Upsilon*) sensitivity.
Publisher, 2026-07-30.

WHY. The 2026-07-29 headline is "a bounded-boost class with B_max <= 6.4 is excluded by
28 SPARC galaxies, under every candidate ceiling normalization." It is convention-immune
by construction (evaluated at the most permissive candidate ceiling Om/Ob = 6.39). It is
NOT yet M/L-immune: f_DM = 1 - vbar^2/vobs^2 depends on the assumed stellar M/L, which
the executed run fixes at the SPARC-standard Upsilon*_disk = 0.5, Upsilon*_bul = 0.7
(Lelli+2016, 3.6um). A referee will ask, and the check is cheap.

Direction of the threat: raising Upsilon* raises vbar, lowers f_DM, lowers the required
boost B_req = 1/(1 - f_DM), and therefore SHRINKS the exclusion. So the adversarial
question is: how large must Upsilon* be before the exclusion disappears, and is there a
subsample for which it never disappears?

Mechanism that makes the answer non-obvious: the gas term vgas|vgas| does NOT scale with
Upsilon*. In gas-dominated dwarfs -- which are exactly the most DM-dominated objects in
SPARC and therefore the ones carrying the exclusion -- f_DM is nearly Upsilon*-independent.
If that holds quantitatively, the exclusion has an M/L-immune core.

NOT the registered ceiling-definition sweep (`Research/proposals/boost_ceiling_provenance_
and_class_exclusion.md`), which also recomputes TEST-09's BTFR slope per definition under a
pre-fixed verdict rule and remains UNRUN. This does not discharge it.

Reuses the executed script's loader and cuts by import, so the sample is identical to the
2026-07-15 run (reproduced 2026-07-29: 106/153 = 69%, f_DM,max = 0.927).
"""
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/scripts")
import test10_dwarf_dm_fraction_ceiling as t10  # noqa: E402

OM, OB = 0.315, 0.0493
B_PERMISSIVE = OM / OB           # 6.39 -- most permissive candidate ceiling
B_SITE = 1.0 / OM                # 3.17 -- the site's choice
SEED = 20260730

# Upsilon*_disk grid. SPARC standard is 0.5 with ~0.1 dex scatter (Lelli+2016 / McGaugh
# & Schombert 2014); 0.7 is the population-synthesis high end; >=1.0 is already beyond
# any defensible 3.6um prior and is included only to locate where the claim would break.
UPSILON_GRID = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 3.0]
NOMINAL = 0.5
BUL_RATIO = 0.7 / 0.5            # hold Upsilon*_bul / Upsilon*_disk at the SPARC ratio


def build_rows(up_disk):
    """t10's pipeline with Upsilon*_disk as a free parameter (bulge scaled at fixed ratio).

    Returns (name, f_DM, gas_share) where gas_share is the gas fraction of vbar^2 at the
    NOMINAL M/L -- a fixed per-galaxy label, so the same galaxies are tagged at every grid
    point.
    """
    up_bul = up_disk * BUL_RATIO
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
        gas2 = vgas * np.abs(vgas)
        star2 = up_disk * vdisk ** 2 + up_bul * vbul ** 2
        vbar2 = gas2 + star2
        # t10 rejects galaxies with any non-positive vbar^2; keep that cut at NOMINAL M/L
        # so the SAMPLE is fixed across the grid and only f_DM moves.
        nom = gas2 + NOMINAL * vdisk ** 2 + NOMINAL * BUL_RATIO * vbul ** 2
        if (nom <= 0).any():
            continue
        w = 1.0 / evobs ** 2
        vobs_o = np.sqrt(np.average(vobs ** 2, weights=w))
        vbar2_o = np.average(vbar2, weights=w)
        f_obs = 1.0 - vbar2_o / vobs_o ** 2
        gas_share = np.average(gas2, weights=w) / np.average(nom, weights=w)
        rows.append((name, f_obs, gas_share))
    return rows


def b_req_of(rows):
    f = np.array([r[1] for r in rows])
    # f -> 1 means B_req -> inf; f <= 0 means baryons already over-account (no boost needed)
    return np.where(f < 1.0, 1.0 / np.clip(1.0 - f, 1e-9, None), np.inf)


def main():
    base = build_rows(NOMINAL)
    names = [r[0] for r in base]
    gas_share = np.array([r[2] for r in base])
    n = len(base)

    # ---- reproduction check against the executed run -------------------------------
    b0 = b_req_of(base)
    print(f"Reproduction @ Upsilon*_disk = {NOMINAL}: N = {n}, "
          f"{(b0 > B_SITE).sum()}/{n} exceed 1/Om = {B_SITE:.2f}, "
          f"f_DM,max = {1 - 1/b0.max():.3f}  (expect 106/153, 0.927)")

    # ---- the sweep -----------------------------------------------------------------
    print(f"\nExceedance vs assumed stellar M/L, at the most permissive ceiling "
          f"B_max = Om/Ob = {B_PERMISSIVE:.2f}")
    print(f"{'Ups*_disk':>10} {'median f_DM':>12} {'K @6.39':>9} {'%':>5} "
          f"{'K @3.17':>9} {'B_max excl by 10':>18}")
    curves = {}
    for u in UPSILON_GRID:
        rows = build_rows(u)
        b = b_req_of(rows)
        curves[u] = b
        f = np.array([r[1] for r in rows])
        k_perm = int((b > B_PERMISSIVE).sum())
        k_site = int((b > B_SITE).sum())
        b10 = np.sort(b)[-10]
        tag = "   <- SPARC standard" if u == NOMINAL else ""
        print(f"{u:>10.1f} {np.median(f):>12.3f} {k_perm:>7}/{n} {100*k_perm/n:>4.0f}% "
              f"{k_site:>7}/{n} {b10:>18.2f}{tag}")

    # ---- where does the claim break? ------------------------------------------------
    print("\nCritical M/L: smallest grid Upsilon*_disk at which the exclusion drops to K")
    for target in (28, 20, 10, 5, 1, 0):
        hit = next((u for u in UPSILON_GRID
                    if int((curves[u] > B_PERMISSIVE).sum()) <= target), None)
        print(f"  K <= {target:>2}  at Upsilon*_disk >= {hit if hit is not None else '> 3.0'}")

    # ---- the M/L-immune core ---------------------------------------------------------
    # Galaxies that exceed the most permissive ceiling at EVERY grid point, i.e. whose
    # verdict no M/L choice in the scanned range can overturn.
    stack = np.vstack([curves[u] > B_PERMISSIVE for u in UPSILON_GRID])
    always = stack.all(axis=0)
    defensible = [u for u in UPSILON_GRID if 0.3 <= u <= 0.7]
    stack_d = np.vstack([curves[u] > B_PERMISSIVE for u in defensible])
    always_d = stack_d.all(axis=0)
    print(f"\nM/L-immune core (exceed B_max = {B_PERMISSIVE:.2f} at EVERY Upsilon* tested):")
    print(f"  over the defensible prior 0.3-0.7 : {int(always_d.sum())}/{n} galaxies")
    print(f"  over the full grid 0.2-3.0        : {int(always.sum())}/{n} galaxies")
    if always.any():
        idx = np.argsort(-curves[3.0])
        core = [names[i] for i in idx if always[i]][:12]
        print(f"  core members (by B_req at Ups*=3.0): {', '.join(core)}")
        print(f"  their median gas share of vbar^2 at nominal M/L: "
              f"{np.median(gas_share[always]):.2f}  "
              f"(sample median {np.median(gas_share):.2f})")

    # ---- the mechanism, stated as a number -------------------------------------------
    hi, lo = gas_share > np.median(gas_share), gas_share <= np.median(gas_share)
    d_hi = np.median(curves[1.0][hi] / curves[0.5][hi])
    d_lo = np.median(curves[1.0][lo] / curves[0.5][lo])
    print(f"\nMechanism check -- doubling Upsilon* (0.5 -> 1.0) multiplies B_req by:")
    print(f"  gas-dominated half (share > {np.median(gas_share):.2f}): x{d_hi:.2f}")
    print(f"  star-dominated half                          : x{d_lo:.2f}")


if __name__ == "__main__":
    main()
