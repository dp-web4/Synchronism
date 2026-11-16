#!/usr/bin/env python3
"""Publisher 2026-11-16: does the dSph exclusion of the Omega_m/Omega_b = 6.39 boost cap
survive the mass-estimator convention?

Source: synchronism-site/maintainer/scripts/cap6p39_vs_dsph.py (PREREG site c7d729c), inputs
copied unchanged from explorer/findings/scripts/pressure_supported_boost_ceiling.py (2026-08-23).
That script uses sigma_N^2 = G*(M/2)/(4*r_half) with r_half labelled 3D. Wolf et al. 2010
(MNRAS 406, 1220, eq. 2) give M_1/2 = 3*<sigma_los^2>*r_1/2/G for the 3D half-light radius;
the factor 4 belongs with the projected R_e (r_1/2 ~ 4/3 R_e). With the 3D radius, k = 4
understates sigma_N by sqrt(3/4), so it overstates the required boost by 4/3 and the z by a
little. This script reruns the 6.39 test at k = 4 (identity control) and k = 3, and adds a
stress test with L_V raised 30% (no catalogue value was read; this only states a margin).
"""
import math

G = 4.300917270e-6  # kpc (km/s)^2 / Msun
SYS = {  # name: L_V, r_half_3D kpc, sigma_obs, err, ML_fid, ML_hi   (2026-08-23 table)
    "Crater II": (1.6e5, 1.42, 2.7, 0.3, 2.0, 4.0),
    "Draco":     (2.2e5, 0.29, 9.1, 1.2, 2.0, 4.0),
}


def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# six-system one-sided equivalent of 3 sigma, computed as the source script does (prints 3.51;
# the unrounded value matters for the identity control's 14.7)
_pS = 1 - (1 - (1 - _phi(3.0))) ** (1 / 6)
_lo, _hi = 0.0, 10.0
for _ in range(80):
    _m = (_lo + _hi) / 2
    (_lo, _hi) = (_m, _hi) if 1 - _phi(_m) > _pS else (_lo, _m)
Z_SIDAK = (_lo + _hi) / 2


def sig_n(M, rh, k):
    return math.sqrt(G * (M / 2.0) / (k * rh))


def z(name, cap, ml, k, lscale=1.0):
    L, rh, s, e, _, _ = SYS[name]
    return (s - sig_n(L * lscale * ml, rh, k) * math.sqrt(cap)) / e


def largest_cap(name, ml, k, zt, lscale=1.0):
    L, rh, s, e, _, _ = SYS[name]
    return ((s - zt * e) / sig_n(L * lscale * ml, rh, k)) ** 2


print("IDENTITY CONTROL (k = 4 must reproduce cap6p39_vs_dsph_output.txt)")
ref = {("Draco", "fid"): 5.7, ("Draco", "hi"): 4.9, ("Crater II", "fid"): 6.1, ("Crater II", "hi"): 4.9}
ok = True
for (n, w), r in ref.items():
    ml = SYS[n][4] if w == "fid" else SYS[n][5]
    v = round(z(n, 6.39, ml, 4), 1)
    ok &= abs(v - r) < 1e-9
    print(f"  {n:<10} M/L {w:<3} z(6.39) = {v:4.1f}  ref {r:4.1f}")
for n, r in (("Draco", 14.7), ("Crater II", 11.2)):
    v = round(largest_cap(n, SYS[n][5], 4, Z_SIDAK), 1)
    ok &= abs(v - r) < 1e-9
    print(f"  {n:<10} largest cap excluded at z = {Z_SIDAK:.2f} (M/L hi) = {v:4.1f}  ref {r:4.1f}")
print(f"  identity control: {'PASS' if ok else 'FAIL'}\n")

print(f"{'case':<34}{'system':<11}{'z fid':>7}{'z M/L=4':>9}{'max cap (trials, M/L=4)':>26}")
for label, k, ls in (("k = 4 (as run)", 4, 1.0),
                     ("k = 3 (Wolf+2010, 3D r_1/2)", 3, 1.0),
                     ("k = 3 and L_V x 1.3 (stress)", 3, 1.3)):
    for n in SYS:
        _, _, _, _, mf, mh = SYS[n]
        print(f"{label:<34}{n:<11}{z(n, 6.39, mf, k, ls):>7.2f}{z(n, 6.39, mh, k, ls):>9.2f}"
              f"{largest_cap(n, mh, k, Z_SIDAK, ls):>26.1f}")
print(f"\nSidak threshold z = {Z_SIDAK:.4f}; 6.39 excluded band-robustly after trials iff z(M/L=4) > that.")
print("B_req at fiducial M/L, k = 3 (vs the 08-23 table's 101.5 / 60.2 at k = 4):")
for n in SYS:
    L, rh, s, e, mf, _ = SYS[n]
    print(f"  {n:<10} {(s / sig_n(L * mf, rh, 3)) ** 2:6.1f}")
