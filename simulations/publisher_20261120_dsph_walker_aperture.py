"""Publisher 2026-11-20: the dSph 6.39-cap exclusion under Walker+2009, with the light matched to Walker's aperture.

Source: synchronism-site/maintainer/scripts/cap6p39_vs_dsph_estimator.py (site maintainer 2026-11-19), routed as
PREDICTIONS.md back-annotation b26f6fa9 (i): "Under Walker+2009's smaller estimator Crater II drops below threshold
(2.9 sigma) and Draco passes it only just (3.7 sigma)".

That script writes Walker as M(R_e) = 2.5 sigma^2 R_e / G with R_e = (3/4) r_3D and then sets M(R_e) equal to the
half-LIGHT mass L*(M/L)/2, i.e. sigma_N^2 = G M_half / (1.875 r_3D). Walker+2009's M is the mass inside a SPHERE of
radius R_e (the projected half-light radius), and that sphere holds less than half the light: for a Plummer profile
(projected half-light radius a), the 3D fraction inside r is r^3/(r^2+a^2)^1.5, which is 2^-1.5 = 0.354 at r = a, and
r_1/2,3D = a / sqrt(2^(2/3) - 1) = 1.305 a. So the baryonic mass Walker's estimator should be fed is 0.354 L(M/L),
not 0.5 L(M/L), and R_e = r_3D / 1.305.

Same z, threshold and inputs as the source. z = (sigma_obs - sigma_N sqrt(B)) / err; Sidak 6-trial threshold 3.5087.
Identity control: the source's k = 4 / 3 / 1.875 rows must reproduce its printed output.
Limits: Plummer assumed for the aperture fraction (exponential/King profiles move it by a few per cent); quoted sigma
errors only; equilibrium estimator.
"""
import math

G = 4.300917270e-6  # kpc (km/s)^2 / Msun
ZS = 3.5087
SYS = [("Crater II", 1.6e5, 1.42, 2.7, 0.3), ("Draco", 2.2e5, 0.29, 9.1, 1.2)]
REF = {  # source output, (system, k): (z at M/L=2, z at M/L=4, escape M/L at 3.51)
    ("Crater II", 4.0): (6.1, 4.9, 7.0), ("Crater II", 3.0): (5.6, 4.2, 5.3), ("Crater II", 1.875): (4.7, 2.9, 3.3),
    ("Draco", 4.0): (5.7, 4.9, 9.2), ("Draco", 3.0): (5.4, 4.5, 6.9), ("Draco", 1.875): (4.8, 3.7, 4.3),
}


def sigma_n(L, ml, r3d, light_frac, denom_r):
    """sigma_N^2 = G * (light_frac * L * ml) / denom_r."""
    return math.sqrt(G * light_frac * L * ml / denom_r)


def row(L, r3d, s, e, light_frac, denom_coeff):
    denom = denom_coeff * r3d
    z = [(s - sigma_n(L, ml, r3d, light_frac, denom) * math.sqrt(6.39)) / e for ml in (2.0, 4.0)]
    s1 = sigma_n(L, 1.0, r3d, light_frac, denom)
    ml_esc = ((s - ZS * e) / math.sqrt(6.39) / s1) ** 2
    ml_esc3 = ((s - 3.0 * e) / math.sqrt(6.39) / s1) ** 2
    bx = ((s - ZS * e) / sigma_n(L, 4.0, r3d, light_frac, denom)) ** 2
    infl = (s - ZS * e) / (sigma_n(L, 4.0, r3d, light_frac, denom) * math.sqrt(6.39))
    return z, ml_esc, ml_esc3, bx, infl


print("IDENTITY CONTROL (half-light mass, sigma_N^2 = G M_half / (k r3D), as the source)")
ok = True
for name, L, r, s, e in SYS:
    for k in (4.0, 3.0, 1.875):
        z, esc, _, _, _ = row(L, r, s, e, 0.5, k)
        ref = REF[(name, k)]
        good = (round(z[0], 1), round(z[1], 1), round(esc, 1)) == ref
        ok &= good
        print(f"  {name:<10} k={k:<5} z(2)={z[0]:4.1f} z(4)={z[1]:4.1f} esc={esc:4.1f}  ref {ref}  {'ok' if good else 'MISMATCH'}")
print(f"  identity control: {'PASS' if ok else 'FAIL'}\n")

frac_re = 2 ** -1.5                       # Plummer light fraction inside a sphere of radius R_e
re_over_r3d = 1 / 1.0 * math.sqrt(2 ** (2 / 3) - 1)  # R_e / r_1/2,3D
print(f"Plummer: light inside sphere r = R_e = {frac_re:.4f}; R_e / r_3D = {re_over_r3d:.4f}")
k_eff = 2.5 * re_over_r3d * 0.5 / frac_re
print(f"Walker with its own aperture == half-light form with k_eff = {k_eff:.3f} (source used 1.875; Wolf 3)\n")

print(f"{'case':<44}{'system':<11}{'z M/L=2':>8}{'z M/L=4':>8}{'esc@3.51':>9}{'esc@3':>7}{'Bexcl':>7}{'infl':>6}")
cases = [("Wolf+2010 (3 r3D, half-light)", 0.5, 3.0),
         ("Walker+2009 as source (M_half, 1.875 r3D)", 0.5, 1.875),
         ("Walker+2009, light inside R_e (Plummer)", frac_re, 2.5 * re_over_r3d)]
for lab, frac, c in cases:
    for name, L, r, s, e in SYS:
        z, esc, esc3, bx, infl = row(L, r, s, e, frac, c)
        print(f"{lab:<44}{name:<11}{z[0]:>8.2f}{z[1]:>8.2f}{esc:>9.1f}{esc3:>7.1f}{bx:>7.1f}{infl:>6.2f}")
print(f"\nSidak threshold {ZS}; 6.39 excluded after trials at M/L = 4 iff z(M/L=4) > {ZS}.")
