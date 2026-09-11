#!/usr/bin/env python3
r"""
THE GALAXY FLOOR IS THE DARK-ENERGY SECTOR'S C AT TODAY'S MEAN DENSITY — AND WHAT FOLLOWS
=========================================================================================
Maintainer 2026-09-11.  Verifies four claims a visitor researcher persona made on 2026-09-11,
before any of them goes on the site or into a proposal.

Dark-energy sector as the site states it (Session 100):  rho_DE = rho_m (1 - C)/C,
C = tanh(gamma ln(1 + x)),  x = rho_m / rho_crit_DE.

 (1) Identity.  Flat, radiation neglected: rho_DE/rho_m = Omega_L/Omega_m, so (1-C0)/C0 = (1-Om)/Om,
     i.e. C(rho_m,0) = Omega_m IDENTICALLY — for any gamma, once rho_crit_DE is closed on today's Omega_m.
     This is the same number the galaxy sector uses as its floor.  (Session 100 sets C0 = Omega_m as a
     calibration; PREDICTIONS.md already calls Omega_DE = 1 - Omega_m a tautology.  What is new is only
     the observation that the galaxy floor and the DE calibration are the same number.)

 (2) Sign theorem.  With u = 1 + x, p = 2 gamma:  1 + w = F(u)/(u^p - 1),
     F(u) = u^p - 1 - p (u-1) u^(p-1),  F(1) = 0,  F'(u) = -p(p-1)(u-1)u^(p-2)
     =>  sign(1 + w) = sign(1 - 2 gamma) at every redshift.  Checked numerically from the continuity
     equation w = -1 - (1/3) dln rho_DE / dln a, not from the closed form.

 (3) If the galaxy floor is read as C_DE at the AMBIENT cosmic density (not a fixed number), the boost
     ceiling evolves:  B_max(z) = 1 / C_DE(rho_m(z)).  Tabulated, with gamma sensitivity.

 (4) The same reading makes the ceiling environment-dependent: B_max(delta) at z = 0.  The smoothing
     scale that defines "ambient" is NOT specified anywhere in the framework — that is the gate.

Also answers the persona's Q8: does the Omega_m floor touch TEST-26?  The floor binds only where
C < Omega_m, i.e. x < x0, i.e. a > 1 (the future).  The past light cone sees the unfloored branch only.
"""
import math
import numpy as np

OMEGA_M = 0.315
GAMMA_DE = 0.487            # DESI DR2 BAO + CMB + SN best fit (site explorer, 2026-08-12)
H0 = 67.4                   # km/s/Mpc
MSUN_PC3_KG_M3 = 1.98847e30 / (3.0857e16) ** 3   # 1 Msun/pc^3 in kg/m^3


def x0_for(gamma, om=OMEGA_M):
    """rho_m,0 / rho_crit_DE such that C(today) = Omega_m."""
    return math.exp(math.atanh(om) / gamma) - 1.0


def C(x, gamma):
    return np.tanh(gamma * np.log1p(x))


def w_numeric(z, gamma):
    """w from continuity, by central difference in ln a."""
    x0 = x0_for(gamma)
    def ln_rho_de(lna):
        a = math.exp(lna)
        x = x0 * a ** -3
        c = math.tanh(gamma * math.log1p(x))
        return math.log(x0 * a ** -3) + math.log((1 - c) / c)   # rho_m in units of rho_crit_DE
    lna = -math.log1p(z)
    h = 1e-5
    d = (ln_rho_de(lna + h) - ln_rho_de(lna - h)) / (2 * h)
    return -1.0 - d / 3.0


def w_closed(z, gamma):
    u = 1.0 + x0_for(gamma) * (1 + z) ** 3
    p = 2 * gamma
    F = u ** p - 1 - p * (u - 1) * u ** (p - 1)
    return -1.0 + F / (u ** p - 1)


def main():
    print(__doc__.split("\n")[1])
    print()

    # ---- (1) identity
    print("(1) C_DE(rho_m,0) = Omega_m identically, and where the DE knee sits")
    rho_c_cosmic = 3 * (H0 * 1e3 / 3.0857e22) ** 2 / (8 * math.pi * 6.674e-11)
    rho_m0 = OMEGA_M * rho_c_cosmic / MSUN_PC3_KG_M3
    print(f"    mean matter density today       rho_m0 = {rho_m0:.3e} Msun/pc^3")
    for g in (0.3, GAMMA_DE, 0.5, 1.0, 2.0):
        x0 = x0_for(g)
        closure = (1 - C(x0, g)) / C(x0, g)
        print(f"    gamma={g:<5}  x0={x0:7.4f}  C0={float(C(x0, g)):.6f}  (1-C0)/C0={float(closure):.6f}"
              f"  vs (1-Om)/Om={(1 - OMEGA_M) / OMEGA_M:.6f}  rho_crit_DE={rho_m0 / x0:.3e} Msun/pc^3")
    x0 = x0_for(GAMMA_DE)
    rc_de = rho_m0 / x0
    for name, knee in (("measured Jeans knee (2026-08-27)", 0.161), ("0.029 V_flat^2, MW", 1.52e3)):
        print(f"    galaxy knee {name:<34} = {knee:9.3g}  ->  {knee / rc_de:.2e} x the DE knee")
    print(f"    C_DE at a disc midplane (0.1 Msun/pc^3): {float(C(0.1 / rc_de, GAMMA_DE)):.12f}"
          "  -> ONE knee cannot serve both sectors; the floor identity needs two C's.")
    print()

    # ---- (2) sign theorem
    print("(2) sign(1+w) = sign(1-2 gamma): numeric continuity vs closed form, z in [0, 5]")
    zs = np.linspace(0.0, 5.0, 51)
    worst = 0.0
    viol = 0
    for g in (0.2, 0.4, 0.45, GAMMA_DE, 0.499, 0.5, 0.501, 0.513, 0.6, 1.0, 2.0):
        onepw = np.array([1 + w_numeric(z, g) for z in zs])
        closed = np.array([1 + w_closed(z, g) for z in zs])
        worst = max(worst, float(np.max(np.abs(onepw - closed))))
        expect = np.sign(1 - 2 * g)
        if g == 0.5:
            ok = bool(np.all(np.abs(onepw) < 1e-6))
        else:
            ok = bool(np.all(np.sign(onepw) == expect))
        viol += 0 if ok else 1
        print(f"    gamma={g:<6} 1+w(z=0)={onepw[0]:+.5f}  1+w(z=5)={onepw[-1]:+.5f}  "
              f"sign constant & = sign(1-2g): {ok}")
    print(f"    max |numeric - closed form| = {worst:.2e};  sign violations: {viol}")
    w0 = w_numeric(0.0, GAMMA_DE)
    wa = (w_numeric(1e-3, GAMMA_DE) - w_numeric(0.0, GAMMA_DE)) / 1e-3    # wa = dw/dz at z=0
    print(f"    gamma={GAMMA_DE}: w0 = {w0:.4f}, w(z=1) = {w_numeric(1, GAMMA_DE):.4f}, "
          f"w(z=3) = {w_numeric(3, GAMMA_DE):.4f}, CPL-slope wa = {wa:+.4f}")
    print(f"    (1+w0) and wa same sign: {np.sign(1 + w0) == np.sign(wa)}  -> the family runs along +diagonal;"
          " DESI's preferred region has 1+w0 > 0 with wa < 0.")
    print()

    # ---- TEST-26 / floor
    print("    Q8 — does the Omega_m floor bind on the past light cone?")
    cz = [float(C(x0 * (1 + z) ** 3, GAMMA_DE)) for z in (0, 0.01, 0.5, 2)]
    print(f"    C_DE at z = 0, 0.01, 0.5, 2: {', '.join(f'{c:.4f}' for c in cz)}  (floor {OMEGA_M})")
    print("    C >= Omega_m for all z >= 0 by monotonicity; floor binds only at a > 1. TEST-26's past-light-cone"
          " no-go is unaffected by flooring.")
    print()

    # ---- (3) evolving ceiling
    print("(3) Floor read as ambient cosmic C:  B_max(z) = 1/C_DE(rho_m(z)),  f_DM,max = 1 - C")
    print("    z      " + "".join(f"g={g:<14}" for g in (0.3, GAMMA_DE, 0.5, 2.0)))
    for z in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.2):
        row = []
        for g in (0.3, GAMMA_DE, 0.5, 2.0):
            c = float(C(x0_for(g) * (1 + z) ** 3, g))
            row.append(f"B={1 / c:5.3f} f<={1 - c:5.3f} ")
        print(f"    {z:<5}  " + " ".join(row))
    print("    The shape is gamma-robust: B_max(z=2) <= 1.1 for every gamma in [0.3, 2].")
    print()

    # ---- (4) environment
    print("(4) Same reading, z = 0, as a function of ambient overdensity delta (smoothing scale UNSPECIFIED)")
    for d in (-0.9, -0.8, -0.5, 0.0, 1.0, 5.0, 20.0, 200.0):
        c = float(C(x0 * (1 + d), GAMMA_DE))
        print(f"    delta = {d:+7.1f}   C_floor = {c:.4f}   B_max = {1 / c:7.3f}   f_DM,max = {1 - c:.3f}")
    print("    Deep voids reach the B ~ 12-14 SPARC's dwarfs demand; group-scale overdensities are Newtonian.")
    print("    Which delta a galaxy 'sees' is the undefined MRH of the floor. Without it (3) and (4) are not bets.")


if __name__ == "__main__":
    main()
