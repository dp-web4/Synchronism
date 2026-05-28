#!/usr/bin/env python3
"""
Session 678: Verify the cluster-bridge impossibility of C(rho) on Coma.

The site-explorer proposal (c_rho_cluster_bridge_impossibility_coma.md) computed
4 natural ansaetze on a Coma beta-model and found all miss observed
M_lens/M_B ~ 4.6 by orders of magnitude under galaxy-anchored rho_crit. I verify
the load-bearing math here (S672/S675 discipline: compute, don't accept the
framing). The decisive piece is structural: A3 has a CODOMAIN UPPER BOUND of 2
for any gamma, any rho_crit, any cluster -- it CAN'T reach 4.6.
"""
import numpy as np

# --- physical constants ----------------------------------------------------------
mp   = 1.6726e-24        # g
mu   = 0.61              # mean molecular weight, fully ionized gas
kpc  = 3.086e21          # cm

# --- Coma beta-model gas (Briel+ 1992) ------------------------------------------
n0_cm3 = 3.4e-3
rho0   = n0_cm3 * mu * mp     # central gas mass density, g/cm^3
rc_cm  = 290 * kpc            # core radius
beta   = 0.65
def rho_gas(r_cm):
    return rho0 * (1.0 + (r_cm/rc_cm)**2)**(-1.5*beta)

# --- framework: C(rho) and gamma -------------------------------------------------
def C(rho_g, rho_crit, gamma):
    # C = tanh(gamma * ln(rho/rho_crit + 1))
    return np.tanh(gamma * np.log(rho_g/rho_crit + 1.0))

# galaxy-anchored rho_crit: a representative galactic mass density,
# 10^-24 to 10^-22 g/cm^3 spans MW DM halo to disk densities; scan a range.
gamma = 2.0                                # single-particle reference (N_corr=1)
r_grid_cm = np.linspace(1*kpc, 1300*kpc, 5001)
r_500_cm  = 1300 * kpc                     # Coma r_500 ~ 1.3 Mpc
rho_r = rho_gas(r_grid_cm)
dV    = 4*np.pi * r_grid_cm**2 * np.gradient(r_grid_cm)
M_B   = np.sum(rho_r * dV)

# --- (A) A3 codomain bound -- the structural impossibility -----------------------
print("=" * 74)
print("(A) A3 codomain bound: M_app/M_B = 1 + integ C*rho dV / M_B")
print("    C in [0,1) -> integ C*rho dV <= M_B  =>  M_app/M_B <= 2, exactly.")
print("=" * 74)
print(f"   Observed Coma M_lens/M_B at r_500: 4.6")
print(f"   A3 upper bound (any gamma, any rho_crit, any cluster): 2.0")
print(f"   => A3 cannot reach 4.6. STRUCTURAL impossibility, codomain-tight.")
# verify numerically by saturating C to 1 (its supremum)
A3_max = 1.0 + np.sum(1.0 * rho_r * dV) / M_B
print(f"   numerical: A3 with C-->1 everywhere = 1 + M_B/M_B = {A3_max:.3f}")
print(f"   => Coma's 4.6 is structurally outside the codomain of A3.")

# --- (B) A2/A3 collapse to Newtonian under galaxy-anchored rho_crit --------------
print()
print("=" * 74)
print("(B) Cluster densities under galaxy-anchored rho_crit -> C ~ 0 -> A2/A3 ~ 1")
print("=" * 74)
print(f"   rho_gas(Coma center) = {rho_gas(0):.2e} g/cm^3")
print(f"   rho_gas(r_500)       = {rho_gas(r_500_cm):.2e} g/cm^3")
print()
print(f"   {'rho_crit (g/cm^3)':<18} {'<C>_vol':>12} {'<C>_mass':>12} {'A2=1/(1-<C>)':>14} {'A3=1+<Cρ>/M_B':>15}")
for rho_crit in [1e-24, 1e-23, 1e-22, 1e-21]:
    Cr = C(rho_r, rho_crit, gamma)
    Cmean_vol  = np.sum(Cr * dV) / np.sum(dV)
    Cmean_mass = np.sum(Cr * rho_r * dV) / M_B
    A2 = 1.0/(1.0 - Cmean_vol)
    A3 = 1.0 + np.sum(Cr * rho_r * dV) / M_B
    print(f"   {rho_crit:<18.0e} {Cmean_vol:>12.2e} {Cmean_mass:>12.2e} {A2:>14.4f} {A3:>15.4f}")
print("   => galaxy-anchored rho_crit (1e-24..1e-22 g/cm^3) gives A2,A3 ~ 1 (Newtonian).")
print("      Cluster needs 4.6; A2/A3 undershoot by factor ~5. Confirms proposal.")

# --- (C) A1/A4 overshoot: 1/C diverges as C -> 0 ---------------------------------
print()
print("=" * 74)
print("(C) A1=<1/C>_vol and A4=1/<C>_mass diverge as C -> 0  (cluster overshoot)")
print("=" * 74)
print(f"   {'rho_crit':<18} {'<1/C>_vol':>14} {'1/<C>_mass':>14}")
for rho_crit in [1e-24, 1e-23, 1e-22, 1e-21]:
    Cr = C(rho_r, rho_crit, gamma)
    Cr_safe = np.maximum(Cr, 1e-30)
    A1 = np.sum((1.0/Cr_safe) * dV) / np.sum(dV)
    A4 = 1.0/(np.sum(Cr * rho_r * dV) / M_B)
    print(f"   {rho_crit:<18.0e} {A1:>14.2e} {A4:>14.2e}")
print("   => A1 and A4 overshoot by 10^3 - 10^5. The same C ~ 0 that pins A2/A3 to")
print("      ~1 makes 1/C huge: there is no rho_crit that hits 4.6 with any ansatz.")

# --- (D) The structural root: one density scale -----------------------------------
print()
print("=" * 74)
print("(D) One-density-scale insufficiency (structural, no computation needed)")
print("=" * 74)
print("   C(rho) = tanh(gamma * ln(rho/rho_crit + 1)) has exactly ONE scale, rho_crit.")
print("   tanh saturates: rho >> rho_crit -> C ~ 1; rho << rho_crit -> C ~ 0.")
print("   To stay in the active transition at BOTH galaxy and cluster densities")
print("   (which differ by ~10^4), you would need rho_crit to sit between them --")
print("   but then C(rho_galaxy) saturates to 1 and gives no galactic effect either.")
print("   The transition has only one regime of action. Verlinde (a0 + r), MOND")
print("   (a0 + a_N) have TWO scales -- the second is what enables bounded")
print("   enhancements at multiple regimes. C(rho) lacks it by construction.")

# --- (E) verdict ------------------------------------------------------------------
print()
print("=" * 74)
print("VERDICT")
print("=" * 74)
print("   - A3 is impossible by CODOMAIN BOUND (<= 2; Coma needs 4.6). Tight, exact,")
print("     independent of gamma, rho_crit, cluster choice.")
print("   - A2 collapses to Newtonian under galaxy-anchored rho_crit (C ~ 0).")
print("   - A1/A4 overshoot by 10^3 - 10^5 (1/C diverges as C -> 0).")
print("   - One-density-scale insufficiency is the structural root: a tanh with one")
print("     scale cannot bridge two density regimes separated by 10^4.")
print("   - Companion to S665 (one scalar field, no vorticity), S667 (one time order,")
print("     no causality + dissipation + oscillation), S677 (no decoherence param):")
print("     C(rho) lacks the structural ingredients (extra scale, kinetic term,")
print("     second time derivative) that real modified-gravity / wave theories use.")
print("   - Confirms the proposal: the modified-gravity track's 'cluster bridge' is")
print("     structurally closed; the C(a) -> C(rho) migration silently cost it.")
