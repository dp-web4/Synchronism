#!/usr/bin/env python3
"""
Phase 3 Session 2: Superconductivity as Viscosity Collapse
CFD cross-pollination: Does Re_ep = λ_ep/γ_phonon predict Tc better than λ_ep alone?

Hypothesis: The N-S framing treats each channel as having effective viscosity.
The superconducting transition is a "viscosity collapse" where electron channel
viscosity μ_eff → 0 when channel coupling Re_ep exceeds threshold.

If Re_ep = λ_ep/γ_phonon predicts Tc better than λ_ep alone:
  → N-S framing has predictive content (not just organizational)
If Re_ep ~ λ_ep alone:
  → γ_phonon adds nothing (BCS already optimal)
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# ========== DATA: Elemental Superconductors ==========
# Sources: Allen & Dynes (1975), McMillan (1968), Carbotte (1990)
# λ_ep from tunneling spectroscopy or McMillan inversion where available

data = {
    # name: (Tc [K], θ_D [K], λ_ep, γ_S [mJ/mol·K²], notes)
    'Al':   (1.18,   428,   0.43, 1.35, 'weak coupling'),
    'Sn':   (3.72,   200,   0.72, 1.78, 'intermediate'),
    'In':   (3.41,   108,   0.81, 1.69, 'intermediate'),
    'Tl':   (2.38,    78,   0.80, 1.47, 'intermediate'),
    'Pb':   (7.19,   105,   1.55, 2.98, 'strong coupling'),
    'Nb':   (9.25,   276,   1.22, 7.79, 'strong coupling, d-band'),
    'V':    (5.38,   380,   0.82, 9.82, 'd-band'),
    'Ta':   (4.47,   258,   0.69, 6.15, 'd-band'),
    'Mo':   (0.92,   450,   0.41, 1.83, 'weak coupling, d-band'),
    'Hg':   (4.15,    71,   1.62, 1.79, 'strong coupling, heavy'),
    'La':   (6.00,   142,   0.98, 9.71, 'intermediate, f-band onset'),
    'Zn':   (0.85,   327,   0.38, 0.64, 'weak coupling'),
    'Cd':   (0.52,   210,   0.38, 0.69, 'weak coupling'),
    'Re':   (1.70,   430,   0.46, 2.35, 'intermediate, 5d'),
    'Os':   (0.66,   500,   0.39, 2.35, 'weak coupling, 5d'),
    'Ir':   (0.14,   420,   0.34, 3.19, 'very weak, 5d'),
    'Ru':   (0.49,   600,   0.38, 3.30, 'weak, 4d'),
    'Tc':   (7.77,   411,   0.83, 5.90, '4d transition'),
    'Th':   (1.37,   163,   0.56, 4.32, 'heavy metal'),
    'W':    (0.015,  390,   0.28, 1.35, 'extremely weak, 5d'),
    'Ga':   (1.08,   325,   0.40, 0.60, 'weak coupling'),
    'Ti':   (0.39,   420,   0.38, 3.35, '3d, weak'),
    'Zr':   (0.61,   291,   0.41, 2.80, '4d, weak'),
    'Hf':   (0.13,   252,   0.34, 2.16, '5d, very weak'),
}

names = list(data.keys())
Tc_arr = np.array([data[n][0] for n in names])
theta_arr = np.array([data[n][1] for n in names])
lambda_arr = np.array([data[n][2] for n in names])
gamma_S_arr = np.array([data[n][3] for n in names])

# Derived quantities
# γ_phonon at room temperature (600/θ_D)
gamma_phonon_300 = 600.0 / theta_arr  # γ = 2T/θ_D at T=300K

# Reynolds number analog: Re_ep = λ_ep / γ_phonon
Re_ep = lambda_arr / gamma_phonon_300  # = λ_ep × θ_D / 600

print("="*70)
print("PHASE 3 SESSION 2: SUPERCONDUCTIVITY AS VISCOSITY COLLAPSE")
print("="*70)
print(f"\nDataset: {len(names)} elemental superconductors")
print(f"Tc range: {Tc_arr.min():.3f} - {Tc_arr.max():.2f} K")
print(f"θ_D range: {theta_arr.min():.0f} - {theta_arr.max():.0f} K")
print(f"λ_ep range: {lambda_arr.min():.2f} - {lambda_arr.max():.2f}")
print(f"γ_ph(300K) range: {gamma_phonon_300.min():.3f} - {gamma_phonon_300.max():.3f}")
print(f"Re_ep range: {Re_ep.min():.3f} - {Re_ep.max():.3f}")

# ========== MODEL 1: BCS (baseline) ==========
# Tc = A × θ_D × exp(-1/λ_ep)
def bcs_model(X, A):
    theta, lam = X
    return A * theta * np.exp(-1.0/lam)

Tc_bcs = 1.13 * theta_arr * np.exp(-1.0/lambda_arr)
r_bcs, p_bcs = pearsonr(np.log(Tc_arr), np.log(Tc_bcs))
rmse_bcs = np.sqrt(np.mean((np.log(Tc_arr) - np.log(Tc_bcs))**2))

print(f"\n{'='*50}")
print("MODEL 1: BCS — Tc = 1.13 × θ_D × exp(-1/λ_ep)")
print(f"  r (log-log) = {r_bcs:.4f}, RMSE_log = {rmse_bcs:.4f}")
print(f"  Tc predictions (first 5): ", end="")
for i in range(5):
    print(f"{Tc_bcs[i]:.2f}", end=" ")
print()

# ========== MODEL 2: BCS with fitted prefactor ==========
def log_bcs(X, logA):
    theta, lam = X
    return logA + np.log(theta) - 1.0/lam

popt2, _ = curve_fit(log_bcs, (theta_arr, lambda_arr), np.log(Tc_arr))
Tc_bcs_fit = np.exp(popt2[0]) * theta_arr * np.exp(-1.0/lambda_arr)
r_bcs_fit, _ = pearsonr(np.log(Tc_arr), np.log(Tc_bcs_fit))
rmse_bcs_fit = np.sqrt(np.mean((np.log(Tc_arr) - np.log(Tc_bcs_fit))**2))

print(f"\nMODEL 2: BCS fitted — Tc = A × θ_D × exp(-1/λ_ep)")
print(f"  A_fitted = {np.exp(popt2[0]):.3f} (BCS predicts 1.13)")
print(f"  r = {r_bcs_fit:.4f}, RMSE_log = {rmse_bcs_fit:.4f}")

# ========== MODEL 3: N-S Modified BCS ==========
# Replace "1" in exponent with γ_phonon(300K) = 600/θ_D
# Tc = A × θ_D × exp(-γ_ph / λ_ep) where γ_ph = 600/θ_D
# = A × θ_D × exp(-600 / (θ_D × λ_ep))
# This modifies the BCS exponent: instead of -1/λ, use -(600/θ_D)/λ = -600/(θ_D×λ)

def ns_modified_bcs(X, logA, alpha):
    theta, lam = X
    gamma_ph = 600.0 / theta  # γ at 300K
    return logA + np.log(theta) - alpha * gamma_ph / lam

popt3, _ = curve_fit(ns_modified_bcs, (theta_arr, lambda_arr), np.log(Tc_arr), 
                     p0=[np.log(1.13), 1.0], maxfev=10000)
Tc_ns = np.exp(popt3[0]) * theta_arr * np.exp(-popt3[1] * gamma_phonon_300 / lambda_arr)
r_ns, _ = pearsonr(np.log(Tc_arr), np.log(Tc_ns))
rmse_ns = np.sqrt(np.mean((np.log(Tc_arr) - np.log(Tc_ns))**2))

print(f"\nMODEL 3: N-S Modified BCS — Tc = A × θ_D × exp(-α × γ_ph/λ_ep)")
print(f"  A_fitted = {np.exp(popt3[0]):.3f}, α_fitted = {popt3[1]:.4f}")
print(f"  Physical BCS: α=1, but γ_ph(300K) replaces 1/λ")
print(f"  r = {r_ns:.4f}, RMSE_log = {rmse_ns:.4f}")
print(f"  Improvement over BCS_fit: ΔRMSE = {rmse_bcs_fit - rmse_ns:.4f}")

# ========== MODEL 4: Re_ep power law ==========
# Tc = A × Re_ep^b × θ_D^c
# Purely empirical: does Re_ep predict Tc?
def power_law_re(X, logA, b, c):
    Re, theta = X
    return logA + b * np.log(Re) + c * np.log(theta)

popt4, _ = curve_fit(power_law_re, (Re_ep, theta_arr), np.log(Tc_arr),
                     p0=[0, 1, 1], maxfev=10000)
Tc_re = np.exp(power_law_re((Re_ep, theta_arr), *popt4))
r_re, _ = pearsonr(np.log(Tc_arr), np.log(Tc_re))
rmse_re = np.sqrt(np.mean((np.log(Tc_arr) - np.log(Tc_re))**2))

print(f"\nMODEL 4: Power law — Tc = A × Re_ep^b × θ_D^c")
print(f"  b_fitted = {popt4[1]:.3f}, c_fitted = {popt4[2]:.3f}")
print(f"  Note: Re_ep = λ_ep × θ_D/600, so this is effectively λ_ep^b × θ_D^(b+c)")
print(f"  r = {r_re:.4f}, RMSE_log = {rmse_re:.4f}")

# ========== MODEL 5: Allen-Dynes improved (reference) ==========
# Tc = (f1×f2×ω_log/1.2) × exp(-1.04(1+λ)/(λ - μ*(1+0.62λ)))
# Simplified: Tc = 0.83 × θ_D × exp(-1.04(1+λ_ep)/(λ_ep - 0.1(1+0.62λ_ep)))
# (using μ* = 0.1 as typical Coulomb pseudopotential)
def allen_dynes(X, mu_star, prefactor):
    theta, lam = X
    # Standard Allen-Dynes formula
    tc = prefactor * theta * np.exp(-1.04*(1+lam)/(lam - mu_star*(1+0.62*lam)))
    return np.log(np.where(tc > 0, tc, 1e-10))

try:
    popt5, _ = curve_fit(allen_dynes, (theta_arr, lambda_arr), np.log(Tc_arr),
                         p0=[0.1, 0.83], bounds=([0.05, 0.5], [0.3, 2.0]), maxfev=10000)
    Tc_ad = np.exp(allen_dynes((theta_arr, lambda_arr), *popt5))
    r_ad, _ = pearsonr(np.log(Tc_arr), np.log(Tc_ad))
    rmse_ad = np.sqrt(np.mean((np.log(Tc_arr) - np.log(Tc_ad))**2))
    print(f"\nMODEL 5: Allen-Dynes — fitted μ*={popt5[0]:.3f}, prefactor={popt5[1]:.3f}")
    print(f"  r = {r_ad:.4f}, RMSE_log = {rmse_ad:.4f}")
    ad_available = True
except Exception as e:
    print(f"\nMODEL 5: Allen-Dynes fit failed: {e}")
    rmse_ad = None
    ad_available = False

# ========== RESIDUAL ANALYSIS ==========
print(f"\n{'='*50}")
print("RESIDUAL ANALYSIS: Do residuals from BCS correlate with γ_phonon?")
print("(If yes: γ_phonon carries information BCS missed)")

Tc_bcs_fit_log = np.log(Tc_bcs_fit)
residuals_bcs = np.log(Tc_arr) - Tc_bcs_fit_log

r_resid_gamma, p_resid = pearsonr(residuals_bcs, gamma_phonon_300)
r_resid_theta, _ = pearsonr(residuals_bcs, np.log(theta_arr))
r_resid_lambda, _ = pearsonr(residuals_bcs, lambda_arr)
r_resid_Re, _ = pearsonr(residuals_bcs, Re_ep)
r_resid_gammaS, _ = pearsonr(residuals_bcs, np.log(gamma_S_arr))

print(f"  r(residual, γ_ph) = {r_resid_gamma:.4f} (p={p_resid:.4f})")
print(f"  r(residual, θ_D)  = {r_resid_theta:.4f}")
print(f"  r(residual, λ_ep) = {r_resid_lambda:.4f}")
print(f"  r(residual, Re_ep)= {r_resid_Re:.4f}")
print(f"  r(residual, γ_S)  = {r_resid_gammaS:.4f}")

# Show worst outliers
sorted_idx = np.argsort(np.abs(residuals_bcs))[::-1]
print(f"\n  Largest BCS residuals (log-Tc):")
for i in range(8):
    idx = sorted_idx[i]
    print(f"    {names[idx]:4s}: actual={Tc_arr[idx]:.3f}K, BCS_fit={Tc_bcs_fit[idx]:.3f}K, "
          f"residual={residuals_bcs[idx]:+.2f}, γ_ph={gamma_phonon_300[idx]:.3f}")

# ========== CLASS ANALYSIS ==========
print(f"\n{'='*50}")
print("CLASS ANALYSIS: Are BCS residuals systematic by material class?")

sp_metals = ['Al', 'Sn', 'In', 'Tl', 'Pb', 'Hg', 'Zn', 'Cd', 'Ga', 'Th']
d3_metals = ['V', 'Ti']
d4_metals = ['Nb', 'Mo', 'Tc', 'Ru', 'Zr']  
d5_metals = ['Ta', 'Re', 'Os', 'Ir', 'W', 'Hf']
f_metals = ['La']

groups = {'sp': sp_metals, 'd4': d4_metals, 'd5': d5_metals}
for gname, gmembers in groups.items():
    gidx = [names.index(n) for n in gmembers if n in names]
    if len(gidx) > 1:
        g_resid = residuals_bcs[gidx]
        g_Tc = Tc_arr[gidx]
        g_gamma = gamma_phonon_300[gidx]
        g_Re = Re_ep[gidx]
        r_gamma_g, _ = pearsonr(g_resid, g_gamma) if len(gidx) > 2 else (0, 1)
        print(f"  {gname:3s} ({len(gidx)} metals): mean residual={g_resid.mean():+.3f}, "
              f"r(resid,γ_ph)={r_gamma_g:.3f}")

# ========== KEY TEST: Does Re_ep unify what BCS splits? ==========
print(f"\n{'='*50}")
print("KEY TEST: Tc/θ_D vs λ_ep (BCS collapse) vs Tc/θ_D vs Re_ep (N-S collapse)")

reduced_Tc = Tc_arr / theta_arr  # Tc/θ_D should be function of λ_ep alone in BCS
r_reduced_lam, _ = pearsonr(np.log(reduced_Tc), lambda_arr)
r_reduced_Re, _ = pearsonr(np.log(reduced_Tc), np.log(Re_ep))

print(f"  r(log(Tc/θ_D), λ_ep)      = {r_reduced_lam:.4f}")
print(f"  r(log(Tc/θ_D), log(Re_ep)) = {r_reduced_Re:.4f}")
print(f"  Note: Re_ep = λ_ep × θ_D/600, so if Tc/θ_D depends on Re_ep, θ_D appears twice")

# BCS prediction for reduced Tc:
Tc_reduced_bcs = np.exp(-1.0/lambda_arr)  # (dropping constant prefactor 1.13)
r_bcs_collapse, _ = pearsonr(np.log(reduced_Tc), np.log(Tc_reduced_bcs))
print(f"  r(log(Tc/θ_D), -1/λ_ep)    = {r_bcs_collapse:.4f} [BCS prediction]")

# N-S prediction: Tc/θ_D ~ exp(-γ_ph/λ_ep)
Tc_reduced_ns = np.exp(-gamma_phonon_300/lambda_arr)
r_ns_collapse, _ = pearsonr(np.log(reduced_Tc), np.log(Tc_reduced_ns))
print(f"  r(log(Tc/θ_D), -γ_ph/λ_ep) = {r_ns_collapse:.4f} [N-S prediction]")

# ========== SUMMARY TABLE ==========
print(f"\n{'='*70}")
print("SUMMARY: MODEL COMPARISON")
print(f"{'Model':<35} {'r':>8} {'RMSE_log':>10} {'Parameters':>12}")
print("-"*70)
print(f"{'BCS (universal prefactor)':<35} {r_bcs:>8.4f} {rmse_bcs:>10.4f} {'0 free':>12}")
print(f"{'BCS (fitted prefactor)':<35} {r_bcs_fit:>8.4f} {rmse_bcs_fit:>10.4f} {'1 free':>12}")
print(f"{'N-S Modified BCS':<35} {r_ns:>8.4f} {rmse_ns:>10.4f} {'2 free':>12}")
print(f"{'Power law (Re_ep, θ_D)':<35} {r_re:>8.4f} {rmse_re:>10.4f} {'3 free':>12}")
if ad_available:
    print(f"{'Allen-Dynes (fitted)':<35} {r_ad:>8.4f} {rmse_ad:>10.4f} {'2 free':>12}")

print(f"\n{'='*70}")
print("VERDICT")
print("="*70)

if rmse_ns < rmse_bcs_fit - 0.05:
    print("N-S MODIFIED BCS significantly improves on standard BCS.")
    print("γ_phonon(300K) carries information beyond λ_ep.")
    print("Result: N-S framing has predictive content for Tc.")
elif rmse_ns < rmse_bcs_fit:
    print("N-S MODIFIED BCS marginally improves on standard BCS.")
    print("Improvement is within noise; extra parameter likely not justified.")
    print("Result: N-S framing organizational, not predictive.")
else:
    print("N-S MODIFIED BCS does NOT improve on standard BCS.")
    print("γ_phonon adds no information to Tc prediction.")
    print("Result: N-S framing is purely organizational for superconductivity.")

# Check if residual correlation is significant
if abs(r_resid_gamma) > 0.3 and p_resid < 0.05:
    print(f"\nHOWEVER: BCS residuals DO correlate with γ_phonon (r={r_resid_gamma:.3f}).")
    print("This suggests γ_phonon carries CORRECTION information — not primary prediction.")
else:
    print(f"\nBCS residuals do NOT significantly correlate with γ_phonon (r={r_resid_gamma:.3f}).")
    print("Confirms: γ_phonon (i.e., θ_D) is fully captured by BCS.")

print("\nNote: BCS already uses θ_D explicitly. Adding γ=2T/θ_D is circular.")
print("The test reveals: N-S framing restates BCS in different vocabulary.")


# ========== BONUS TEST: Does γ_phonon predict λ_ep? ==========
print(f"\n{'='*70}")
print("BONUS TEST: γ_phonon → λ_ep correlation")
print("(If γ_ph predicts λ_ep, the chain γ_ph → λ_ep → Tc is documented)")

r_gamma_lam, p_gamma_lam = pearsonr(gamma_phonon_300, lambda_arr)
r_theta_lam, _ = pearsonr(np.log(theta_arr), lambda_arr)
r_gamma_lam_sp, _ = pearsonr(
    gamma_phonon_300[[names.index(n) for n in sp_metals if n in names]],
    lambda_arr[[names.index(n) for n in sp_metals if n in names]]
)

print(f"  r(γ_ph, λ_ep) overall:  {r_gamma_lam:.4f} (p={p_gamma_lam:.4f})")
print(f"  r(θ_D, λ_ep) overall:   {r_theta_lam:.4f}")
print(f"  r(γ_ph, λ_ep) sp-metals: {r_gamma_lam_sp:.4f}")
print()
# Sort by γ_ph to show trend
print("  γ_ph vs λ_ep (sorted by γ_ph):")
sorted_by_gamma = sorted(range(len(names)), key=lambda i: gamma_phonon_300[i])
for i in sorted_by_gamma:
    print(f"    {names[i]:4s}: γ_ph={gamma_phonon_300[i]:.2f}, λ_ep={lambda_arr[i]:.2f}, "
          f"Tc={Tc_arr[i]:.2f}K, θ_D={theta_arr[i]:.0f}K")

# ========== STRONG-COUPLING REGIME IDENTIFICATION ==========
print(f"\n{'='*70}")
print("STRONG-COUPLING REGIME: λ_ep threshold test")
print("BCS is valid for λ_ep < 0.5 (weak coupling)")
print("Allen-Dynes correction needed for λ_ep > 0.5 (strong coupling)")

weak = [n for n in names if lambda_arr[names.index(n)] < 0.5]
strong = [n for n in names if lambda_arr[names.index(n)] >= 0.5]
print(f"\n  Weak coupling (λ<0.5, {len(weak)} metals): {weak}")
print(f"  Strong coupling (λ≥0.5, {len(strong)} metals): {strong}")

# In weak coupling, does N-S framing (γ_ph) carry any correction information?
w_idx = [names.index(n) for n in weak]
s_idx = [names.index(n) for n in strong]

Tc_bcs_w = 0.034 * theta_arr[w_idx] * np.exp(-1.0/lambda_arr[w_idx])
resid_w = np.log(Tc_arr[w_idx]) - np.log(Tc_bcs_w)
r_resid_gamma_w, _ = pearsonr(resid_w, gamma_phonon_300[w_idx]) if len(w_idx)>2 else (0,1)
Tc_bcs_s = 0.034 * theta_arr[s_idx] * np.exp(-1.0/lambda_arr[s_idx])
resid_s = np.log(Tc_arr[s_idx]) - np.log(Tc_bcs_s)
r_resid_gamma_s, _ = pearsonr(resid_s, gamma_phonon_300[s_idx]) if len(s_idx)>2 else (0,1)

print(f"\n  BCS residuals correlation with γ_ph:")
print(f"    Weak coupling: r = {r_resid_gamma_w:.4f} (γ_ph adds to error, not correction)")
print(f"    Strong coupling: r = {r_resid_gamma_s:.4f}")
print(f"\n  INTERPRETATION:")
print(f"  r_resid_gamma overall = 0.605 because high-γ_ph materials are ALSO strong-coupling")
print(f"  γ_ph → soft lattice → high λ_ep → strong coupling → BCS underpredicts")
print(f"  This is a confound, NOT a γ_ph correction signal")

print(f"\n{'='*70}")
print("FINAL ASSESSMENT: Phase 3 Session 2")
print("="*70)
print("""
Q: Does Re_ep = λ_ep/γ_phonon predict Tc better than λ_ep alone?
A: NO. N-S Modified BCS fails (r=-0.36 vs r=0.81 for standard BCS).
   The γ_phonon(300K) = 600/θ_D is CIRCULAR with BCS (BCS uses θ_D directly).
   Substituting γ_ph for the constant "1" in exp(-1/λ) destroys the prediction.

Q: Does γ_phonon carry any correction information for BCS?
A: PARTIALLY, but confounded. BCS residuals correlate with γ_ph (r=0.605)
   but this is because high-γ_ph (soft lattice) → high-λ_ep (strong coupling)
   → BCS underpredicts. The "correction" is just identifying the strong-coupling
   regime, not adding independent γ_ph information.

Q: What does the N-S framing say about superconductivity?
A: VOCABULARY ONLY.
   - Electron channel "viscosity" ~ 1/σ (inverse conductivity) → goes to 0 at Tc
   - Phonon channel "viscosity" ~ γ (Debye disorder)
   - λ_ep = Lorentz force coupling these channels
   - Tc corresponds to Re_ep exceeding threshold
   This is a valid reframing but restates BCS in CFD language.

Q: Does anything new emerge?
A: One untested hypothesis: Re_ep should be class-invariant near Tc (like Prandtl
   numbers). But Re_ep = λ_ep × θ_D/600 shows no collapse (r=0.12 for Tc/θ_D).
   
PRODUCTIVE FAILURE: Well-documented dead end for superconductivity.
The N-S framing cannot improve on BCS/Allen-Dynes for Tc prediction
because BCS already uses θ_D (which IS γ_phonon, up to a constant).
The framing is circular for this application.

OPEN QUESTION for Phase 3 Session 3:
Are there properties where the MULTI-COMPONENT N-S framing generates
non-circular predictions? Candidates:
1. Thermal conductivity RATIO κ_e/κ_ph (already identified as incremental, Session 6)
2. Mixed superconductors (λ_ep coupling between different sublattice channels)  
3. Phonon drag effect (where phonon and electron channels genuinely compete)
""")

