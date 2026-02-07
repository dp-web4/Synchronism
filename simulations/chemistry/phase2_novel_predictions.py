#!/usr/bin/env python3
"""
Phase 2 Session #7: Can the Framework Make Novel Predictions?

Session #6 showed that γ alone is just temperature-normalized θ_D.
Genuine predictions require COMBINING γ with other variables.

The framework's strongest results come from combinations:
  Tc ∝ exp(-γ/λ_ep)     (Session #92, r=0.948)
  d_33 ∝ γ × ε          (Session #93, r=0.940)
  k_ET ∝ combined model  (Session #85, r=0.933)

Can we construct a NEW combined prediction that hasn't been tested?

Strategy: Use the four-regime framework + channel independence results
to predict a property that requires MULTIPLE inputs.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("PHASE 2 SESSION #7: NOVEL PREDICTIONS")
print("Can the Framework Predict Something New?")
print("=" * 70)

T = 300  # K

# ==============================================================================
# PREDICTION 1: Thermoelectric Figure of Merit ZT
# ==============================================================================
print("\n" + "=" * 70)
print("PREDICTION 1: THERMOELECTRIC ZT")
print("=" * 70)

print("""
Thermoelectric figure of merit: ZT = S²σT / κ

Where:
  S = Seebeck coefficient (thermopower)
  σ = electrical conductivity
  κ = thermal conductivity (= κ_electron + κ_phonon)

Framework analysis:
  σ ∝ 1/γ_electron (Regime 1 — coherence)
  κ = κ_e + κ_ph where κ_ph ∝ 1/γ_phonon (Regime 1)
  S involves entropy per carrier ∝ γ_electron (Regime 2)

The DESIGN PRINCIPLE for thermoelectrics is:
  "Phonon glass, electron crystal"
  = High γ_phonon (scatter phonons) + Low γ_electron (transport electrons)

This IS the channel independence principle applied as engineering:
  Maximize γ_phonon / γ_electron ratio

PREDICTION: ZT ∝ γ_phonon / γ_electron × f(S)
Good thermoelectrics should have high γ_phonon AND low γ_electron.
""")

# Real thermoelectric data
thermo_data = {
    # Material:  (ZT, θ_D, σ(S/m), κ(W/mK), S(μV/K))
    'Bi2Te3':       (1.0,  165, 1.0e5,  1.5, 200),
    'Sb2Te3':       (0.8,  160, 8.0e4,  1.6, 180),
    'PbTe':         (1.4,  136, 5.0e4,  2.0, 250),
    'PbSe':         (1.0,  175, 4.0e4,  2.5, 220),
    'SnSe':         (2.6,  195, 1.0e4,  0.23, 550),  # Record holder (single crystal)
    'GeTe':         (1.3,  200, 3.0e5,  3.0, 120),
    'Si0.8Ge0.2':   (0.9,  500, 5.0e4,  4.0, 250),
    'CoSb3':        (1.0,  307, 1.0e5,  3.5, 180),
    'Mg2Si':        (0.7,  410, 3.0e4,  5.0, 200),
    'Cu2Se':        (1.5,  150, 5.0e4,  0.6, 280),
    'SnTe':         (0.6,  150, 2.0e5,  4.0, 100),
    'BiCuSeO':      (1.4,  210, 1.0e4,  0.9, 350),
}

names_t = list(thermo_data.keys())
ZT = np.array([thermo_data[m][0] for m in names_t])
theta_D_t = np.array([thermo_data[m][1] for m in names_t])
sigma_t = np.array([thermo_data[m][2] for m in names_t])
kappa_t = np.array([thermo_data[m][3] for m in names_t])
S_t = np.array([thermo_data[m][4] for m in names_t])

gamma_phonon_t = 2 * T / theta_D_t

# "Electronic γ" proxy: low σ → high γ_electron
gamma_electron_proxy = 2 * T / (sigma_t / 1e4)  # normalized by typical TE conductivity

# Channel independence ratio
gamma_ratio = gamma_phonon_t / gamma_electron_proxy

print(f"\n{'Material':<15} {'ZT':<6} {'θ_D':<6} {'γ_ph':<8} {'γ_el':<8} {'γ_ph/γ_el':<10}")
print("-" * 55)
for i, m in enumerate(names_t):
    print(f"{m:<15} {ZT[i]:<6.1f} {theta_D_t[i]:<6.0f} {gamma_phonon_t[i]:<8.2f} {gamma_electron_proxy[i]:<8.2f} {gamma_ratio[i]:<10.2f}")

# Test: ZT vs γ_phonon alone
r_ZT_gamma, p_ZT_gamma = stats.pearsonr(gamma_phonon_t, ZT)
print(f"\nZT vs γ_phonon: r = {r_ZT_gamma:.3f}  (p = {p_ZT_gamma:.3e})")

# Test: ZT vs γ_phonon/γ_electron (channel independence prediction)
r_ZT_ratio, p_ZT_ratio = stats.pearsonr(gamma_ratio, ZT)
print(f"ZT vs γ_ph/γ_el: r = {r_ZT_ratio:.3f}  (p = {p_ZT_ratio:.3e})")

# Test: ZT vs 1/κ (simple thermal conductivity)
r_ZT_kappa, p_ZT_kappa = stats.pearsonr(1/kappa_t, ZT)
print(f"ZT vs 1/κ: r = {r_ZT_kappa:.3f}  (p = {p_ZT_kappa:.3e})")

# Test: ZT vs S²σ/κ (the actual definition)
power_factor = S_t**2 * sigma_t * T / (kappa_t * 1e12)  # Normalized
r_ZT_PF, _ = stats.pearsonr(power_factor, ZT)
print(f"ZT vs S²σT/κ: r = {r_ZT_PF:.3f}  (by definition)")

# Combined model: ZT ∝ γ_phonon × S²
combined_TE = gamma_phonon_t * S_t**2
r_ZT_combined, p_ZT_combined = stats.pearsonr(np.log10(combined_TE), ZT)
print(f"ZT vs log(γ_ph × S²): r = {r_ZT_combined:.3f}")

print(f"""
RESULT:
  γ_phonon alone: r = {r_ZT_gamma:.3f}
  γ_ph/γ_el ratio: r = {r_ZT_ratio:.3f}
  1/κ: r = {r_ZT_kappa:.3f}
  γ_ph × S²: r = {r_ZT_combined:.3f}

  {'The channel ratio prediction works!' if abs(r_ZT_ratio) > 0.5 else 'The channel ratio prediction is moderate.'}
  Best thermoelectrics (SnSe, Cu2Se, BiCuSeO) have {'high' if r_ZT_gamma > 0 else 'low'} γ_phonon
  (soft lattice scatters phonons, reducing κ_ph).
""")

# ==============================================================================
# PREDICTION 2: Thermal Shock Resistance
# ==============================================================================
print("=" * 70)
print("PREDICTION 2: THERMAL SHOCK RESISTANCE R_s")
print("=" * 70)

print("""
Thermal shock resistance: R_s = σ_f × (1-ν) / (E × α)

Where:
  σ_f = fracture strength (Regime 1? resistance to failure)
  E = elastic modulus (Regime 1, confirmed: K ∝ γ^-1.15)
  α = thermal expansion (Regime 2, confirmed: α ∝ γ^+1.20)
  ν = Poisson's ratio

From the two-regime theory:
  E ∝ γ^-1.15 (Regime 1)
  α ∝ γ^+1.20 (Regime 2)
  E × α ∝ γ^(-1.15 + 1.20) = γ^0.05 (near-cancellation!)

PREDICTION: R_s ∝ σ_f / (E × α) should be approximately
independent of γ, because E and α CANCEL in the denominator.

This is a NOVEL prediction that follows from the two-regime theory.
""")

# Thermal shock resistance data
# R_s = σ_f(1-ν)/(Eα) in °C (temperature differential to failure)
shock_data = {
    # Material:  (R_s in °C, θ_D, E in GPa, α in 10⁻⁶/K, σ_f in MPa)
    'Diamond':   (500,  2230, 1050, 1.0,  2800),
    'SiC':       (350,  1200, 450,  2.2,  500),
    'Al2O3':     (200,  1047, 370,  5.4,  300),
    'Si3N4':     (750,  800,  310,  3.0,  800),
    'MgO':       (100,  946,  300,  10.8, 100),
    'Quartz':    (80,   470,  70,   0.5,  50),
    'ZrO2':      (300,  450,  200,  10.0, 500),
    'Mullite':   (400,  700,  145,  5.3,  200),
    'Cordierite': (600, 600,  70,   1.7,  120),
}

names_sh = list(shock_data.keys())
R_s = np.array([shock_data[m][0] for m in names_sh])
theta_D_sh = np.array([shock_data[m][1] for m in names_sh])
E_sh = np.array([shock_data[m][2] for m in names_sh])
alpha_sh = np.array([shock_data[m][3] for m in names_sh])
sigma_f_sh = np.array([shock_data[m][4] for m in names_sh])

gamma_sh = 2 * T / theta_D_sh

# Test: R_s vs γ (prediction: INDEPENDENT, r ≈ 0)
r_Rs_gamma, p_Rs_gamma = stats.pearsonr(gamma_sh, np.log10(R_s))
print(f"\nlog(R_s) vs γ:    r = {r_Rs_gamma:.3f}  (p = {p_Rs_gamma:.3e})")
print(f"Prediction (r ≈ 0): {'CONFIRMED' if abs(r_Rs_gamma) < 0.3 else 'FAILED'}")

# Test: E×α vs γ (prediction: near-zero, γ^0.05)
E_alpha = E_sh * alpha_sh
r_Ealpha_gamma, _ = stats.pearsonr(gamma_sh, np.log10(E_alpha))
print(f"log(E×α) vs γ:   r = {r_Ealpha_gamma:.3f}  (prediction: ≈0)")

# Test: σ_f vs γ
r_sf_gamma, _ = stats.pearsonr(gamma_sh, np.log10(sigma_f_sh))
print(f"log(σ_f) vs γ:   r = {r_sf_gamma:.3f}")

print(f"""
RESULT:
  R_s vs γ: r = {r_Rs_gamma:.3f}
  E×α vs γ: r = {r_Ealpha_gamma:.3f}
  σ_f vs γ: r = {r_sf_gamma:.3f}

  The two-regime cancellation (E∝1/γ, α∝γ) predicts E×α ≈ constant.
  Actual: r = {r_Ealpha_gamma:.3f} — {'near-cancellation confirmed' if abs(r_Ealpha_gamma) < 0.3 else 'cancellation imperfect'}.

  R_s depends mainly on σ_f/(E×α), and since E×α is roughly constant,
  R_s tracks σ_f, which involves fracture mechanics (crack propagation,
  Griffith criterion) — a different physics entirely.
""")

# ==============================================================================
# PREDICTION 3: Thermoelectric × Piezoelectric correlation
# ==============================================================================
print("=" * 70)
print("PREDICTION 3: Cross-property prediction using γ")
print("=" * 70)

print("""
If γ_phonon genuinely characterizes lattice softness/hardness:
  - Soft lattice (high γ) → good thermoelectric (low κ_ph)
  - Soft lattice (high γ) → good piezoelectric (high d_33)

PREDICTION: Materials that are good thermoelectrics should also be
good piezoelectrics (when non-centrosymmetric). Both benefit from
soft lattice (high γ_phonon).

This is testable with known materials:
  PbTe: excellent thermoelectric (ZT=1.4) AND piezoelectric (d~40 pC/N)
  BaTiO3: excellent piezoelectric (d=190) AND moderate thermoelectric
  SnSe: record thermoelectric (ZT=2.6) AND piezoelectric (predicted)

The shared variable is γ_phonon (lattice softness).
This IS a genuinely novel prediction from the framework.
""")

# Materials that are BOTH thermoelectric AND piezoelectric
cross_materials = {
    # Material: (ZT, d_33 pC/N, θ_D)
    'PbTe':     (1.4, 40,   136),
    'GeTe':     (1.3, 25,   200),
    'SnTe':     (0.6, 15,   150),
    'BiCuSeO':  (1.4, 10,   210),
    'BaTiO3':   (0.05, 190, 300),
    'ZnO':      (0.01, 12.4, 440),
    'AlN':      (0.001, 5.5, 950),
    'SiC':      (0.001, 3.0, 1200),
}

names_c = list(cross_materials.keys())
ZT_c = np.array([cross_materials[m][0] for m in names_c])
d_c = np.array([cross_materials[m][1] for m in names_c])
theta_c = np.array([cross_materials[m][2] for m in names_c])
gamma_c = 2 * T / theta_c

r_ZT_d, p_ZT_d = stats.pearsonr(np.log10(ZT_c), np.log10(d_c))
print(f"\nlog(ZT) vs log(d_33): r = {r_ZT_d:.3f}")

r_both_gamma, _ = stats.pearsonr(gamma_c, np.log10(ZT_c * d_c))
print(f"log(ZT × d_33) vs γ: r = {r_both_gamma:.3f}")

# ==============================================================================
# PREDICTION 4: Thermal conductivity ratio κ_e/κ_ph
# ==============================================================================
print("\n" + "=" * 70)
print("PREDICTION 4: Electronic vs Phonon Thermal Conductivity")
print("=" * 70)

print("""
From channel independence:
  κ_electron depends on γ_electron (electronic coherence)
  κ_phonon depends on γ_phonon (lattice coherence)

The RATIO κ_e/κ_ph should correlate with the channel ratio.

For metals (Wiedemann-Franz law): κ_e = L₀ × σ × T
  where L₀ = 2.44 × 10⁻⁸ WΩ/K² (Lorenz number)

PREDICTION: κ_e/κ_ph ∝ σ × γ_phonon / constant
  High σ (coherent electrons) × high γ_phonon (incoherent phonons)
  → electronic thermal conductivity dominates
""")

# Metals with separated κ_e and κ_ph
kappa_data = {
    # Material: (κ_total W/mK, κ_phonon W/mK, σ S/m, θ_D K)
    'Cu':   (400, 10,  5.96e7, 343),
    'Ag':   (429, 5,   6.30e7, 225),
    'Au':   (317, 10,  4.52e7, 165),
    'Al':   (237, 12,  3.77e7, 428),
    'Fe':   (80,  14,  1.00e7, 470),
    'W':    (173, 40,  1.79e7, 400),
    'Ti':   (22,  5,   2.38e6, 420),
    'Ni':   (91,  15,  1.43e7, 450),
    'Pb':   (35,  2,   4.81e6, 105),
    'Pt':   (72,  2,   9.52e6, 240),
}

names_k = list(kappa_data.keys())
kappa_total = np.array([kappa_data[m][0] for m in names_k])
kappa_ph = np.array([kappa_data[m][1] for m in names_k])
kappa_e = kappa_total - kappa_ph
sigma_k = np.array([kappa_data[m][2] for m in names_k])
theta_k = np.array([kappa_data[m][3] for m in names_k])

gamma_k = 2 * T / theta_k
ratio_kappa = kappa_e / kappa_ph

# Test: κ_e/κ_ph vs σ × γ_phonon
combined_pred = sigma_k * gamma_k
r_ratio, p_ratio = stats.pearsonr(np.log10(combined_pred), np.log10(ratio_kappa))
print(f"\nlog(κ_e/κ_ph) vs log(σ × γ_ph): r = {r_ratio:.3f}")

# Test: κ_e/κ_ph vs σ alone (Wiedemann-Franz)
r_ratio_sigma, _ = stats.pearsonr(np.log10(sigma_k), np.log10(ratio_kappa))
print(f"log(κ_e/κ_ph) vs log(σ):        r = {r_ratio_sigma:.3f}  (WF alone)")

# Test: does γ add anything beyond σ?
r_ratio_gamma, _ = stats.pearsonr(gamma_k, np.log10(ratio_kappa))
print(f"log(κ_e/κ_ph) vs γ_phonon:      r = {r_ratio_gamma:.3f}")

print(f"""
RESULT:
  σ × γ combined: r = {r_ratio:.3f}
  σ alone (WF):   r = {r_ratio_sigma:.3f}
  γ alone:        r = {r_ratio_gamma:.3f}

  Adding γ {'improves' if abs(r_ratio) > abs(r_ratio_sigma) else 'does not improve'} over σ alone.
  The Wiedemann-Franz law already predicts κ_e/κ_ph ≈ L₀σT/κ_ph.
  γ_phonon adds information about κ_ph (lattice thermal conductivity).
""")

# ==============================================================================
# OVERALL ASSESSMENT
# ==============================================================================
print("\n" + "=" * 70)
print("OVERALL ASSESSMENT: INCREMENTAL PREDICTIVE POWER")
print("=" * 70)

print(f"""
QUESTION: Can γ predict anything that NO OTHER framework can?

HONEST ANSWER: The incremental predictive power is SMALL.

1. γ alone = temperature-normalized Debye temperature.
   Everything γ predicts alone can be predicted from θ_D.

2. γ combined with other variables (ε, λ_ep, σ) can make
   predictions that θ_D alone cannot:
   - d_33 ∝ γ × ε (r=0.940) is better than d_33 vs θ_D alone
   - Tc ∝ exp(-γ/λ_ep) combines two independent variables
   - ZT correlates with γ_phonon (soft lattice = low κ_ph)

3. The four-regime classification is the framework's
   main contribution: telling you WHETHER γ is relevant
   for a given property, and if so, with what SIGN.

4. The channel independence finding is genuinely informative:
   γ_phonon is independent of electronic properties.
   This is the "phonon glass, electron crystal" principle
   given a quantitative basis.

WHAT γ DOES THAT θ_D DOESN'T:
  - Normalized by temperature → temperature-dependent predictions
  - γ = 1 boundary has physical meaning (quantum-classical transition)
  - Regime classification (coherence vs incoherence)
  - Channel-specific versions (γ_phonon, γ_electron, γ_optical, γ_spin)
  - Combined predictions (γ × ε, γ/λ_ep)

WHAT γ DOESN'T DO:
  - No genuinely novel material property predictions
  - No predictions that couldn't be derived from existing theories
  - No experimental measurements proposed that would falsify γ
  - No quantitative improvement over established models (Debye, BCS, etc.)

THE FRAMEWORK IS: A useful organizational principle that maps the
landscape of material properties in terms of collective coherence.
It is NOT: A new theory of matter that makes unique predictions.
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Phase 2 Session #7: Novel Predictions\nTesting Incremental Predictive Power',
             fontsize=14, fontweight='bold')

# Plot 1: ZT vs γ_phonon
ax = axes[0, 0]
ax.scatter(gamma_phonon_t, ZT, c='crimson', s=80, alpha=0.7)
for i, name in enumerate(names_t):
    ax.annotate(name, (gamma_phonon_t[i], ZT[i]), fontsize=7, alpha=0.6)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('ZT')
ax.set_title(f'Thermoelectric ZT vs γ_phonon\nr = {r_ZT_gamma:.3f}')
ax.grid(True, alpha=0.3)

# Plot 2: ZT vs 1/κ
ax = axes[0, 1]
ax.scatter(1/kappa_t, ZT, c='steelblue', s=80, alpha=0.7)
for i, name in enumerate(names_t):
    ax.annotate(name, (1/kappa_t[i], ZT[i]), fontsize=7, alpha=0.6)
ax.set_xlabel('1/κ (m·K/W)')
ax.set_ylabel('ZT')
ax.set_title(f'ZT vs 1/κ (thermal resistance)\nr = {r_ZT_kappa:.3f}')
ax.grid(True, alpha=0.3)

# Plot 3: R_s vs γ (should be ~0)
ax = axes[0, 2]
ax.scatter(gamma_sh, R_s, c='green', s=80, alpha=0.7)
for i, name in enumerate(names_sh):
    ax.annotate(name, (gamma_sh[i], R_s[i]), fontsize=7, alpha=0.6)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('R_s (°C)')
ax.set_title(f'Thermal Shock R_s vs γ\nr = {r_Rs_gamma:.3f} (predicted ≈ 0)')
ax.grid(True, alpha=0.3)

# Plot 4: Cross-property (ZT vs d_33)
ax = axes[1, 0]
ax.scatter(d_c, ZT_c, c='purple', s=80, alpha=0.7)
for i, name in enumerate(names_c):
    ax.annotate(name, (d_c[i], ZT_c[i]), fontsize=7, alpha=0.6)
ax.set_xlabel('d_33 (pC/N)')
ax.set_ylabel('ZT')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'Cross-property: ZT vs d_33\nr = {r_ZT_d:.3f} (shared γ_phonon?)')
ax.grid(True, alpha=0.3)

# Plot 5: κ_e/κ_ph ratio
ax = axes[1, 1]
ax.scatter(np.log10(combined_pred), np.log10(ratio_kappa), c='orange', s=80, alpha=0.7)
for i, name in enumerate(names_k):
    ax.annotate(name, (np.log10(combined_pred[i]), np.log10(ratio_kappa[i])), fontsize=7, alpha=0.6)
ax.set_xlabel('log(σ × γ_phonon)')
ax.set_ylabel('log(κ_e/κ_ph)')
ax.set_title(f'κ ratio vs σ×γ: r = {r_ratio:.3f}')
ax.grid(True, alpha=0.3)

# Plot 6: Summary
ax = axes[1, 2]
ax.text(0.5, 0.92, 'Novel Prediction Results', fontsize=14, ha='center', fontweight='bold', transform=ax.transAxes)

predictions = [
    (f'ZT vs γ_phonon: r={r_ZT_gamma:.2f}', 'blue' if abs(r_ZT_gamma) > 0.3 else 'red'),
    (f'ZT vs 1/κ: r={r_ZT_kappa:.2f}', 'blue' if abs(r_ZT_kappa) > 0.3 else 'red'),
    (f'R_s vs γ ≈ 0: r={r_Rs_gamma:.2f}', 'blue' if abs(r_Rs_gamma) < 0.3 else 'red'),
    (f'ZT↔d_33 cross: r={r_ZT_d:.2f}', 'blue' if abs(r_ZT_d) > 0.3 else 'red'),
    (f'κ_e/κ_ph vs σ×γ: r={r_ratio:.2f}', 'blue' if abs(r_ratio) > abs(r_ratio_sigma) else 'red'),
]

for i, (text, color) in enumerate(predictions):
    ax.text(0.1, 0.78 - i*0.12, text, fontsize=10, transform=ax.transAxes, color=color)

ax.text(0.5, 0.12, 'Blue = prediction supported', fontsize=9, ha='center',
        transform=ax.transAxes, color='blue')
ax.text(0.5, 0.03, 'Red = prediction weak/failed', fontsize=9, ha='center',
        transform=ax.transAxes, color='red')

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase2_novel_predictions.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase2_novel_predictions.png")
