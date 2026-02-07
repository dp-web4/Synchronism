#!/usr/bin/env python3
"""
Phase 2 Session #6: Testing the Four-Regime Framework

The framework predicts which regime each property belongs to:
  Regime 0 (Neutral): counting → γ irrelevant
  Regime 1 (Coherence): propagation → P ∝ 1/γ
  Regime 2 (Incoherence): response → P ∝ γ
  Regime 3 (Barrier): activated → P ∝ exp(-E/kT)

Test with properties NOT previously analyzed in the framework.
This is the "stop cataloguing, start testing" recommendation.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("PHASE 2 SESSION #6: TESTING FOUR-REGIME PREDICTIONS")
print("New Properties Not Previously Analyzed")
print("=" * 70)

T = 300  # K

# ==============================================================================
# TEST 1: Speed of Sound (Predicted: Regime 1 — Coherence)
# ==============================================================================
print("\n" + "=" * 70)
print("TEST 1: SPEED OF SOUND")
print("Prediction: Regime 1 (propagation → v_s ∝ 1/γ or v_s ∝ γ^(-n))")
print("=" * 70)

# Speed of sound (m/s) and Debye temperature (K) for various materials
sound_data = {
    'Diamond':  (12000, 2230),
    'SiC':      (9800,  1200),
    'Al2O3':    (10000, 1047),
    'Si':       (8433,  640),
    'Fe':       (5120,  470),
    'Cu':       (3810,  343),
    'Al':       (6320,  428),
    'W':        (5220,  400),
    'Ag':       (2600,  225),
    'Au':       (2030,  165),
    'Pb':       (1260,  105),
    'Na':       (3200,  158),
    'K':        (2000,  91),
    'Ti':       (6070,  420),
    'Ni':       (4970,  450),
    'Mg':       (5770,  400),
    'Zn':       (3850,  327),
    'NaCl':     (3530,  321),
    'KCl':      (3150,  235),
    'MgO':      (9200,  946),
}

names_s = list(sound_data.keys())
v_s = np.array([sound_data[m][0] for m in names_s])
theta_D_s = np.array([sound_data[m][1] for m in names_s])
gamma_s = 2 * T / theta_D_s

log_v = np.log10(v_s)
log_g = np.log10(gamma_s)

r_v_gamma, p_v_gamma = stats.pearsonr(gamma_s, log_v)
slope_v, intercept_v = np.polyfit(log_g, log_v, 1)
r_v_ll, p_v_ll = stats.pearsonr(log_g, log_v)

print(f"\nlog(v_s) vs γ:      r = {r_v_gamma:.3f}  (p = {p_v_gamma:.3e})")
print(f"log(v_s) vs log(γ): r = {r_v_ll:.3f}  → v_s ∝ γ^{slope_v:.2f}")
print(f"Predicted regime: {'CONFIRMED (Regime 1)' if slope_v < -0.3 else 'UNEXPECTED'}")

# Note: v_s should be related to θ_D by v_D = k_B θ_D / (ℏ(6π²n)^(1/3))
# so v_s ∝ θ_D ∝ 1/γ is expected from the DEFINITION of Debye temperature
# This is partially tautological! Flag this.
print(f"""
CAVEAT: v_s ∝ θ_D is partially tautological because the Debye temperature
is DEFINED from the speed of sound: θ_D = (ℏ/k_B)(6π²n)^(1/3) × v_Debye.
So γ = 2T/θ_D = 2T k_B / (ℏ(6π²n)^(1/3) × v_D) inherently contains v_s.
This test confirms internal consistency but does NOT test predictive power.
""")

# ==============================================================================
# TEST 2: Hardness (Predicted: Regime 1 — Coherence)
# ==============================================================================
print("=" * 70)
print("TEST 2: VICKERS HARDNESS")
print("Prediction: Regime 1 (resistance to deformation → H ∝ 1/γ)")
print("=" * 70)

# Vickers hardness (GPa) and Debye temperature
hardness_data = {
    'Diamond':  (115.0, 2230),
    'SiC':      (28.0,  1200),
    'Al2O3':    (20.0,  1047),
    'Si':       (11.0,  640),
    'MgO':      (8.0,   946),
    'Fe':       (1.5,   470),
    'Cu':       (0.37,  343),
    'Al':       (0.17,  428),
    'W':        (3.4,   400),
    'Ag':       (0.25,  225),
    'Au':       (0.22,  165),
    'Pb':       (0.04,  105),
    'Na':       (0.007, 158),
    'Ti':       (1.0,   420),
    'Ni':       (0.64,  450),
    'Zn':       (0.39,  327),
    'Cr':       (1.1,   630),
}

names_h = list(hardness_data.keys())
H_v = np.array([hardness_data[m][0] for m in names_h])
theta_D_h = np.array([hardness_data[m][1] for m in names_h])
gamma_h = 2 * T / theta_D_h

log_H = np.log10(H_v)
log_g_h = np.log10(gamma_h)

r_H_gamma, p_H_gamma = stats.pearsonr(gamma_h, log_H)
slope_H, intercept_H = np.polyfit(log_g_h, log_H, 1)
r_H_ll, p_H_ll = stats.pearsonr(log_g_h, log_H)

print(f"\nlog(H) vs γ:      r = {r_H_gamma:.3f}  (p = {p_H_gamma:.3e})")
print(f"log(H) vs log(γ): r = {r_H_ll:.3f}  → H ∝ γ^{slope_H:.2f}")
print(f"Predicted regime: {'CONFIRMED (Regime 1)' if slope_H < -0.3 else 'UNEXPECTED'}")

# Hardness is famously NOT just about lattice stiffness — it involves
# dislocation mobility, which is a complex function of crystal structure
# Let's check if there's a material-class effect
print(f"""
NOTE: Hardness involves dislocation mobility, grain boundaries, and
crystal structure — not just lattice stiffness. The correlation with
θ_D may be partially spurious (stiff lattice → high θ_D AND hard).
Ceramics (Diamond, SiC, Al2O3) dominate both θ_D and H ranges.
""")

# Test: remove ceramics
metal_mask = np.array([m not in ['Diamond', 'SiC', 'Al2O3', 'MgO'] for m in names_h])
if np.sum(metal_mask) >= 4:
    r_H_metal, p_H_metal = stats.pearsonr(gamma_h[metal_mask], log_H[metal_mask])
    slope_H_metal, _ = np.polyfit(np.log10(gamma_h[metal_mask]), log_H[metal_mask], 1)
    print(f"Metals only (N={np.sum(metal_mask)}): r = {r_H_metal:.3f}, H ∝ γ^{slope_H_metal:.2f}")

# ==============================================================================
# TEST 3: Ductility/Elongation (Predicted: Regime 2 — Incoherence)
# ==============================================================================
print("\n" + "=" * 70)
print("TEST 3: DUCTILITY (% Elongation at Break)")
print("Prediction: Regime 2 (ease of deformation → ductility ∝ γ)")
print("=" * 70)

# % elongation at break and Debye temperature for pure metals
ductility_data = {
    'Au':  (45, 165),    # Very ductile
    'Ag':  (45, 225),    # Very ductile
    'Cu':  (42, 343),    # Ductile
    'Al':  (40, 428),    # Ductile
    'Pb':  (50, 105),    # Very ductile (soft)
    'Na':  (70, 158),    # Extremely ductile
    'Fe':  (25, 470),    # Moderate
    'Ni':  (30, 450),    # Moderate
    'Ti':  (20, 420),    # Low ductility
    'W':   (2,  400),    # Brittle (at RT)
    'Cr':  (0.5, 630),   # Very brittle
    'Mg':  (8,  400),    # Low ductility (hcp)
    'Zn':  (10, 327),    # Low ductility (hcp)
}

names_d = list(ductility_data.keys())
elong = np.array([ductility_data[m][0] for m in names_d])
theta_D_d = np.array([ductility_data[m][1] for m in names_d])
gamma_d = 2 * T / theta_D_d

r_elong_gamma, p_elong_gamma = stats.pearsonr(gamma_d, elong)
r_elong_invg, _ = stats.pearsonr(1/gamma_d, elong)

print(f"\nElongation vs γ:   r = {r_elong_gamma:.3f}  (p = {p_elong_gamma:.3e})")
print(f"Elongation vs 1/γ: r = {r_elong_invg:.3f}")
print(f"Predicted regime: {'CONFIRMED (Regime 2)' if r_elong_gamma > 0.3 else 'Regime 1' if r_elong_invg > 0.3 else 'UNCERTAIN'}")

print(f"""
NOTE: Ductility depends heavily on crystal structure (fcc >> bcc > hcp),
not just lattice softness. Au, Cu, Ag (fcc, ductile) vs Cr, W (bcc, brittle)
vs Mg, Zn (hcp, poor ductility). The prediction is weak because crystal
structure, not coherence, is the primary determinant.
""")

# ==============================================================================
# TEST 4: Melting Point (Predicted: Regime 1 — Coherence/Stability)
# ==============================================================================
print("=" * 70)
print("TEST 4: MELTING POINT")
print("Prediction: Regime 1 (stability → T_m ∝ 1/γ)")
print("=" * 70)

# Melting point (K) and Debye temperature
melt_data = {
    'W':      (3695, 400),
    'Re':     (3459, 430),
    'Mo':     (2896, 450),
    'Ta':     (3290, 258),
    'Fe':     (1811, 470),
    'Ni':     (1728, 450),
    'Cu':     (1358, 343),
    'Al':     (933,  428),
    'Ag':     (1235, 225),
    'Au':     (1337, 165),
    'Pb':     (601,  105),
    'Na':     (371,  158),
    'K':      (337,  91),
    'Cs':     (302,  38),
    'Ti':     (1941, 420),
    'Cr':     (2180, 630),
    'Si':     (1687, 640),
    'Ge':     (1211, 374),
    'Zn':     (693,  327),
    'Mg':     (923,  400),
}

names_m = list(melt_data.keys())
T_m = np.array([melt_data[m][0] for m in names_m])
theta_D_m = np.array([melt_data[m][1] for m in names_m])
gamma_m = 2 * T / theta_D_m

log_Tm = np.log10(T_m)

r_Tm_gamma, p_Tm_gamma = stats.pearsonr(gamma_m, log_Tm)
r_Tm_theta, p_Tm_theta = stats.pearsonr(theta_D_m, log_Tm)

print(f"\nlog(T_m) vs γ:    r = {r_Tm_gamma:.3f}  (p = {p_Tm_gamma:.3e})")
print(f"log(T_m) vs θ_D:  r = {r_Tm_theta:.3f}  (p = {p_Tm_theta:.3e})")

# The Lindemann criterion: T_m ∝ θ_D² × M × a² / (C × ℏ²)
# So T_m should correlate with θ_D² not θ_D
r_Tm_theta2, _ = stats.pearsonr(theta_D_m**2, T_m)
print(f"T_m vs θ_D²:      r = {r_Tm_theta2:.3f}  (Lindemann criterion)")

# γ at melting point: γ_m = 2T_m/θ_D
gamma_at_melt = 2 * T_m / theta_D_m
print(f"\nγ at melting point: mean = {gamma_at_melt.mean():.2f}, std = {gamma_at_melt.std():.2f}")
print(f"  Range: {gamma_at_melt.min():.2f} - {gamma_at_melt.max():.2f}")

# Is γ_melt approximately constant? That would be the Lindemann criterion
# in coherence framework language
print(f"""
KEY INSIGHT: The Lindemann melting criterion states that melting occurs
when atomic displacements reach ~10% of nearest-neighbor distance.

In the coherence framework, this translates to:
  γ_melt = 2T_m/θ_D ≈ constant

Our data: γ_melt = {gamma_at_melt.mean():.2f} ± {gamma_at_melt.std():.2f}

This is {'approximately constant' if gamma_at_melt.std()/gamma_at_melt.mean() < 0.4 else 'NOT constant'} (CV = {gamma_at_melt.std()/gamma_at_melt.mean():.2f}).

The Lindemann criterion IS the coherence framework for melting:
"A solid melts when its phonon coherence parameter exceeds a
universal threshold γ_melt ≈ {gamma_at_melt.mean():.1f}."

This is NOT a prediction of the framework — it's a REDISCOVERY
of Lindemann (1910) in γ language. It confirms internal consistency
but adds nothing new.
""")

# ==============================================================================
# TEST 5: Dielectric Loss tan(δ) (Predicted: Regime 2 — Incoherence)
# ==============================================================================
print("=" * 70)
print("TEST 5: DIELECTRIC LOSS tan(δ)")
print("Prediction: Regime 2 (energy absorption → tan(δ) ∝ γ)")
print("=" * 70)

# Dielectric loss at 1 MHz and Debye temperature for ceramics/insulators
dielectric_loss_data = {
    'Diamond':    (1e-5,   2230),   # Ultra-low loss
    'Sapphire':   (1e-5,   1047),   # Ultra-low loss
    'Quartz':     (1e-4,   470),    # Very low loss
    'MgO':        (5e-4,   946),    # Low loss
    'Si3N4':      (3e-3,   800),    # Moderate
    'AlN':        (2e-3,   950),    # Low loss
    'BN':         (5e-4,   1900),   # Low loss
    'SiO2_glass': (1e-3,   300),    # Amorphous → higher loss
    'BaTiO3':     (0.05,   300),    # High loss (ferroelectric)
    'PZT':        (0.02,   240),    # High loss (ferroelectric)
    'TiO2':       (5e-3,   700),    # Moderate
    'ZrO2':       (0.01,   450),    # Moderate-high
}

names_dl = list(dielectric_loss_data.keys())
tan_delta = np.array([dielectric_loss_data[m][0] for m in names_dl])
theta_D_dl = np.array([dielectric_loss_data[m][1] for m in names_dl])
gamma_dl = 2 * T / theta_D_dl

log_tan = np.log10(tan_delta)

r_tan_gamma, p_tan_gamma = stats.pearsonr(gamma_dl, log_tan)
r_tan_ll, _ = stats.pearsonr(np.log10(gamma_dl), log_tan)
slope_tan, _ = np.polyfit(np.log10(gamma_dl), log_tan, 1)

print(f"\nlog(tan δ) vs γ:      r = {r_tan_gamma:.3f}  (p = {p_tan_gamma:.3e})")
print(f"log(tan δ) vs log(γ): r = {r_tan_ll:.3f}  → tan δ ∝ γ^{slope_tan:.2f}")
print(f"Predicted regime: {'CONFIRMED (Regime 2)' if r_tan_gamma > 0.3 else 'UNCERTAIN'}")

# ==============================================================================
# TEST 6: Creep Rate (Predicted: Regime 2+3 hybrid)
# ==============================================================================
print("\n" + "=" * 70)
print("TEST 6: CREEP RATE")
print("Prediction: Regime 2+3 (deformation + activation → complex)")
print("=" * 70)

# Steady-state creep rate at fixed stress and temperature
# Creep rate = A × σ^n × exp(-Q/RT)
# This is a HYBRID: barrier (exp(-Q/RT)) + response (σ^n)
# Can γ predict the pre-exponential part?

print("""
Creep is a HYBRID regime property:
  ε̇ = A × σⁿ × exp(-Q/RT)

  The activation energy Q dominates (Regime 3)
  But the stress exponent n and pre-factor A relate to
  dislocation dynamics which may involve coherence.

  Prediction: γ alone will show weak correlation (r ~ 0.3-0.5)
  because Q dominates but A has some coherence content.

  Testing this properly requires holding Q constant (within
  the same material class) and looking at A variation.
  With our cross-material dataset, Q varies too much to isolate γ.

  SKIP: This test requires controlled experiments, not cross-material data.
""")

# ==============================================================================
# Summary of all tests
# ==============================================================================
print("\n" + "=" * 70)
print("SUMMARY: FOUR-REGIME PREDICTION TEST RESULTS")
print("=" * 70)

results = [
    ('Speed of sound v_s', 'Regime 1', f'γ^{slope_v:.2f}', f'{r_v_ll:.3f}', 'CONFIRMED (partially tautological)'),
    ('Hardness H_v', 'Regime 1', f'γ^{slope_H:.2f}', f'{r_H_ll:.3f}', 'CONFIRMED (confounded by material class)'),
    ('Ductility (elongation)', 'Regime 2', f'r = {r_elong_gamma:.3f}', f'{r_elong_gamma:.3f}',
     'WEAK' if abs(r_elong_gamma) < 0.3 else 'MODERATE'),
    ('Melting point T_m', 'Regime 1', f'r = {r_Tm_gamma:.3f}', f'{r_Tm_gamma:.3f}',
     'Lindemann rediscovery'),
    ('Dielectric loss tan(δ)', 'Regime 2', f'γ^{slope_tan:.2f}', f'{r_tan_ll:.3f}',
     'CONFIRMED' if r_tan_gamma > 0.3 else 'UNCERTAIN'),
]

print(f"\n{'Property':<25} {'Predicted':<12} {'Scaling':<15} {'|r|':<8} {'Status'}")
print("-" * 75)
for prop, regime, scaling, r_val, status in results:
    print(f"{prop:<25} {regime:<12} {scaling:<15} {r_val:<8} {status}")

# Count outcomes
confirmed = sum(1 for _, _, _, _, s in results if 'CONFIRMED' in s)
total = len(results)
print(f"\nConfirmed: {confirmed}/{total}")
print(f"Partially tautological: 2/{total} (v_s and T_m relate to θ_D by definition)")

# ==============================================================================
# What we actually learned
# ==============================================================================
print("\n" + "=" * 70)
print("WHAT WE ACTUALLY LEARNED")
print("=" * 70)

print(f"""
1. HONEST ASSESSMENT: The regime classification is USEFUL but not all
   predictions are genuinely novel:

   a) Speed of sound: v_s ∝ θ_D by definition → TAUTOLOGICAL
   b) Melting point: Lindemann criterion rediscovered → NOT NEW
   c) Hardness: Confirmed but confounded by material class
   d) Dielectric loss: {'Genuinely predictive' if r_tan_gamma > 0.3 else 'Inconclusive'}
   e) Ductility: Crystal structure dominates over γ

2. THE FRAMEWORK'S DOMAIN IS NARROWER THAN CLAIMED
   Many "predictions" are really restatements of known relationships
   in γ language. Genuine predictive power requires:
   - Property NOT derivable from θ_D alone
   - Correlation NOT confounded by material class
   - Physical mechanism NOT tautological

3. THE FRAMEWORK'S GENUINE STRENGTHS
   From Era 1 sessions, the framework genuinely predicts:
   - Superconductivity: Tc ∝ exp(-γ/λ) — combines γ with coupling
   - Piezoelectricity: d ∝ γ × ε — two-regime with material-specific ε
   - Optical Huang-Rhys: S vs γ_optical — non-trivial prediction
   - Electron transfer: k_ET combined model — non-trivial prediction

   These work because they combine γ with OTHER variables (λ, ε, donor-
   acceptor geometry), creating predictions that go beyond θ_D alone.

4. THE KEY INSIGHT
   γ = 2T/θ_D is useful when it COMBINES with other parameters.
   Alone, it's just a temperature-normalized Debye temperature.
   Its power lies in the COMBINATION: γ × ε, γ/λ_ep, γ × CF.
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Phase 2 Session #6: Testing Four-Regime Predictions\nNew Properties Not Previously Analyzed',
             fontsize=14, fontweight='bold')

# Plot 1: Speed of sound
ax = axes[0, 0]
ax.scatter(gamma_s, v_s, c='steelblue', s=60, alpha=0.7)
for i, name in enumerate(names_s):
    ax.annotate(name, (gamma_s[i], v_s[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('Speed of Sound (m/s)')
ax.set_title(f'v_s ∝ γ^{slope_v:.2f} (r={r_v_ll:.3f})\n*Partially tautological')
ax.grid(True, alpha=0.3)

# Plot 2: Hardness
ax = axes[0, 1]
ax.scatter(gamma_h, H_v, c='crimson', s=60, alpha=0.7)
for i, name in enumerate(names_h):
    ax.annotate(name, (gamma_h[i], H_v[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('Vickers Hardness (GPa)')
ax.set_yscale('log')
ax.set_title(f'H ∝ γ^{slope_H:.2f} (r={r_H_ll:.3f})')
ax.grid(True, alpha=0.3)

# Plot 3: Ductility
ax = axes[0, 2]
ax.scatter(gamma_d, elong, c='green', s=60, alpha=0.7)
for i, name in enumerate(names_d):
    ax.annotate(name, (gamma_d[i], elong[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('Elongation at Break (%)')
ax.set_title(f'Ductility vs γ: r={r_elong_gamma:.3f}\n(crystal structure dominates)')
ax.grid(True, alpha=0.3)

# Plot 4: Melting point / γ_melt
ax = axes[1, 0]
ax.scatter(theta_D_m, gamma_at_melt, c='orange', s=60, alpha=0.7)
for i, name in enumerate(names_m):
    ax.annotate(name, (theta_D_m[i], gamma_at_melt[i]), fontsize=6, alpha=0.6)
ax.axhline(y=gamma_at_melt.mean(), color='red', linestyle='--', alpha=0.5,
           label=f'Mean γ_melt = {gamma_at_melt.mean():.1f}')
ax.set_xlabel('θ_D (K)')
ax.set_ylabel('γ at melting point')
ax.set_title(f'Lindemann: γ_melt = {gamma_at_melt.mean():.1f} ± {gamma_at_melt.std():.1f}')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 5: Dielectric loss
ax = axes[1, 1]
ax.scatter(gamma_dl, tan_delta, c='purple', s=60, alpha=0.7)
for i, name in enumerate(names_dl):
    ax.annotate(name, (gamma_dl[i], tan_delta[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('tan(δ) at 1 MHz')
ax.set_yscale('log')
ax.set_title(f'tan(δ) vs γ: r={r_tan_gamma:.3f}')
ax.grid(True, alpha=0.3)

# Plot 6: Scorecard
ax = axes[1, 2]
ax.text(0.5, 0.92, 'Prediction Scorecard', fontsize=14, ha='center', fontweight='bold', transform=ax.transAxes)

tests = [
    ('Speed of sound', 'Regime 1', 'Tautological', 'gray'),
    ('Hardness', 'Regime 1', 'Confirmed*', 'blue'),
    ('Ductility', 'Regime 2', 'Weak', 'red'),
    ('Melting point', 'Regime 1', 'Lindemann', 'gray'),
    ('Dielectric loss', 'Regime 2', f'r={r_tan_gamma:.2f}', 'blue' if r_tan_gamma > 0.3 else 'red'),
]

for i, (name, regime, result, color) in enumerate(tests):
    y = 0.78 - i * 0.14
    ax.text(0.05, y, f'{name}:', fontsize=10, transform=ax.transAxes)
    ax.text(0.50, y, regime, fontsize=10, transform=ax.transAxes)
    ax.text(0.72, y, result, fontsize=10, transform=ax.transAxes, color=color, fontweight='bold')

ax.text(0.5, 0.08, 'Gray = tautological, Blue = confirmed, Red = failed',
        fontsize=9, ha='center', transform=ax.transAxes, alpha=0.6)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase2_regime_predictions.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase2_regime_predictions.png")
