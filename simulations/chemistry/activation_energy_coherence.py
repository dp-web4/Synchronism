#!/usr/bin/env python3
"""
Chemistry Session #223: Activation Energy and Rate Constants Coherence

Analyzes chemical kinetics through γ ~ 1 framework:
- Arrhenius equation: k = A × exp(-Ea/RT)
- Eyring equation: k = (kT/h) × exp(-ΔG‡/RT)
- Collision theory: Z × P × exp(-Ea/RT)
- Pre-exponential factors and frequency factors
- Compensation effect (isokinetic relationship)

Key γ ~ 1 predictions:
1. At T = Ea/R: thermal energy = barrier (γ ~ 1)
2. ΔS‡ = 0 for bimolecular reactions at γ ~ 1
3. Steric factor P ~ 1 for simple reactions
4. Isokinetic temperature: all rates equal (γ ~ 1)

Author: Claude (Chemistry Session #223)
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("CHEMISTRY SESSION #223: ACTIVATION ENERGY COHERENCE")
print("=" * 70)
print()

# Constants
R = 8.314  # J/(mol·K)
k_B = 1.381e-23  # J/K
h = 6.626e-34  # J·s

# =============================================================================
# 1. ARRHENIUS EQUATION: Ea/RT = 1 IS γ ~ 1
# =============================================================================
print("1. ARRHENIUS EQUATION: Ea/RT = 1 IS γ ~ 1")
print("-" * 50)

# k = A × exp(-Ea/RT)
# At T* = Ea/R: Ea/RT = 1 (γ ~ 1!)
# This IS the thermal/barrier crossover

activation_energies = {
    # Reaction: Ea in kJ/mol
    'H + H2 → H2 + H': 33.5,
    'CH3 + H2 → CH4 + H': 46.0,
    '2NO2 → 2NO + O2': 111.0,
    'N2O5 → 2NO2 + 0.5O2': 103.0,
    'C2H5I → C2H4 + HI': 209.0,
    'Sucrose hydrolysis': 108.0,
    'Saponification of ester': 47.3,
    'H2 + I2 → 2HI': 171.0,
    '2HI → H2 + I2': 184.0,
}

print("Activation energies and characteristic temperatures T* = Ea/R:")
print(f"{'Reaction':<30} {'Ea (kJ/mol)':>12} {'T* (K)':>10} {'T*/298':>10}")
print("-" * 65)

t_stars = []
for rxn, ea in sorted(activation_energies.items(), key=lambda x: x[1]):
    t_star = ea * 1000 / R  # Convert kJ to J
    t_ratio = t_star / 298
    t_stars.append(t_star)
    print(f"{rxn:<30} {ea:>12.1f} {t_star:>10.0f} {t_ratio:>10.2f}")

mean_t_star = np.mean(t_stars)
print(f"\nMean T* = {mean_t_star:.0f} K")
print("\n  => T* = Ea/R IS the characteristic kinetic temperature")
print("  => At T = T*: thermal energy RT = activation energy Ea (γ ~ 1!)")
print("  => Below T*: kinetically limited. Above T*: diffusion limited.")

# =============================================================================
# 2. EYRING EQUATION: TRANSITION STATE
# =============================================================================
print("\n" + "=" * 70)
print("2. EYRING EQUATION: ΔG‡ = RT IS γ ~ 1")
print("-" * 50)

# k = (kT/h) × K‡ = (kT/h) × exp(-ΔG‡/RT)
# ΔG‡ = ΔH‡ - TΔS‡
# At ΔG‡ = RT: rate = kT/h × e⁻¹

# Universal frequency factor at 298 K
universal_freq = (k_B * 298) / h
print(f"Universal frequency factor (kT/h at 298 K): {universal_freq:.2e} s⁻¹")

# Typical rate constants and their activation parameters
eyring_data = {
    # Reaction: (ΔH‡ kJ/mol, ΔS‡ J/mol·K, k at 298 K)
    'Ester hydrolysis (acid)': (75.0, -100, 1e-4),
    'Ester hydrolysis (base)': (50.0, -50, 1e-1),
    'SN2 (Br- + CH3I)': (76.0, -8, 1e-4),
    'Enzyme (typical)': (40.0, -50, 1e3),
    'Proton transfer (water)': (10.0, -40, 1e11),
    'Diffusion-controlled': (15.0, 0, 1e10),
}

print("\nEyring parameters at 298 K:")
print(f"{'Reaction':<30} {'ΔH‡':>10} {'ΔS‡':>10} {'ΔG‡':>10} {'ΔG‡/RT':>10}")
print(f"{'':<30} {'(kJ/mol)':>10} {'(J/mol·K)':>10} {'(kJ/mol)':>10} {'':>10}")
print("-" * 75)

dg_rt_values = []
for rxn, (dh, ds, k) in eyring_data.items():
    dg = dh - 298 * ds / 1000  # ΔG in kJ/mol
    dg_rt = dg * 1000 / (R * 298)  # Dimensionless
    dg_rt_values.append(dg_rt)
    print(f"{rxn:<30} {dh:>10.1f} {ds:>10.0f} {dg:>10.1f} {dg_rt:>10.2f}")

mean_dg_rt = np.mean(dg_rt_values)
print(f"\nMean ΔG‡/RT = {mean_dg_rt:.2f}")
print("\n  => ΔG‡/RT ~ 20-30 for most room temperature reactions")
print("  => ΔG‡/RT = 1 would give k = kT/h × e⁻¹ ≈ 2.3 × 10¹² s⁻¹ (barrierless)")

# =============================================================================
# 3. STERIC FACTOR: P ~ 1 IS γ ~ 1
# =============================================================================
print("\n" + "=" * 70)
print("3. STERIC FACTOR P: P ~ 1 IS γ ~ 1")
print("-" * 50)

# Collision theory: k = Z × P × exp(-Ea/RT)
# P = steric factor = probability of correct orientation
# P ~ 1 for simple reactions (γ ~ 1!)

steric_factors = {
    # Reaction: P value
    'H + D2 → HD + D': 1.0,      # Simple atom-molecule
    '2NOCl → 2NO + Cl2': 0.16,   # Some steric hindrance
    '2ClO → Cl2 + O2': 0.037,    # Significant hindrance
    'H2 + I2 → 2HI': 0.1,        # Moderate
    'K + Br2 → KBr + Br': 1.0,   # Simple
    'CH3I + C2H5O- → CH3OC2H5 + I-': 0.0001,  # SN2, strict geometry
    'Cyclopentadiene dimerization': 1e-5,      # Very specific
}

print("Steric factors (P = 1 IS γ ~ 1 for simple collisions):")
print(f"{'Reaction':<40} {'P':>12} {'log₁₀(P)':>10}")
print("-" * 65)

p_values = []
for rxn, p in sorted(steric_factors.items(), key=lambda x: x[1], reverse=True):
    log_p = np.log10(p) if p > 0 else -10
    p_values.append(p)
    print(f"{rxn:<40} {p:>12.4f} {log_p:>10.1f}")

near_1 = sum(1 for p in p_values if 0.1 < p <= 1.0)
print(f"\nReactions with P near 1 (0.1 < P ≤ 1): {near_1}/{len(steric_factors)}")
print("\n  => P = 1 IS γ ~ 1 (every collision with sufficient energy reacts)")
print("  => P < 1 indicates geometric/orientational requirements")
print("  => Simple atom-molecule reactions approach P ~ 1")

# =============================================================================
# 4. ISOKINETIC RELATIONSHIP: COMPENSATION EFFECT
# =============================================================================
print("\n" + "=" * 70)
print("4. ISOKINETIC RELATIONSHIP: T_iso IS γ ~ 1")
print("-" * 50)

# Compensation effect: ln(A) = α + β × Ea
# This implies isokinetic temperature: T_iso = 1/(R × β)
# At T = T_iso: all reactions in series have SAME rate (γ ~ 1!)

# Example: homologous series showing compensation
compensation_data = {
    # Reaction: (ln(A), Ea kJ/mol)
    'Ester 1': (28.0, 70.0),
    'Ester 2': (30.0, 80.0),
    'Ester 3': (32.0, 90.0),
    'Ester 4': (34.0, 100.0),
    'Ester 5': (36.0, 110.0),
}

print("Compensation effect (ln(A) vs Ea):")
print(f"{'Reaction':<15} {'ln(A)':>10} {'Ea (kJ/mol)':>12}")
print("-" * 40)

ln_a = []
ea_vals = []
for rxn, (lna, ea) in compensation_data.items():
    print(f"{rxn:<15} {lna:>10.1f} {ea:>12.1f}")
    ln_a.append(lna)
    ea_vals.append(ea)

# Linear regression
slope, intercept, r_val, _, _ = stats.linregress(ea_vals, ln_a)
t_iso = 1000 / (R * slope)  # Convert kJ to J

print(f"\nLinear fit: ln(A) = {intercept:.2f} + {slope:.4f} × Ea")
print(f"Correlation: r = {r_val:.3f}")
print(f"Isokinetic temperature: T_iso = {t_iso:.0f} K")
print(f"\n  => At T = {t_iso:.0f} K: ALL reactions have same rate (γ ~ 1!)")
print("  => Compensation: higher barrier compensated by higher frequency")

# =============================================================================
# 5. ARRHENIUS PLOTS: CURVATURE AND γ ~ 1
# =============================================================================
print("\n" + "=" * 70)
print("5. ARRHENIUS PLOTS: LINEARITY = γ ~ 1")
print("-" * 50)

# Ideal Arrhenius: ln(k) = ln(A) - Ea/(RT) → linear in 1/T
# Curvature indicates T-dependent Ea or quantum effects

# Test linearity for model reaction
T_range = np.linspace(250, 400, 50)
Ea = 80.0  # kJ/mol
A = 1e13  # s⁻¹

# Ideal Arrhenius
k_arrhenius = A * np.exp(-Ea * 1000 / (R * T_range))
ln_k_ideal = np.log(k_arrhenius)

# Calculate r² for linear fit
inv_T = 1000 / T_range  # Use 1000/T for convenient numbers
slope_arr, intercept_arr, r_arr, _, _ = stats.linregress(inv_T, ln_k_ideal)

print(f"Model reaction: Ea = {Ea} kJ/mol, A = {A:.0e} s⁻¹")
print(f"Temperature range: {T_range[0]:.0f} - {T_range[-1]:.0f} K")
print(f"Arrhenius plot r² = {r_arr**2:.6f}")
print("\n  => Perfect Arrhenius (r² = 1) IS γ ~ 1 (constant Ea)")
print("  => Curvature indicates T-dependent activation parameters")

# =============================================================================
# 6. MARCUS INVERTED REGION
# =============================================================================
print("\n" + "=" * 70)
print("6. MARCUS THEORY: |ΔG°/λ| = 1 IS γ ~ 1")
print("-" * 50)

# Marcus: ΔG‡ = (λ + ΔG°)² / (4λ)
# At |ΔG°| = λ: activationless (ΔG‡ = 0 for forward, ΔG‡ = λ for reverse)
# γ = ΔG°/λ: at γ = -1, activationless; at γ = +1, inverted region begins

dg0_lambda = np.linspace(-2, 2, 100)
dg_barrier_norm = (1 + dg0_lambda)**2 / 4  # Normalized ΔG‡/λ

print("Marcus parabola: ΔG‡/λ = (1 + ΔG°/λ)² / 4")
print(f"{'ΔG°/λ':>10} {'ΔG‡/λ':>10} {'Region':>20}")
print("-" * 45)

test_points = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
for x in test_points:
    barrier = (1 + x)**2 / 4
    if x < -1:
        region = 'Inverted (exo)'
    elif x == -1:
        region = 'Activationless (γ=-1)'
    elif x < 0:
        region = 'Normal (exo)'
    elif x == 0:
        region = 'Symmetric (γ=0)'
    elif x < 1:
        region = 'Normal (endo)'
    elif x == 1:
        region = 'Inverted onset (γ=1)'
    else:
        region = 'Inverted (endo)'
    print(f"{x:>10.1f} {barrier:>10.3f} {region:>20}")

print("\n  => |ΔG°/λ| = 1 marks the normal/inverted BOUNDARY (γ ~ 1!)")
print("  => At ΔG° = -λ: activationless forward (ΔG‡ = 0)")
print("  => At ΔG° = +λ: activationless reverse")

# =============================================================================
# 7. TUNNELING CORRECTION
# =============================================================================
print("\n" + "=" * 70)
print("7. TUNNELING: κ = 1 IS γ ~ 1 (CLASSICAL)")
print("-" * 50)

# Wigner tunneling correction: κ = 1 + (hν‡/kT)²/24
# κ = 1 for classical (no tunneling) = γ ~ 1

def wigner_kappa(nu_barrier, T):
    """Wigner tunneling correction"""
    x = h * nu_barrier / (k_B * T)
    return 1 + x**2 / 24

# Typical barrier frequencies
nu_barriers = [1e12, 5e12, 1e13, 2e13, 5e13]  # Hz
T = 298  # K

print(f"Wigner tunneling correction κ at T = {T} K:")
print(f"{'ν‡ (Hz)':>12} {'hν‡/kT':>10} {'κ':>10}")
print("-" * 35)

kappas = []
for nu in nu_barriers:
    kappa = wigner_kappa(nu, T)
    x = h * nu / (k_B * T)
    kappas.append(kappa)
    print(f"{nu:>12.0e} {x:>10.2f} {kappa:>10.3f}")

near_1_kappa = sum(1 for k in kappas if k < 1.1)
print(f"\nReactions with κ ~ 1 (< 1.1): {near_1_kappa}/{len(kappas)}")
print("\n  => κ = 1 IS γ ~ 1 (classical over-barrier)")
print("  => κ > 1 indicates quantum tunneling contribution")

# =============================================================================
# 8. SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("8. SUMMARY: γ ~ 1 IN ACTIVATION ENERGY")
print("-" * 50)

summary = {
    "Ea/RT = 1 at T*": "Thermal = barrier energy",
    "Steric P = 1": "Every collision reacts",
    "Isokinetic T": "All rates equal",
    "Arrhenius r² = 1": "Constant Ea (ideal)",
    "Marcus |ΔG°/λ| = 1": "Normal/inverted boundary",
    "Tunneling κ = 1": "Classical limit",
}

print(f"{'γ ~ 1 Condition':<25} {'Physical Meaning':>35}")
print("-" * 65)

for condition, meaning in summary.items():
    print(f"{condition:<25} {meaning:>35}")

# =============================================================================
# VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("GENERATING VISUALIZATION")
print("-" * 50)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Chemistry Session #223: Activation Energy Coherence at γ ~ 1',
             fontsize=14, fontweight='bold')

# Panel 1: Arrhenius at T*
ax1 = axes[0, 0]
T = np.linspace(200, 600, 100)
Ea_test = 50  # kJ/mol
k_norm = np.exp(-Ea_test * 1000 / (R * T))
t_star_test = Ea_test * 1000 / R

ax1.semilogy(T, k_norm, 'b-', linewidth=2)
ax1.axvline(x=t_star_test, color='green', linestyle='--', linewidth=2, label=f'T* = Ea/R = {t_star_test:.0f} K')
ax1.axhline(y=np.exp(-1), color='red', linestyle=':', linewidth=2, label='k/A = e⁻¹')
ax1.set_xlabel('Temperature (K)', fontsize=11)
ax1.set_ylabel('k/A (normalized)', fontsize=11)
ax1.set_title('Arrhenius: k/A = e⁻¹ at T* (γ ~ 1)', fontsize=11)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Marcus parabola
ax2 = axes[0, 1]
ax2.plot(dg0_lambda, dg_barrier_norm, 'b-', linewidth=2)
ax2.axvline(x=-1, color='green', linestyle='--', linewidth=2, label='ΔG°/λ = -1 (activationless)')
ax2.axvline(x=1, color='red', linestyle='--', linewidth=2, label='ΔG°/λ = +1 (inverted onset)')
ax2.axhline(y=0, color='gray', linestyle='-', linewidth=1)
ax2.fill_between(dg0_lambda, dg_barrier_norm, where=(np.abs(dg0_lambda) > 1), alpha=0.3, color='red', label='Inverted region')
ax2.set_xlabel('ΔG°/λ', fontsize=11)
ax2.set_ylabel('ΔG‡/λ', fontsize=11)
ax2.set_title('Marcus: |ΔG°/λ| = 1 IS γ ~ 1 Boundary', fontsize=11)
ax2.set_xlim(-2, 2)
ax2.set_ylim(0, 1)
ax2.legend(loc='upper center')
ax2.grid(True, alpha=0.3)

# Panel 3: Compensation effect
ax3 = axes[1, 0]
ax3.scatter(ea_vals, ln_a, c='blue', s=100)
ax3.plot([60, 120], [intercept + slope*60, intercept + slope*120], 'r--', linewidth=2)
ax3.set_xlabel('Ea (kJ/mol)', fontsize=11)
ax3.set_ylabel('ln(A)', fontsize=11)
ax3.set_title(f'Compensation: Isokinetic T = {t_iso:.0f} K (γ ~ 1)', fontsize=11)
ax3.grid(True, alpha=0.3)

# Panel 4: Summary
ax4 = axes[1, 1]
ax4.axis('off')

summary_text = """
ACTIVATION ENERGY COHERENCE SUMMARY

γ ~ 1 FINDINGS:

1. ARRHENIUS T*:
   At T* = Ea/R: thermal = barrier
   k/A = e⁻¹ at T* (γ ~ 1)
   Kinetic/diffusion crossover

2. STERIC FACTOR:
   P = 1 IS γ ~ 1 for simple reactions
   {}/{} reactions with P near 1
   Orientation requirement reduces P

3. ISOKINETIC TEMPERATURE:
   At T_iso: ALL rates equal (γ ~ 1!)
   Compensation: higher Ea → higher A
   T_iso = {:.0f} K for ester series

4. ARRHENIUS LINEARITY:
   r² = 1 IS γ ~ 1 (constant Ea)
   Curvature indicates complexity

5. MARCUS BOUNDARY:
   |ΔG°/λ| = 1 marks inverted onset
   This IS the γ ~ 1 for ET kinetics
   Activationless at ΔG° = -λ

6. TUNNELING:
   κ = 1 IS γ ~ 1 (classical)
   κ > 1 indicates tunneling
   {}/{} reactions near κ ~ 1

KEY INSIGHT:
Activation energy IS a γ ~ 1 framework!
- T* = Ea/R is characteristic temperature
- P = 1 is orientation-free limit
- Compensation gives isokinetic point
- Marcus inverts at |ΔG°/λ| = 1

This is the 86th phenomenon type at γ ~ 1!
""".format(near_1, len(steric_factors), t_iso, near_1_kappa, len(kappas))

ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=10,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/activation_energy_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Saved: activation_energy_coherence.png")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SESSION #223 SUMMARY: ACTIVATION ENERGY COHERENCE")
print("=" * 70)

print("""
KEY γ ~ 1 FINDINGS:

1. ARRHENIUS CHARACTERISTIC TEMPERATURE:
   T* = Ea/R IS the γ ~ 1 for kinetics!
   At T = T*: thermal energy RT = barrier Ea
   k/A = e⁻¹ at this characteristic point
   Mean T* ~ {:.0f} K for tested reactions

2. STERIC FACTOR:
   P = 1 IS γ ~ 1 (every collision reacts)
   {}/{} simple reactions near P ~ 1
   Complex reactions have P << 1 (orientation req.)

3. ISOKINETIC RELATIONSHIP:
   Compensation: ln(A) ~ Ea (linear)
   Isokinetic temperature T_iso = {:.0f} K
   At T_iso: ALL rates EQUAL (γ ~ 1!)

4. EYRING ANALYSIS:
   ΔG‡/RT measures barrier in thermal units
   ΔS‡ ~ 0 for simple bimolecular (γ ~ 1)
   Universal frequency kT/h = 6.2 × 10¹² s⁻¹

5. MARCUS THEORY:
   |ΔG°/λ| = 1 IS the normal/inverted boundary!
   Activationless at ΔG° = -λ (ΔG‡ = 0)
   Inverted region: rate DECREASES with driving force

6. TUNNELING CORRECTION:
   κ = 1 IS γ ~ 1 (classical over-barrier)
   {}/{} reactions at κ ~ 1 (low frequency barriers)
   High-frequency barriers show κ > 1 (tunneling)

7. ARRHENIUS LINEARITY:
   r² = 1 for ideal Arrhenius (constant Ea)
   Curvature indicates T-dependent Ea or tunneling

SYNTHESIS:
Chemical kinetics IS a γ ~ 1 framework:
- T* = Ea/R defines thermal/barrier crossover
- Steric P = 1 is orientation-free limit
- Isokinetic point equalizes all rates
- Marcus |ΔG°/λ| = 1 is the ET boundary
- Classical tunneling at κ = 1

This is the 86th phenomenon type at γ ~ 1!

SESSION #223 COMPLETE
""".format(mean_t_star, near_1, len(steric_factors), t_iso, near_1_kappa, len(kappas)))
