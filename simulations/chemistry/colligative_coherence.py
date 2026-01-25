"""
Chemistry Session #198: Colligative Properties at γ ~ 1
Analyzing colligative phenomena through the coherence framework.

Colligative properties depend only on NUMBER of particles, not identity.
This is a perfect coherence scenario - counting correlations.

Key γ ~ 1 parameters:
- van't Hoff factor i (observed/expected particles)
- Raoult's law deviations
- Osmotic coefficient φ
- Activity coefficients
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("CHEMISTRY SESSION #198: COLLIGATIVE PROPERTIES AT γ ~ 1")
print("=" * 60)

# ============================================================
# 1. VAN'T HOFF FACTOR
# ============================================================
print("\n" + "=" * 60)
print("1. VAN'T HOFF FACTOR i")
print("=" * 60)

print("""
van't Hoff factor i:
  i = (actual property) / (expected for non-electrolyte)

For ideal behavior:
  Non-electrolytes: i = 1
  Strong electrolytes: i = ν (number of ions)

γ_vH = i_observed / i_ideal

At γ = 1: ideal behavior (no ion pairing, complete dissociation)
γ < 1: incomplete dissociation or ion pairing
γ > 1: rare (super-dissociation)

Ion pairing reduces effective particle count → γ < 1
""")

# van't Hoff factor data (solute, i_ideal, i_observed at 0.1 M)
vH_data = {
    # Non-electrolytes
    'Glucose': (1, 1.00),
    'Sucrose': (1, 1.00),
    'Urea': (1, 1.00),
    'Ethanol': (1, 1.00),
    # Strong electrolytes (1:1)
    'NaCl': (2, 1.87),
    'KCl': (2, 1.85),
    'NaBr': (2, 1.86),
    'KBr': (2, 1.84),
    'HCl': (2, 1.90),
    'NaNO3': (2, 1.83),
    # Strong electrolytes (1:2 or 2:1)
    'CaCl2': (3, 2.70),
    'MgCl2': (3, 2.68),
    'Na2SO4': (3, 2.60),
    'K2SO4': (3, 2.55),
    # Strong electrolytes (2:2)
    'MgSO4': (2, 1.21),
    'CaSO4': (2, 1.15),
    # Weak electrolytes
    'Acetic acid': (2, 1.04),
    'NH3': (2, 1.01),
}

print("\nvan't Hoff Factor Analysis:")
print("-" * 65)
print(f"{'Solute':<15} {'i_ideal':<10} {'i_obs':<10} {'γ = i_obs/i_ideal'}")
print("-" * 65)

gamma_values = []
for solute, (i_ideal, i_obs) in vH_data.items():
    gamma = i_obs / i_ideal
    gamma_values.append(gamma)
    print(f"{solute:<15} {i_ideal:<10} {i_obs:<10.2f} {gamma:.3f}")

print("-" * 65)
mean_gamma = np.mean(gamma_values)
std_gamma = np.std(gamma_values)
print(f"Mean γ = {mean_gamma:.3f} ± {std_gamma:.3f}")

# Separate by type
non_elec = [vH_data[s][1]/vH_data[s][0] for s in ['Glucose', 'Sucrose', 'Urea', 'Ethanol']]
strong_11 = [vH_data[s][1]/vH_data[s][0] for s in ['NaCl', 'KCl', 'NaBr', 'KBr', 'HCl', 'NaNO3']]
strong_12 = [vH_data[s][1]/vH_data[s][0] for s in ['CaCl2', 'MgCl2', 'Na2SO4', 'K2SO4']]
strong_22 = [vH_data[s][1]/vH_data[s][0] for s in ['MgSO4', 'CaSO4']]
weak = [vH_data[s][1]/vH_data[s][0] for s in ['Acetic acid', 'NH3']]

print(f"\nBy electrolyte type:")
print(f"  Non-electrolytes: γ = {np.mean(non_elec):.3f} (ideal!)")
print(f"  1:1 strong: γ = {np.mean(strong_11):.3f} (slight ion pairing)")
print(f"  1:2/2:1 strong: γ = {np.mean(strong_12):.3f} (more pairing)")
print(f"  2:2 strong: γ = {np.mean(strong_22):.3f} (significant pairing)")
print(f"  Weak: γ = {np.mean(weak):.3f} (incomplete dissociation)")

# ============================================================
# 2. RAOULT'S LAW
# ============================================================
print("\n" + "=" * 60)
print("2. RAOULT'S LAW AND ACTIVITY")
print("=" * 60)

print("""
Raoult's law:
  P = x × P°  (ideal)
  P = a × P° = γ × x × P° (real)

where γ is activity coefficient, a = γx is activity.

At γ = 1: ideal solution (Raoult's law)
γ > 1: positive deviation (less favorable mixing)
γ < 1: negative deviation (more favorable mixing)

γ_Raoult IS the activity coefficient!
""")

# Activity coefficient data for binary mixtures at x = 0.5
activity_data = {
    # (solute, solvent): (γ_solute, γ_solvent)
    ('Ethanol', 'Water'): (1.8, 2.5),
    ('Acetone', 'Water'): (4.5, 5.0),
    ('Methanol', 'Water'): (1.3, 1.4),
    ('Benzene', 'Toluene'): (1.05, 1.03),
    ('n-Hexane', 'n-Heptane'): (1.01, 1.01),
    ('Ethanol', 'Benzene'): (3.5, 2.8),
    ('Acetone', 'Chloroform'): (0.65, 0.70),
    ('Chloroform', 'Acetone'): (0.70, 0.65),
}

print("\nActivity Coefficients (γ) at x ≈ 0.5:")
print("-" * 55)
print(f"{'System':<30} {'γ₁':<10} {'γ₂':<10} {'Deviation'}")
print("-" * 55)

for (sol1, sol2), (g1, g2) in activity_data.items():
    avg_g = (g1 + g2) / 2
    if avg_g > 1.1:
        dev = "Positive"
    elif avg_g < 0.9:
        dev = "Negative"
    else:
        dev = "Ideal (γ~1)"
    print(f"{sol1 + '/' + sol2:<30} {g1:<10.2f} {g2:<10.2f} {dev}")

print("-" * 55)
print("Benzene/Toluene and Hexane/Heptane: γ ~ 1 (similar molecules)")
print("Acetone/Chloroform: γ < 1 (H-bonding → favorable)")

# ============================================================
# 3. FREEZING POINT DEPRESSION
# ============================================================
print("\n" + "=" * 60)
print("3. FREEZING POINT DEPRESSION")
print("=" * 60)

print("""
Freezing point depression:
  ΔT_f = K_f × m × i

K_f is the cryoscopic constant:
  K_f = R × T_f² × M / (1000 × ΔH_f)

The ratio ΔT_f_obs / ΔT_f_ideal = γ_fp

At γ = 1: ideal colligative behavior
γ < 1: ion pairing or incomplete dissociation

Water: K_f = 1.86 K·kg/mol
Benzene: K_f = 5.12 K·kg/mol
""")

# K_f values
Kf_data = {
    'Water': 1.86,
    'Benzene': 5.12,
    'Cyclohexane': 20.0,
    'Camphor': 40.0,
    'Naphthalene': 6.94,
    'Acetic acid': 3.90,
}

print("\nCryoscopic Constants K_f:")
print("-" * 40)
for solvent, Kf in Kf_data.items():
    print(f"{solvent:<20} K_f = {Kf:.2f} K·kg/mol")

# ============================================================
# 4. OSMOTIC COEFFICIENT
# ============================================================
print("\n" + "=" * 60)
print("4. OSMOTIC COEFFICIENT φ")
print("=" * 60)

print("""
Osmotic coefficient φ:
  φ = (π_obs) / (π_ideal) = (ln a₁) / (ln x₁)

Related to activity coefficient by:
  φ = 1 - (1/m) × ∫₀^m (1 - γ±) dm

At φ = 1: ideal solution
φ < 1: ion pairing (common for electrolytes)
φ > 1: enhanced activity

This is IDENTICAL to γ_vH interpretation!
""")

# Osmotic coefficients at 0.1 M
phi_data = {
    'NaCl': 0.932,
    'KCl': 0.926,
    'CaCl2': 0.856,
    'MgSO4': 0.575,
    'H2SO4': 0.676,
    'Sucrose': 1.004,
    'Glucose': 1.000,
    'Urea': 0.998,
}

print("\nOsmotic Coefficients φ at 0.1 M:")
print("-" * 45)
print(f"{'Solute':<15} {'φ':<10} {'Status'}")
print("-" * 45)

phi_values = []
for solute, phi in phi_data.items():
    phi_values.append(phi)
    if 0.95 <= phi <= 1.05:
        status = "γ ~ 1 (ideal)"
    elif phi < 0.95:
        status = "Ion pairing"
    else:
        status = "Enhanced"
    print(f"{solute:<15} {phi:<10.3f} {status}")

print("-" * 45)
print(f"Mean φ = {np.mean(phi_values):.3f} ± {np.std(phi_values):.3f}")

# ============================================================
# 5. EBULLIOSCOPIC AND CRYOSCOPIC RATIOS
# ============================================================
print("\n" + "=" * 60)
print("5. K_b / K_f RATIO")
print("=" * 60)

print("""
Ebullioscopic (K_b) and Cryoscopic (K_f) constants:

K_b = R × T_b² × M / (1000 × ΔH_vap)
K_f = R × T_f² × M / (1000 × ΔH_f)

The ratio:
γ_K = K_b / K_f = (T_b² / T_f²) × (ΔH_f / ΔH_vap)

For water:
  K_b = 0.512, K_f = 1.86
  γ_K = 0.512 / 1.86 = 0.275

This reflects the thermodynamic asymmetry:
  ΔH_vap >> ΔH_f typically
""")

# Kb/Kf data
Kb_Kf_data = {
    # (K_b, K_f)
    'Water': (0.512, 1.86),
    'Benzene': (2.53, 5.12),
    'Acetic acid': (3.07, 3.90),
    'Cyclohexane': (2.79, 20.0),
    'Naphthalene': (5.65, 6.94),
}

print("\nK_b / K_f Ratios:")
print("-" * 45)
print(f"{'Solvent':<15} {'K_b':<10} {'K_f':<10} {'K_b/K_f'}")
print("-" * 45)

for solvent, (Kb, Kf) in Kb_Kf_data.items():
    ratio = Kb / Kf
    print(f"{solvent:<15} {Kb:<10.3f} {Kf:<10.2f} {ratio:.3f}")

# ============================================================
# 6. VISUALIZATION
# ============================================================
print("\n" + "=" * 60)
print("6. GENERATING VISUALIZATION")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: van't Hoff factor γ distribution
ax1 = axes[0, 0]
ax1.hist(gamma_values, bins=10, edgecolor='black', alpha=0.7, color='steelblue')
ax1.axvline(x=1, color='red', linestyle='--', linewidth=2, label='γ = 1 (ideal)')
ax1.axvline(x=mean_gamma, color='green', linestyle='-', linewidth=2,
            label=f'Mean = {mean_gamma:.3f}')
ax1.set_xlabel('γ = i_obs / i_ideal', fontsize=12)
ax1.set_ylabel('Count', fontsize=12)
ax1.set_title("van't Hoff Factor: γ ~ 1 Reference", fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: By electrolyte type
ax2 = axes[0, 1]
categories = ['Non-\nelectrolyte', '1:1\nStrong', '1:2/2:1\nStrong', '2:2\nStrong', 'Weak\nelectrolyte']
means = [np.mean(non_elec), np.mean(strong_11), np.mean(strong_12),
         np.mean(strong_22), np.mean(weak)]
stds = [np.std(non_elec) if len(non_elec) > 1 else 0,
        np.std(strong_11), np.std(strong_12),
        np.std(strong_22) if len(strong_22) > 1 else 0,
        np.std(weak) if len(weak) > 1 else 0]

bars = ax2.bar(categories, means, yerr=stds, capsize=5,
               color=['green', 'steelblue', 'coral', 'purple', 'orange'],
               edgecolor='black', alpha=0.7)
ax2.axhline(y=1, color='red', linestyle='--', linewidth=2, label='γ = 1')
ax2.set_ylabel('γ = i_obs / i_ideal', fontsize=12)
ax2.set_title('γ by Electrolyte Type', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_ylim(0, 1.2)

# Plot 3: Osmotic coefficient
ax3 = axes[1, 0]
solutes = list(phi_data.keys())
phi_vals = list(phi_data.values())
colors = ['green' if 0.95 <= p <= 1.05 else ('red' if p < 0.95 else 'blue') for p in phi_vals]
ax3.barh(range(len(solutes)), phi_vals, color=colors, edgecolor='black', alpha=0.7)
ax3.axvline(x=1, color='black', linestyle='--', linewidth=2)
ax3.set_yticks(range(len(solutes)))
ax3.set_yticklabels(solutes)
ax3.set_xlabel('Osmotic Coefficient φ', fontsize=12)
ax3.set_title('Osmotic Coefficient: φ = 1 is Ideal', fontsize=14)
ax3.grid(True, alpha=0.3, axis='x')
ax3.set_xlim(0.5, 1.1)

# Plot 4: Activity coefficient regimes
ax4 = axes[1, 1]
x_range = np.linspace(0, 1, 100)
# Ideal
y_ideal = np.ones_like(x_range)
# Positive deviation (Margules)
y_pos = np.exp(0.5 * (1-x_range)**2)
# Negative deviation
y_neg = np.exp(-0.3 * (1-x_range)**2)

ax4.plot(x_range, y_ideal, 'k-', linewidth=2, label='Ideal (γ = 1)')
ax4.plot(x_range, y_pos, 'r-', linewidth=2, label='Positive deviation')
ax4.plot(x_range, y_neg, 'b-', linewidth=2, label='Negative deviation')
ax4.axhline(y=1, color='gray', linestyle=':', alpha=0.5)
ax4.fill_between(x_range, 0.95, 1.05, alpha=0.2, color='green', label='γ ~ 1 region')
ax4.set_xlabel('Mole fraction x', fontsize=12)
ax4.set_ylabel('Activity coefficient γ', fontsize=12)
ax4.set_title('Activity Coefficient: γ = 1 is Ideal', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 1)
ax4.set_ylim(0.5, 2)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/colligative_coherence.png', dpi=150)
print("Saved: colligative_coherence.png")
plt.close()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SESSION #198 SUMMARY: COLLIGATIVE PROPERTIES AT γ ~ 1")
print("=" * 60)

print(f"""
KEY FINDINGS:

1. VAN'T HOFF FACTOR
   γ_vH = i_obs / i_ideal
   Mean γ = {mean_gamma:.3f} ± {std_gamma:.3f}
   Non-electrolytes: γ = 1.00 (ideal!)
   1:1 electrolytes: γ ~ 0.93 (slight ion pairing)
   2:2 electrolytes: γ ~ 0.60 (significant pairing)

2. ACTIVITY COEFFICIENTS
   Raoult's law: γ = 1 is ideal solution
   Similar molecules (benzene/toluene): γ ~ 1.0
   Dissimilar molecules: γ >> 1 (positive deviation)
   H-bonding pairs: γ < 1 (negative deviation)

3. OSMOTIC COEFFICIENT
   φ = 1 is ideal behavior
   Non-electrolytes: φ = 1.00 exactly
   Strong 1:1: φ ~ 0.93
   2:2 (MgSO4): φ ~ 0.58 (strong pairing)

4. FUNDAMENTAL PATTERN
   γ = 1 IS the ideal solution/ideal electrolyte reference
   Deviations measure:
   - Ion pairing (γ < 1 for electrolytes)
   - Non-ideality (γ ≠ 1 for mixtures)

5. THERMODYNAMIC MEANING
   γ = 1 means:
   - Complete dissociation
   - Ideal mixing (no excess properties)
   - Predicted colligative effects

CENTRAL INSIGHT:
Colligative properties ARE a coherence count:
- γ = 1 means all particles contribute fully
- Ion pairing reduces effective count (γ < 1)
- Activity coefficient γ measures deviation from ideal

The "ideal solution" is THE γ ~ 1 reference state
for all thermodynamic mixture properties.

This is the 61st phenomenon type at γ ~ 1!
""")

print("=" * 60)
print("SESSION #198 COMPLETE")
print("=" * 60)
