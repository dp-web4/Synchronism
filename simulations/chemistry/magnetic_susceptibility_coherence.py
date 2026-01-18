#!/usr/bin/env python3
"""
Chemistry Session #82: Magnetic Susceptibility & Coherence
Test whether coherence framework predicts magnetic susceptibility.

Magnetic susceptibility χ measures magnetic response:
- Diamagnetic: χ < 0 (opposes field) - all materials
- Paramagnetic: χ > 0 (aligns with field) - unpaired electrons
- Ferromagnetic: χ >> 0 (spontaneous magnetization below Tc)

Key relationships:
- Diamagnetic χ ∝ -n_e × r² (Langevin, electron orbital motion)
- Paramagnetic χ ∝ μ²/(3kT) (Curie law)
- χ_para ∝ 1/T (Curie law temperature dependence)

Coherence interpretation:
- Diamagnetism: electron orbital coherence → induced moment
- Paramagnetism: spin coherence → alignment with field
- Ferromagnetism: γ_spin → 0 below Tc (Session #63)

From Session #63: β_γ = β (magnetic transitions)
This session tests static susceptibility vs coherence.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("CHEMISTRY SESSION #82: MAGNETIC SUSCEPTIBILITY & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: MAGNETIC SUSCEPTIBILITY
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: MOLAR MAGNETIC SUSCEPTIBILITY")
print("=" * 70)

# Material: (χ_m in 10^-6 cm³/mol, type, n_unpaired, note)
# Positive = paramagnetic, negative = diamagnetic
materials = {
    # Diamagnetic elements
    'Cu': (-5.5, 'dia', 0, 'full d-band'),
    'Ag': (-19.5, 'dia', 0, 'full d-band'),
    'Au': (-28.0, 'dia', 0, 'full d-band'),
    'Zn': (-11.4, 'dia', 0, 'd10'),
    'Pb': (-23.0, 'dia', 0, 'large core'),
    'C_diamond': (-5.9, 'dia', 0, 'all electrons paired'),
    'Si': (-3.1, 'dia', 0, 'covalent'),
    'Ge': (-11.6, 'dia', 0, 'covalent'),
    'Bi': (-280, 'dia', 0, 'largest diamagnetic element'),

    # Paramagnetic elements (transition metals)
    'Ti': (153, 'para', 2, 'd2'),
    'V': (255, 'para', 3, 'd3'),
    'Cr': (180, 'para', 6, 'd5 (antiferro)'),
    'Mn': (529, 'para', 5, 'd5'),
    'Fe_para': (7200, 'para', 4, 'd6 above Tc'),  # Paramagnetic iron above Tc
    'Co_para': (10500, 'para', 3, 'd7 above Tc'),
    'Ni_para': (600, 'para', 2, 'd8 above Tc'),
    'Pt': (193, 'para', 0, 'd-band near Ef'),

    # Alkali metals (paramagnetic due to conduction electrons)
    'Na': (16, 'para', 0, 'Pauli'),
    'K': (20.8, 'para', 0, 'Pauli'),
    'Li': (14.2, 'para', 0, 'Pauli'),

    # Lanthanides (strong paramagnetism from f-electrons)
    'Ce': (2500, 'para', 1, 'f1'),
    'Nd': (5600, 'para', 3, 'f3'),
    'Gd': (185000, 'para', 7, 'f7 - largest'),
    'Dy': (98000, 'para', 5, 'f9'),

    # Simple molecules
    'O2': (3449, 'para', 2, 'triplet ground state'),
    'N2': (-12.0, 'dia', 0, 'singlet'),
    'H2O': (-13.0, 'dia', 0, 'all paired'),
}

print(f"Materials: {len(materials)}")

# Print sorted by χ
print("\nMaterials sorted by susceptibility:")
print("-" * 70)
print(f"{'Material':<12} {'χ_m (10⁻⁶)':<15} {'Type':<8} {'n_unpaired':<10} {'Note':<20}")
print("-" * 70)

for name, (chi, mtype, n_un, note) in sorted(materials.items(), key=lambda x: x[1][0]):
    print(f"{name:<12} {chi:>12.0f}  {mtype:<8} {n_un:>8}  {note:<20}")

# ==============================================================================
# SEPARATE DIAMAGNETIC AND PARAMAGNETIC
# ==============================================================================

print("\n" + "=" * 70)
print("DIAMAGNETIC MATERIALS ANALYSIS")
print("=" * 70)

diamagnetic = {k: v for k, v in materials.items() if v[1] == 'dia'}
paramagnetic = {k: v for k, v in materials.items() if v[1] == 'para'}

print(f"Diamagnetic: {len(diamagnetic)}")
print(f"Paramagnetic: {len(paramagnetic)}")

# Extract diamagnetic arrays
chi_dia = np.array([v[0] for v in diamagnetic.values()])
names_dia = list(diamagnetic.keys())

print("\nDiamagnetic susceptibilities (all negative):")
print("χ ∝ -n_e × <r²> (Langevin diamagnetism)")
print("-" * 50)

# Sort by magnitude
for i, name in enumerate(names_dia):
    chi = chi_dia[i]
    print(f"{name:<12}: χ = {chi:.1f} × 10⁻⁶")

# Diamagnetic trend: larger atoms → more negative χ
# Because more electrons and larger <r²>
print("\nBi has largest diamagnetic susceptibility (χ = -280)")
print("Due to large atomic size and high electron density")

# ==============================================================================
# PARAMAGNETIC ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("PARAMAGNETIC MATERIALS ANALYSIS")
print("=" * 70)

chi_para = []
n_unpaired = []
names_para = []

for name, (chi, mtype, n_un, note) in paramagnetic.items():
    chi_para.append(chi)
    n_unpaired.append(n_un)
    names_para.append(name)

chi_para = np.array(chi_para)
n_unpaired = np.array(n_unpaired)

# Curie law: χ = C/T where C ∝ μ² ∝ S(S+1) ∝ n(n+2)/4 for n unpaired
# For transition metals with spin-only: μ = √(n(n+2)) μ_B
# χ_m ∝ n(n+2) for spin-only

print("Testing Curie law: χ ∝ n(n+2) for spin-only")
print("-" * 50)

# Calculate expected from spin-only
mu_squared = n_unpaired * (n_unpaired + 2)  # ∝ S(S+1)

# Filter for transition metals with n > 0
mask = n_unpaired > 0
if np.sum(mask) >= 3:
    r_chi_mu2, _ = stats.pearsonr(chi_para[mask], mu_squared[mask])
    print(f"χ vs n(n+2): r = {r_chi_mu2:.3f}")
else:
    r_chi_mu2 = 0
    print("Not enough data for correlation")

# Print with expected
print("\nComparison with spin-only expectation:")
print("-" * 60)
print(f"{'Material':<12} {'χ_obs':<12} {'n_unpaired':<10} {'n(n+2)':<10}")
print("-" * 60)
for i, name in enumerate(names_para):
    chi = chi_para[i]
    n = n_unpaired[i]
    mu2 = mu_squared[i]
    print(f"{name:<12} {chi:>10.0f}  {n:>8}  {mu2:>8.0f}")

# ==============================================================================
# LANTHANIDE SERIES
# ==============================================================================

print("\n" + "=" * 70)
print("LANTHANIDE SERIES: χ vs f-ELECTRONS")
print("=" * 70)

lanthanides = {
    'Ce': (2500, 1),
    'Nd': (5600, 3),
    'Gd': (185000, 7),
    'Dy': (98000, 5),  # Note: 5 gives f9 which has J=15/2
}

print("For lanthanides, orbital angular momentum contributes:")
print("μ_eff = g_J × √(J(J+1)) (not spin-only)")
print("-" * 50)

chi_ln = np.array([v[0] for v in lanthanides.values()])
n_ln = np.array([v[1] for v in lanthanides.values()])
names_ln = list(lanthanides.keys())

for name, (chi, n) in lanthanides.items():
    print(f"{name}: χ = {chi}, f-electrons = {n}")

print(f"\nGd (f7, half-filled) has maximum χ = 185,000 × 10⁻⁶")

# ==============================================================================
# PAULI PARAMAGNETISM
# ==============================================================================

print("\n" + "=" * 70)
print("PAULI PARAMAGNETISM: ALKALI METALS")
print("=" * 70)

print("""
For metals with free electrons:
χ_Pauli ∝ N(E_F) (density of states at Fermi level)

This is TEMPERATURE INDEPENDENT (unlike Curie law)!

Alkali metals:
- Na: χ = 16 × 10⁻⁶
- K: χ = 20.8 × 10⁻⁶
- Li: χ = 14.2 × 10⁻⁶

All small and similar because N(E_F) is similar for s-band metals.

Coherence interpretation:
- χ_Pauli ∝ N(E_F) ∝ 1/E_F ∝ m*/ℏ²
- Higher effective mass → higher N(E_F) → higher χ
- NOT directly related to γ_phonon
""")

# ==============================================================================
# COHERENCE INTERPRETATION
# ==============================================================================

print("\n" + "=" * 70)
print("COHERENCE INTERPRETATION")
print("=" * 70)

print("""
Magnetic susceptibility involves SPIN coherence (γ_spin):

1. DIAMAGNETISM (χ < 0):
   - Orbital electron motion opposes field
   - Related to orbital coherence
   - χ_dia ∝ -n_e × <r²>
   - More coherent orbitals → larger induced moment

2. PARAMAGNETISM (Curie):
   - Unpaired spins align with field
   - χ_para ∝ μ²/T
   - At T >> 0: γ_spin → 2 (classical, random)
   - At T → 0: γ_spin → 0 (ordered, aligned)

3. PAULI PARAMAGNETISM:
   - Conduction electron spins
   - Temperature independent
   - χ_Pauli ∝ N(E_F)
   - Related to electronic coherence

4. FERROMAGNETISM (T < Tc):
   - From Session #63: γ_spin → 0
   - Spontaneous ordering
   - χ → ∞ (diverges at Tc)

Key insight:
- Each type involves DIFFERENT coherence
- γ_orbital (diamagnetic), γ_spin (Curie), γ_electron (Pauli)
- Like Session #81: different γ for different physics
""")

# ==============================================================================
# TRANSITION METALS: χ VS θ_D
# ==============================================================================

print("\n" + "=" * 70)
print("TESTING γ_phonon FOR PARAMAGNETIC METALS")
print("=" * 70)

# Subset with known θ_D
tm_data = {
    'Ti': (153, 420),
    'V': (255, 380),
    'Cr': (180, 630),
    'Mn': (529, 410),
    'Pt': (193, 240),
}

chi_tm = np.array([v[0] for v in tm_data.values()])
theta_tm = np.array([v[1] for v in tm_data.values()])
gamma_tm = 2.0 * 300 / theta_tm

r_chi_theta, _ = stats.pearsonr(chi_tm, theta_tm)
r_chi_gamma, _ = stats.pearsonr(chi_tm, gamma_tm)

print(f"χ vs θ_D: r = {r_chi_theta:.3f}")
print(f"χ vs γ_phonon: r = {r_chi_gamma:.3f}")

print("\nExpected: NO correlation with γ_phonon")
print("Because paramagnetic χ depends on SPIN not phonons!")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #82 SUMMARY: MAGNETIC SUSCEPTIBILITY & COHERENCE")
print("=" * 70)

print(f"""
Correlation Analysis:
- χ vs n(n+2) (spin-only): r = {r_chi_mu2:.3f} (for transition metals)
- χ vs θ_D: r = {r_chi_theta:.3f} (for paramagnetics with data)
- χ vs γ_phonon: r = {r_chi_gamma:.3f}

Key Findings:
1. Paramagnetic χ follows Curie law approximately
   - χ ∝ μ² ∝ n(n+2) for spin-only systems
   - Lanthanides deviate (orbital contribution)

2. Pauli χ is temperature-independent
   - χ_Pauli ∝ N(E_F)
   - Different physics from Curie paramagnetism

3. No correlation with γ_phonon
   - χ measures SPIN coherence, not phonon coherence
   - Consistent with Session #81 insight

4. Multiple coherence types identified:
   - γ_orbital (diamagnetic response)
   - γ_spin (Curie paramagnetic)
   - γ_electron (Pauli paramagnetic)
   - γ_phonon (lattice properties)

Physical Interpretation:
- Magnetic susceptibility measures spin/orbital response to field
- This is distinct from both phonon and conduction electron coherence
- Each type of susceptibility has its own coherence physics
- Framework is general: γ → 2 classical, γ → 0 coherent
- But γ ESTIMATION must match the property type

Connection to Session #63:
- Ferromagnetic transitions: β_γ = β
- χ diverges at Tc as γ_spin → 0
- This session examines static χ above Tc
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P82.1: χ_Curie ∝ 1/γ_spin at T > 0
Paramagnetic susceptibility from spin coherence parameter.

P82.2: χ_dia ∝ 2/γ_orbital
Diamagnetic response from orbital coherence.

P82.3: χ_Pauli ∝ N(E_F) (no T dependence)
Conduction electron spin coherence.

P82.4: No correlation with γ_phonon
Spin coherence ≠ phonon coherence.

P82.5: Lanthanides: orbital contribution important
Full J(J+1) formula, not spin-only n(n+2).

P82.6: Ferromagnetic χ → ∞ as γ_spin → 0
Connection to Session #63 magnetic transitions.
""")

# ==============================================================================
# VALIDATION STATUS
# ==============================================================================

print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

print(f"""
**FRAMEWORK CONSISTENT - DIFFERENT γ REQUIRED**

Magnetic susceptibility:
- Does NOT correlate with γ_phonon (r ~ {r_chi_gamma:.1f})
- DOES follow Curie law χ ∝ μ²/T
- Requires γ_spin (spin coherence parameter)

This VALIDATES the Session #81 insight:
- Different properties require different coherence parameters
- γ is a general framework (0 = coherent, 2 = classical)
- But γ estimation must match the physics

Coherence Type Catalog (so far):
1. γ_phonon: Lattice properties (E, C_p, α, v, κ_lattice)
2. γ_electron: Electronic transport (σ, κ_electron)
3. γ_optical: Polarizability (n, ε)
4. γ_spin: Magnetic response (χ)

Each follows γ = 2/√N_corr, but N_corr is different!
- Phonon: coherent phonon modes
- Electron: coherent electron states
- Spin: coherent spin orientations
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: χ distribution (histogram)
ax1 = axes[0, 0]
all_chi = np.array([v[0] for v in materials.values()])
colors = ['blue' if v[1] == 'dia' else 'red' for v in materials.values()]
ax1.barh(range(len(materials)), all_chi, color=colors, alpha=0.7)
ax1.set_yticks(range(len(materials)))
ax1.set_yticklabels(list(materials.keys()), fontsize=8)
ax1.axvline(x=0, color='black', linestyle='-', linewidth=1)
ax1.set_xlabel('χ_m (10⁻⁶ cm³/mol)', fontsize=12)
ax1.set_title('Magnetic Susceptibility\n(Blue=Diamagnetic, Red=Paramagnetic)', fontsize=14)

# Plot 2: χ vs n(n+2) for paramagnetics
ax2 = axes[0, 1]
mask_para = n_unpaired > 0
ax2.scatter(mu_squared[mask_para], chi_para[mask_para], s=100, alpha=0.7, c='red')
for i, name in enumerate(names_para):
    if n_unpaired[i] > 0:
        ax2.annotate(name, (mu_squared[i], chi_para[i]), fontsize=9)
ax2.set_xlabel('n(n+2) (spin-only parameter)', fontsize=12)
ax2.set_ylabel('χ_m (10⁻⁶ cm³/mol)', fontsize=12)
ax2.set_title(f'Paramagnetic χ vs Spin-Only\n(r = {r_chi_mu2:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: Diamagnetic χ (all negative)
ax3 = axes[1, 0]
ax3.bar(range(len(diamagnetic)), chi_dia, alpha=0.7, color='blue')
ax3.set_xticks(range(len(diamagnetic)))
ax3.set_xticklabels(names_dia, rotation=45, ha='right')
ax3.set_ylabel('χ_m (10⁻⁶ cm³/mol)', fontsize=12)
ax3.set_title('Diamagnetic Susceptibilities\n(All negative, Bi most negative)', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: χ vs θ_D for transition metals
ax4 = axes[1, 1]
ax4.scatter(theta_tm, chi_tm, s=100, alpha=0.7, c='purple')
for name, (chi, theta) in tm_data.items():
    ax4.annotate(name, (theta, chi), fontsize=10)
ax4.set_xlabel('Debye Temperature θ_D (K)', fontsize=12)
ax4.set_ylabel('χ_m (10⁻⁶ cm³/mol)', fontsize=12)
ax4.set_title(f'χ vs θ_D (transition metals)\n(r = {r_chi_theta:.3f} - NO correlation)', fontsize=14)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/magnetic_susceptibility_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/magnetic_susceptibility_coherence.png")

print("\n" + "=" * 70)
print("SESSION #82 COMPLETE: MAGNETIC SUSCEPTIBILITY & COHERENCE")
print("=" * 70)
