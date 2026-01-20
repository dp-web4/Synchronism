#!/usr/bin/env python3
"""
Chemistry Session #142: Quantum Criticality and Coherence

Quantum Critical Points (QCPs): Phase transitions at T = 0.

Key physics:
- Tuning parameter g (pressure, field, doping)
- Critical value g_c where order vanishes at T = 0
- Quantum fluctuations replace thermal fluctuations
- "Strange metal" behavior above QCP

Coherence interpretation:
- QCP is where γ = 1 extends down to T = 0
- Quantum critical fan: γ ~ T/T_0 (unusual T dependence)
- NFL behavior from critical fluctuations

Connects to:
- Session #139: Kondo (T_K can be tuned to 0)
- Session #140: Mott (U/W can be tuned through 1)
- Session #141: SC dome (QCP may underlie dome)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Physical constants
kB = 1.381e-23    # J/K
kB_meV = 0.0862   # meV/K

# =============================================================================
# DATASET: Quantum Critical Systems
# =============================================================================

# Collect QCP data from various systems
# Sources: Si, Gegenwart, Sachdev reviews

systems = [
    # Heavy fermion QCPs (pressure-tuned)
    {"name": "CePd2Si2", "type": "heavy_fermion",
     "g_c": 2.8, "g_unit": "GPa", "T_FL": 0.5, "T_star": 10,
     "z": 2, "nu": 0.5, "alpha_rho": 1.0},

    {"name": "CeIn3", "type": "heavy_fermion",
     "g_c": 2.6, "g_unit": "GPa", "T_FL": 0.3, "T_star": 15,
     "z": 3, "nu": 0.5, "alpha_rho": 1.5},

    {"name": "CeRhIn5", "type": "heavy_fermion",
     "g_c": 2.4, "g_unit": "GPa", "T_FL": 0.8, "T_star": 20,
     "z": 2, "nu": 0.5, "alpha_rho": 1.0},

    {"name": "YbRh2Si2", "type": "heavy_fermion",
     "g_c": 0.066, "g_unit": "T", "T_FL": 0.02, "T_star": 25,
     "z": 2, "nu": 0.7, "alpha_rho": 1.0},

    # Itinerant ferromagnet QCPs
    {"name": "ZrZn2", "type": "itinerant_FM",
     "g_c": 21, "g_unit": "GPa", "T_FL": 2.0, "T_star": 30,
     "z": 3, "nu": 0.5, "alpha_rho": 5/3},

    {"name": "MnSi", "type": "itinerant_FM",
     "g_c": 1.5, "g_unit": "GPa", "T_FL": 1.0, "T_star": 30,
     "z": 3, "nu": 0.5, "alpha_rho": 5/3},

    # Cuprate QCPs (doping-tuned)
    {"name": "YBCO (p*)", "type": "cuprate",
     "g_c": 0.19, "g_unit": "doping", "T_FL": 50, "T_star": 300,
     "z": 2, "nu": 0.67, "alpha_rho": 1.0},

    {"name": "LSCO (p*)", "type": "cuprate",
     "g_c": 0.19, "g_unit": "doping", "T_FL": 30, "T_star": 150,
     "z": 2, "nu": 0.67, "alpha_rho": 1.0},

    # Antiferromagnetic QCPs
    {"name": "TlCuCl3", "type": "AFM",
     "g_c": 0.5, "g_unit": "GPa", "T_FL": 5.0, "T_star": 40,
     "z": 1, "nu": 0.5, "alpha_rho": 2.0},

    {"name": "Sr3Ru2O7", "type": "metamagnetic",
     "g_c": 8.0, "g_unit": "T", "T_FL": 1.0, "T_star": 10,
     "z": 2, "nu": 0.5, "alpha_rho": 1.5},
]

# Convert to numpy arrays
names = [s["name"] for s in systems]
types = np.array([s["type"] for s in systems])
g_c = np.array([s["g_c"] for s in systems])
T_FL = np.array([s["T_FL"] for s in systems])  # Fermi liquid scale
T_star = np.array([s["T_star"] for s in systems])  # Crossover scale
z = np.array([s["z"] for s in systems])  # Dynamic exponent
nu = np.array([s["nu"] for s in systems])  # Correlation length exponent
alpha_rho = np.array([s["alpha_rho"] for s in systems])  # Resistivity exponent

print("=" * 70)
print("CHEMISTRY SESSION #142: Quantum Criticality and Coherence")
print("=" * 70)
print(f"\nSystems analyzed: {len(systems)}")
print(f"Types: {np.unique(types)}")
print(f"z (dynamic exponent) range: {z.min():.0f} - {z.max():.0f}")

# =============================================================================
# COHERENCE AT QCP
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: Coherence at Quantum Critical Point")
print("=" * 70)

# At QCP: γ(T, g=g_c) ~ (T/T_0)^(1/zν)
# This is different from classical: γ(T) ~ T/T_c

# Define quantum critical coherence
# γ_QC = (T/T_star)^(1/zν)
def gamma_QC(T, T_star, z, nu):
    """Quantum critical coherence parameter"""
    return (T / T_star) ** (1 / (z * nu))

# At T = T_star: γ_QC = 1 (crossover)
# At T << T_star: γ_QC << 1 (quantum critical)
# At T >> T_star: γ_QC >> 1 (classical)

print(f"\nγ_QC = (T/T*)^(1/zν)")
print(f"\nAt QCP, coherence scales with QUANTUM critical exponents:")

for i, name in enumerate(names):
    znu = z[i] * nu[i]
    print(f"  {name}: z = {z[i]}, ν = {nu[i]:.2f}, zν = {znu:.2f}")

# =============================================================================
# FERMI LIQUID CROSSOVER
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: Fermi Liquid Crossover")
print("=" * 70)

# T_FL marks where FL behavior sets in
# For T < T_FL: ρ ~ T² (FL)
# For T > T_FL: ρ ~ T^α (NFL)

# T_FL should scale with distance from QCP: T_FL ∝ |g - g_c|^zν

print(f"\nFermi liquid scale T_FL:")
for i, name in enumerate(names):
    ratio = T_star[i] / T_FL[i]
    print(f"  {name}: T_FL = {T_FL[i]:.2f} K, T*/T_FL = {ratio:.0f}")

# Correlation
r1, p1 = stats.pearsonr(z * nu, np.log(T_star / T_FL))
print(f"\nzν vs log(T*/T_FL): r = {r1:.3f}, p = {p1:.4f}")

# =============================================================================
# NON-FERMI LIQUID BEHAVIOR
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: Non-Fermi Liquid Behavior")
print("=" * 70)

# NFL resistivity: ρ ~ T^α where α ≠ 2
# Hertz-Millis theory: α = (d + z - 2) / z for AFM QCP in d dimensions

# Check if α_rho follows Hertz-Millis
d = 3  # Assuming 3D
alpha_HM = (d + z - 2) / z

print(f"\nResistivity exponent α:")
print(f"  System | α_exp | α_HM = (d+z-2)/z | Δα")
print("-" * 50)
for i, name in enumerate(names):
    delta = alpha_rho[i] - alpha_HM[i]
    print(f"  {name[:15]:15} | {alpha_rho[i]:.2f}  | {alpha_HM[i]:.2f}           | {delta:+.2f}")

# Correlation
r2, p2 = stats.pearsonr(alpha_HM, alpha_rho)
print(f"\nα_HM vs α_exp: r = {r2:.3f}, p = {p2:.4f}")

# =============================================================================
# PHASE DIAGRAM TOPOLOGY
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: Quantum Critical Fan")
print("=" * 70)

# The quantum critical fan: T > T_FL for |g - g_c| < δg
# Width scales as: δg ~ T^(1/zν)

print(f"\nQuantum critical fan width ∝ T^(1/zν):")
for i, name in enumerate(names):
    exponent = 1 / (z[i] * nu[i])
    print(f"  {name}: fan width ∝ T^{exponent:.2f}")

# Classify by fan shape
narrow_fan = z * nu > 1.5  # Narrow fan (weak coupling)
wide_fan = z * nu < 1.0    # Wide fan (strong coupling)

print(f"\nFan classification:")
print(f"  Narrow (zν > 1.5): {np.sum(narrow_fan)} systems")
print(f"  Wide (zν < 1.0): {np.sum(wide_fan)} systems")

# =============================================================================
# TYPE COMPARISON
# =============================================================================
print("\n" + "=" * 70)
print("PART 5: System Type Comparison")
print("=" * 70)

type_stats = {}
for t in np.unique(types):
    mask = types == t
    type_stats[t] = {
        'n': np.sum(mask),
        'z': np.mean(z[mask]),
        'nu': np.mean(nu[mask]),
        'alpha': np.mean(alpha_rho[mask]),
        'T_star': np.mean(T_star[mask]),
    }
    print(f"\n{t}:")
    print(f"  n = {type_stats[t]['n']}")
    print(f"  <z> = {type_stats[t]['z']:.1f}")
    print(f"  <ν> = {type_stats[t]['nu']:.2f}")
    print(f"  <α> = {type_stats[t]['alpha']:.2f}")
    print(f"  <T*> = {type_stats[t]['T_star']:.0f} K")

# Compare heavy fermion to itinerant FM
hf_mask = types == "heavy_fermion"
fm_mask = types == "itinerant_FM"

if np.sum(hf_mask) > 1 and np.sum(fm_mask) > 1:
    t_stat, p_val = stats.ttest_ind(z[hf_mask], z[fm_mask])
    print(f"\nHeavy Fermion vs Itinerant FM (z): p = {p_val:.4f}")

# =============================================================================
# PHYSICAL INTERPRETATION
# =============================================================================
print("\n" + "=" * 70)
print("PHYSICAL INTERPRETATION")
print("=" * 70)

print("""
Quantum Criticality from Coherence:

1. QUANTUM CRITICAL POINT (QCP):
   The QCP is where γ = 1 extends to T = 0.
   No energy scale separates quantum from classical.

2. COHERENCE SCALING:
   At QCP: γ_QC = (T/T*)^(1/zν)

   Different from classical: γ_classical = T/T_c

   The exponent 1/(zν) < 1 typically, so:
   - Coherence grows SLOWER with cooling
   - Quantum fluctuations persist longer
   - "Strange metal" behavior emerges

3. FERMI LIQUID CROSSOVER:
   Below T_FL: γ → FL value (coherent quasi-particles)
   Above T_FL: γ grows continuously (NFL)

   T_FL ∝ |g - g_c|^zν (vanishes at QCP)

4. NON-FERMI LIQUID:
   ρ ~ T^α where α = (d + z - 2) / z

   This is COHERENCE-LIMITED scattering:
   - Scattering rate Γ ∝ γ
   - γ ∝ T^(1/zν)
   - ρ ∝ Γ ∝ T^(1/zν) × (other factors)

5. QUANTUM CRITICAL FAN:
   Fan width ∝ T^(1/zν)
   Inside fan: NFL (γ ~ 1, critical)
   Outside fan: FL or ordered (γ ≠ 1)

6. CONNECTION TO SESSIONS:
   - Session #139: Kondo QCP when T_K → 0
   - Session #140: Mott QCP at bandwidth-controlled MIT
   - Session #141: Cuprate QCP under SC dome

   ALL mark γ = 1 boundary extending to T = 0.
""")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Resistivity exponent comparison
ax1 = axes[0, 0]
ax1.scatter(alpha_HM, alpha_rho, c=z, cmap='viridis', s=100)
ax1.plot([0.5, 2], [0.5, 2], 'k--', label='α_exp = α_HM')
ax1.set_xlabel('α_HM = (d+z-2)/z')
ax1.set_ylabel('α_exp')
ax1.set_title(f'Hertz-Millis Test: r = {r2:.3f}')
cbar1 = plt.colorbar(ax1.collections[0], ax=ax1)
cbar1.set_label('z')
for i, name in enumerate(names):
    ax1.annotate(name[:8], (alpha_HM[i], alpha_rho[i]), fontsize=7)

# Plot 2: z vs ν
ax2 = axes[0, 1]
for t in np.unique(types):
    mask = types == t
    ax2.scatter(z[mask], nu[mask], label=t, s=100, alpha=0.7)
ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('z (dynamic exponent)')
ax2.set_ylabel('ν (correlation length exponent)')
ax2.set_title('Critical Exponents by Type')
ax2.legend(fontsize=8)

# Plot 3: Schematic phase diagram
ax3 = axes[1, 0]
g_plot = np.linspace(-1, 1, 200)
T_plot = np.linspace(0, 1, 200)
G, T = np.meshgrid(g_plot, T_plot)

# Phase boundaries
T_FL_line = 0.1 * np.abs(g_plot)**1.5  # FL boundary
T_order = 0.3 * (1 - np.abs(g_plot)**0.5)  # Ordered phase
T_order = np.where(g_plot < 0, T_order, 0)

ax3.fill_between(g_plot, 0, T_order, alpha=0.3, color='blue', label='Ordered')
ax3.fill_between(g_plot, T_FL_line, 1, alpha=0.2, color='green', label='QC fan')
ax3.plot(g_plot, T_FL_line, 'k-', linewidth=2)
ax3.plot(g_plot[g_plot < 0], T_order[g_plot < 0], 'b-', linewidth=2)
ax3.axvline(x=0, color='red', linestyle='--', linewidth=2, label='QCP')
ax3.set_xlabel('g - g_c')
ax3.set_ylabel('T / T*')
ax3.set_title('Quantum Critical Phase Diagram')
ax3.legend()
ax3.set_ylim(0, 1)
ax3.text(0.05, 0.5, 'NFL\nγ~1', fontsize=10)
ax3.text(-0.7, 0.2, 'Ordered\nγ<1', fontsize=10)
ax3.text(0.5, 0.05, 'FL\nγ>1', fontsize=10)

# Plot 4: T* vs T_FL
ax4 = axes[1, 1]
ax4.scatter(T_FL, T_star, c=z*nu, cmap='coolwarm', s=100)
ax4.plot([0, 60], [0, 600], 'k--', alpha=0.5)  # Guide line
ax4.set_xlabel('T_FL (K)')
ax4.set_ylabel('T* (K)')
ax4.set_title('Energy Scales')
cbar4 = plt.colorbar(ax4.collections[0], ax=ax4)
cbar4.set_label('zν')
ax4.set_xscale('log')
ax4.set_yscale('log')
for i, name in enumerate(names):
    ax4.annotate(name[:8], (T_FL[i], T_star[i]), fontsize=7)

plt.suptitle('Session #142: Quantum Criticality and Coherence', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/quantum_criticality_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nVisualization saved to quantum_criticality_coherence.png")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SESSION #142 SUMMARY")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. Quantum critical coherence:
   γ_QC = (T/T*)^(1/zν)
   Different from classical γ = T/T_c

2. Hertz-Millis test:
   α_HM = (d+z-2)/z vs α_exp
   r = {r2:.3f}, p = {p2:.4f}

3. Dynamic exponent z by type:
""")

for t in ['heavy_fermion', 'itinerant_FM', 'cuprate']:
    if t in type_stats:
        print(f"   {t}: <z> = {type_stats[t]['z']:.1f}")

print(f"""
4. Energy scale ratio:
   T*/T_FL ranges from {np.min(T_star/T_FL):.0f} to {np.max(T_star/T_FL):.0f}

5. QCP = γ = 1 at T = 0:
   The quantum critical point is where the
   coherent-incoherent boundary extends to absolute zero.

FRAMEWORK EXTENSION:
QCP is where γ = 1 extends to T = 0.
Quantum critical fan marks region where γ ~ 1.
NFL behavior = coherence-limited scattering.
Exponent α reflects quantum critical dynamics.
""")

print("\nSESSION #142 COMPLETE")
print("=" * 70)
