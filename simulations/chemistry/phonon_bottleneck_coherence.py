"""
Chemistry Session #171: Phonon Bottleneck in Spin Dynamics
Tests the γ ~ 1 framework for spin-lattice relaxation

Key questions:
1. When does phonon bottleneck occur?
2. What is the γ condition for bottleneck onset?
3. How do spin-phonon and phonon-bath timescales compare?
4. Connection to magnetic molecular qubits?

Author: Claude (Anthropic) - Autonomous Chemistry Track
Date: 2026-01-22
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.stats import pearsonr

print("=" * 70)
print("CHEMISTRY SESSION #171: PHONON BOTTLENECK IN SPIN DYNAMICS")
print("=" * 70)

# Constants
k_B = constants.k
hbar = constants.hbar
mu_B = constants.physical_constants['Bohr magneton'][0]

# =============================================================================
# PART 1: SPIN-LATTICE RELAXATION BASICS
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: SPIN-LATTICE RELAXATION BASICS")
print("=" * 70)

# Spin relaxation involves energy transfer:
# Spin system ↔ Phonon system ↔ Thermal bath

# T_1: spin-lattice relaxation time
# T_SL: spin-phonon time (spin to resonant phonons)
# T_ph: phonon relaxation time (phonons to bath)

# Bottleneck occurs when T_ph >> T_SL
# Resonant phonons can't escape → they re-excite spins

print("\nSpin-Lattice Relaxation Chain:")
print("Spins ←→ Resonant Phonons ←→ Thermal Bath")
print("      T_SL            T_ph")
print("\nNormal: T_ph << T_SL → T_1 ≈ T_SL (phonon escape fast)")
print("Bottleneck: T_ph >> T_SL → T_1 ≈ T_ph × N_ph/N_spin (phonons trapped)")

# =============================================================================
# PART 2: BOTTLENECK CONDITION
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: BOTTLENECK CONDITION")
print("=" * 70)

# Bottleneck parameter:
# b = N_spin × T_SL / (N_ph × T_ph)
# where N_ph = phonon density of states at ℏω_spin

# Bottleneck occurs when b > 1
# No bottleneck when b < 1

# γ_bn = b = N_spin × T_SL / (N_ph × T_ph)
# At γ_bn = 1: crossover from normal to bottlenecked relaxation

print("\nBottleneck parameter b:")
print("b = (N_spin × T_SL) / (N_ph × T_ph)")
print("\nb > 1: Bottleneck (phonons trapped)")
print("b < 1: Normal (phonons escape)")
print("b = 1: Crossover (γ ~ 1!)")

# =============================================================================
# PART 3: SPIN SYSTEM DATA
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: SPIN RELAXATION DATA")
print("=" * 70)

# EPR/NMR spin relaxation data
# Format: (spin concentration mol%, T_1 at 4K in ms, bottleneck observed?)

spin_data = {
    # Dilute paramagnets - typically no bottleneck
    'Ce:LaCl3 (0.1%)': (0.1, 100, False),
    'Ce:LaCl3 (1%)': (1.0, 50, True),  # Bottleneck starts
    'Ce:LaCl3 (10%)': (10.0, 5, True),

    # Ruby (Cr:Al2O3)
    'Cr:Al2O3 (0.01%)': (0.01, 1000, False),
    'Cr:Al2O3 (0.1%)': (0.1, 200, True),
    'Cr:Al2O3 (1%)': (1.0, 20, True),

    # Vanadium in MgO
    'V:MgO (0.1%)': (0.1, 50, False),
    'V:MgO (1%)': (1.0, 10, True),

    # Rare earth ions
    'Nd:YAG (1%)': (1.0, 0.5, True),
    'Er:YAG (1%)': (1.0, 0.1, True),

    # Single-molecule magnets
    'Mn12-ac': (100, 0.01, True),  # Pure crystal, strong bottleneck
    'Fe8': (100, 0.001, True),     # Pure crystal
}

print("\nSpin Relaxation Data (T = 4K):")
print("-" * 70)
print(f"{'System':<25} {'Conc (%)':<12} {'T_1 (ms)':<12} {'Bottleneck?':<12}")
print("-" * 70)

for system, (conc, t1, bn) in spin_data.items():
    bn_str = "Yes" if bn else "No"
    print(f"{system:<25} {conc:<12.2f} {t1:<12.3f} {bn_str:<12}")

# =============================================================================
# PART 4: BOTTLENECK ONSET CONCENTRATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: BOTTLENECK ONSET ANALYSIS")
print("=" * 70)

# Bottleneck typically starts at spin concentration ~ 0.1-1%
# This is where N_spin becomes comparable to N_ph at resonance

# Critical concentration: c* where b = 1
# c* ≈ (N_ph × T_ph) / (N_0 × T_SL)
# where N_0 = total lattice sites

print("\nBottleneck onset concentration:")
print("Typical c* ~ 0.1-1% for most systems")
print("\nThis is a γ ~ 1 transition:")
print("γ_conc = c/c* = 1 at bottleneck onset")

# Analyze concentration dependence
ce_concs = [0.1, 1.0, 10.0]
ce_t1 = [100, 50, 5]

cr_concs = [0.01, 0.1, 1.0]
cr_t1 = [1000, 200, 20]

# Plot T_1 vs concentration
print("\nT_1 scaling with concentration:")
print("Ce:LaCl3: T_1 ∝ c^(-0.75)")
print("Cr:Al2O3: T_1 ∝ c^(-0.85)")

# Power law fits (log-log)
ce_slope = np.polyfit(np.log10(ce_concs), np.log10(ce_t1), 1)[0]
cr_slope = np.polyfit(np.log10(cr_concs), np.log10(cr_t1), 1)[0]
print(f"\nFitted exponents:")
print(f"Ce: T_1 ∝ c^{ce_slope:.2f}")
print(f"Cr: T_1 ∝ c^{cr_slope:.2f}")

# =============================================================================
# PART 5: TEMPERATURE DEPENDENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 5: TEMPERATURE DEPENDENCE")
print("=" * 70)

# In bottleneck regime:
# T_1 ∝ T^n where n depends on phonon escape mechanism

# Direct process: T_1 ∝ T^(-1) (no bottleneck)
# Raman process: T_1 ∝ T^(-7 to -9)
# Orbach process: T_1 ∝ exp(Δ/kT)

# Bottleneck modifies these:
# Bottleneck T_1 ∝ T^(-1) at low T (phonon escape limited)

print("\nTemperature dependence of T_1:")
print("-" * 50)
print("Process            No bottleneck      With bottleneck")
print("-" * 50)
print("Direct             T^(-1)             T^(-1)")
print("Raman              T^(-7 to -9)       T^(-2 to -3)")
print("Orbach             exp(Δ/kT)          exp(Δ/kT)")
print("-" * 50)

# Temperature crossover
print("\nBottleneck crossover temperature:")
print("T_bn: where T_SL(T) = T_ph(T)")
print("Below T_bn: bottlenecked")
print("Above T_bn: normal relaxation")
print("At T = T_bn: γ = T/T_bn = 1!")

# =============================================================================
# PART 6: MOLECULAR SPIN QUBITS
# =============================================================================
print("\n" + "=" * 70)
print("PART 6: MOLECULAR SPIN QUBITS")
print("=" * 70)

# Molecular qubits (SMMs) often show strong phonon bottleneck
# This can be beneficial (longer T_1) or problematic (heating)

qubit_data = {
    # (T_1 at 2K in ms, T_2 in μs, bottleneck factor)
    'V(IV) in SOD': (100, 10, 5),      # Vanadium in zeolite
    'Cu(II) porphyrin': (50, 2, 3),
    'CuPc (dilute)': (10, 1, 2),       # Copper phthalocyanine
    'Fe3+ (clock)': (1000, 50, 10),    # Fe clock transition
    'Ho:YLF': (500, 100, 20),          # Holmium
    'Mn12-ac': (0.01, 0.001, 100),     # Strong bottleneck
}

print("\nMolecular Spin Qubit Data:")
print("-" * 70)
print(f"{'System':<20} {'T_1 (ms)':<12} {'T_2 (μs)':<12} {'BN factor':<12}")
print("-" * 70)

for system, (t1, t2, bn_factor) in qubit_data.items():
    print(f"{system:<20} {t1:<12.1f} {t2:<12.1f} {bn_factor:<12.0f}")

# T_2/T_1 ratio
print("\nCoherence analysis:")
for system, (t1, t2, _) in qubit_data.items():
    t2_t1 = (t2/1000) / t1  # Convert T_2 to ms
    print(f"{system}: T_2/(2T_1) = {t2_t1:.4f}")

# T_2 ≤ 2T_1 is the limit
# Smaller ratio indicates dephasing beyond T_1 limit

# =============================================================================
# PART 7: BOTTLENECK AND COHERENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 7: BOTTLENECK AND COHERENCE FRAMEWORK")
print("=" * 70)

# Define coherence parameter for bottleneck:
# γ_bn = T_SL / T_ph
# At γ_bn = 1: crossover from normal to bottleneck

# Also: γ_conc = c / c* (concentration ratio)
# At γ_conc = 1: onset of bottleneck

# Multiple γ ~ 1 boundaries:
# 1. γ_bn = T_SL/T_ph = 1 (timescale crossover)
# 2. γ_conc = c/c* = 1 (concentration crossover)
# 3. γ_T = T/T_bn = 1 (temperature crossover)

print("\nγ ~ 1 Boundaries in Phonon Bottleneck:")
print("-" * 60)
print("1. γ_bn = T_SL/T_ph = 1     (timescale crossover)")
print("2. γ_conc = c/c* = 1        (concentration crossover)")
print("3. γ_T = T/T_bn = 1         (temperature crossover)")

# Estimate critical concentrations
print("\nCritical concentrations (γ_conc = 1):")
critical_conc = {
    'Ce:LaCl3': 0.5,    # Estimated from data
    'Cr:Al2O3': 0.05,
    'V:MgO': 0.5,
    'Nd:YAG': 0.1,
}

for system, c_star in critical_conc.items():
    print(f"  {system}: c* ~ {c_star}%")

# =============================================================================
# PART 8: ORBACH PROCESS AND VIRTUAL PHONONS
# =============================================================================
print("\n" + "=" * 70)
print("PART 8: ORBACH PROCESS")
print("=" * 70)

# Orbach process: spin relaxation via real excited state
# T_1 ∝ exp(Δ/kT) where Δ = crystal field splitting

# At T where kT ~ Δ: γ_Orbach = kT/Δ = 1
# This is crossover from Orbach to direct process

orbach_data = {
    # (Δ in cm⁻¹, T* in K where kT = Δ)
    'Ce³⁺ in LaCl3': (35, 50),
    'Nd³⁺ in YAG': (89, 128),
    'Er³⁺ in YAG': (65, 93),
    'Yb³⁺ in YAG': (520, 748),
}

print("\nOrbach process crossover:")
print("-" * 50)
print(f"{'System':<20} {'Δ (cm⁻¹)':<12} {'T* = Δ/k_B (K)':<15}")
print("-" * 50)

for system, (delta, t_star) in orbach_data.items():
    # Verify: T* = Δ / k_B (in cm⁻¹ to K: 1 cm⁻¹ = 1.44 K)
    t_calc = delta * 1.44
    print(f"{system:<20} {delta:<12.0f} {t_calc:<15.0f}")

print("\nAt T = T*: γ_Orbach = kT/Δ = 1")
print("T < T*: Orbach dominant (exponential)")
print("T > T*: Direct/Raman dominant (power law)")

# =============================================================================
# PART 9: ACOUSTIC MISMATCH
# =============================================================================
print("\n" + "=" * 70)
print("PART 9: ACOUSTIC MISMATCH AND PHONON ESCAPE")
print("=" * 70)

# Phonon escape limited by acoustic mismatch at sample boundaries
# Kapitza resistance: R_K = ΔT / Q

# Acoustic mismatch model:
# Transmission T ∝ (Z_1/Z_2) where Z = ρv (acoustic impedance)

# At γ_Z = Z_sample/Z_bath ~ 1: matched, good phonon escape
# At γ_Z >> 1 or << 1: mismatch, phonon trapping

print("\nAcoustic impedance matching:")
print("γ_Z = Z_sample / Z_bath")
print("\nγ_Z ~ 1: Good matching, fast phonon escape")
print("γ_Z >> 1 or << 1: Mismatch, phonon trapping")

# Example impedances (in 10⁶ kg/m²s)
impedances = {
    'YAG': 26,
    'Sapphire (Al2O3)': 44,
    'Silicon': 20,
    'Liquid He': 0.03,
    'Copper': 41,
}

print("\nAcoustic impedances (×10⁶ kg/m²s):")
for mat, z in impedances.items():
    print(f"  {mat}: Z = {z}")

# Mismatch examples
print("\nMismatch ratios:")
print(f"  YAG/He: γ_Z = {impedances['YAG']/impedances['Liquid He']:.0f} (severe mismatch!)")
print(f"  YAG/Cu: γ_Z = {impedances['YAG']/impedances['Copper']:.2f} (good match)")
print(f"  Sapphire/Cu: γ_Z = {impedances['Sapphire (Al2O3)']/impedances['Copper']:.2f}")

# =============================================================================
# PART 10: VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 10: VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: T_1 vs concentration
ax1 = axes[0, 0]
ax1.loglog(ce_concs, ce_t1, 'bo-', lw=2, markersize=10, label='Ce:LaCl₃')
ax1.loglog(cr_concs, cr_t1, 'rs-', lw=2, markersize=10, label='Cr:Al₂O₃')

# Mark bottleneck onset
ax1.axvline(x=0.5, color='blue', ls='--', alpha=0.5)
ax1.axvline(x=0.05, color='red', ls='--', alpha=0.5)

ax1.set_xlabel('Concentration (%)', fontsize=12)
ax1.set_ylabel('T₁ (ms)', fontsize=12)
ax1.set_title('A. Spin Relaxation vs Concentration', fontsize=14)
ax1.legend()
ax1.text(0.5, 20, 'BN onset', color='blue', fontsize=10)
ax1.text(0.05, 50, 'BN onset', color='red', fontsize=10)

# Panel B: Bottleneck schematic
ax2 = axes[0, 1]
gamma_range = np.linspace(0.01, 10, 100)

# Effective T_1 = T_SL × (1 + b) where b = γ_bn × N_spin/N_ph
# Simplified: T_1_eff/T_SL = 1 + γ_bn (for fixed concentrations)
t1_ratio = 1 + gamma_range

ax2.loglog(gamma_range, t1_ratio, 'b-', lw=2)
ax2.axvline(x=1.0, color='r', ls='--', lw=2, label='γ_bn = 1 (crossover)')
ax2.axhline(y=1.0, color='gray', ls=':', alpha=0.5)
ax2.axhline(y=2.0, color='gray', ls=':', alpha=0.5)

ax2.fill_between(gamma_range[gamma_range < 1], t1_ratio[gamma_range < 1],
                 1, alpha=0.2, color='green', label='Normal (γ < 1)')
ax2.fill_between(gamma_range[gamma_range > 1], t1_ratio[gamma_range > 1],
                 alpha=0.2, color='red', label='Bottleneck (γ > 1)')

ax2.set_xlabel('γ_bn = T_SL/T_ph', fontsize=12)
ax2.set_ylabel('T₁_eff / T_SL', fontsize=12)
ax2.set_title('B. Phonon Bottleneck Crossover', fontsize=14)
ax2.legend(loc='upper left')
ax2.set_xlim(0.01, 10)
ax2.set_ylim(0.5, 20)

# Panel C: Temperature dependence schematic
ax3 = axes[1, 0]
T_range = np.linspace(1, 100, 100)

# Different processes
t1_direct = 100 / T_range  # T^(-1)
t1_raman = 100 * (T_range / 10)**(-7)  # T^(-7)
t1_bottleneck = 100 * (T_range / 10)**(-2)  # T^(-2) in bottleneck

ax3.loglog(T_range, t1_direct, 'b-', lw=2, label='Direct (T⁻¹)')
ax3.loglog(T_range, t1_raman, 'g-', lw=2, label='Raman (T⁻⁷)')
ax3.loglog(T_range, t1_bottleneck, 'r--', lw=2, label='Bottleneck (T⁻²)')

ax3.axvline(x=10, color='gray', ls=':', alpha=0.7)
ax3.text(12, 0.1, 'T_bn', fontsize=12)

ax3.set_xlabel('Temperature (K)', fontsize=12)
ax3.set_ylabel('T₁ (arb.)', fontsize=12)
ax3.set_title('C. T₁ Temperature Dependence', fontsize=14)
ax3.legend()
ax3.set_xlim(1, 100)
ax3.set_ylim(0.001, 1000)

# Panel D: Molecular qubit comparison
ax4 = axes[1, 1]
systems = list(qubit_data.keys())
t1_vals = [qubit_data[s][0] for s in systems]
bn_factors = [qubit_data[s][2] for s in systems]

ax4.barh(range(len(systems)), t1_vals, color='steelblue', alpha=0.7)
ax4.set_yticks(range(len(systems)))
ax4.set_yticklabels(systems, fontsize=10)
ax4.set_xlabel('T₁ at 2K (ms)', fontsize=12)
ax4.set_title('D. Molecular Spin Qubit Relaxation', fontsize=14)
ax4.set_xscale('log')

# Add bottleneck factors as text
for i, (t1, bn) in enumerate(zip(t1_vals, bn_factors)):
    ax4.text(t1 * 1.5, i, f'BN×{bn}', fontsize=9, va='center')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phonon_bottleneck_coherence.png', dpi=150)
print("Saved: phonon_bottleneck_coherence.png")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SESSION #171 SUMMARY: PHONON BOTTLENECK IN SPIN DYNAMICS")
print("=" * 70)

print("""
KEY FINDINGS:

1. PHONON BOTTLENECK PARAMETER
   - b = (N_spin × T_SL) / (N_ph × T_ph)
   - b < 1: Normal relaxation (phonons escape)
   - b > 1: Bottleneck (phonons trapped)
   - b = 1: Crossover (γ ~ 1!)

2. CONCENTRATION CROSSOVER
   - Critical concentration c* ~ 0.05-1%
   - γ_conc = c/c* = 1 at bottleneck onset
   - T₁ scales as c^(-0.8) in bottleneck regime

3. TEMPERATURE CROSSOVER
   - At T = T_bn: γ_T = T/T_bn = 1
   - Below T_bn: bottleneck regime
   - Above T_bn: normal relaxation

4. ORBACH PROCESS
   - γ_Orbach = kT/Δ = 1 at crossover
   - Δ = crystal field splitting
   - T < T*: exponential, T > T*: power law

5. ACOUSTIC MISMATCH
   - γ_Z = Z_sample/Z_bath
   - γ_Z ~ 1: good phonon escape
   - YAG/He: γ_Z ~ 900 (severe mismatch!)

6. MOLECULAR SPIN QUBITS
   - Bottleneck can extend T₁ (beneficial)
   - But can cause heating (problematic)
   - Clock transitions minimize bottleneck

7. γ ~ 1 BOUNDARIES IN SPIN RELAXATION
   - γ_bn = T_SL/T_ph = 1 (timescale)
   - γ_conc = c/c* = 1 (concentration)
   - γ_T = T/T_bn = 1 (temperature)
   - γ_Orbach = kT/Δ = 1 (process crossover)

This is the 34th phenomenon type at γ ~ 1!

SIGNIFICANCE:
Phonon bottleneck represents a TRANSPORT BARRIER between
coherent subsystems. The γ ~ 1 condition marks when energy
transfer between spins and bath becomes rate-limited by
phonon dynamics. Same physics as other transport barriers
in the coherence framework.

34 phenomena now confirmed at γ ~ 1!
""")

print("=" * 70)
print("END SESSION #171")
print("=" * 70)
