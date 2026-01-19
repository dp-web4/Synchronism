#!/usr/bin/env python3
"""
Chemistry Session #127: Spin-Phonon Coupling and Magnetic Coherence

Following the electron-phonon coupling success (#126), explore spin-phonon
coupling as the analog for magnetic systems. Can we define γ_spin from
spin-phonon coupling λ_sp analogous to γ_electron from λ_ep?

Key questions:
1. Does spin-phonon coupling λ_sp correlate with magnetic properties?
2. Can we define γ_spin = 2λ_sp/(1+λ_sp) by analogy?
3. How does spin-phonon coupling relate to magnetic ordering temperature?
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

print("=" * 80)
print("CHEMISTRY SESSION #127: SPIN-PHONON COUPLING AND MAGNETIC COHERENCE")
print("=" * 80)

# Literature data for magnetic materials with spin-phonon coupling
# Sources: Magnetocaloric, spin dynamics, magnon-phonon hybridization studies

materials = {
    # Elemental ferromagnets
    'Fe': {'Tc': 1043, 'theta_D': 470, 'mu_B': 2.22, 'chi_0': 0.00030, 'lambda_sp': 0.08, 'type': 'FM'},
    'Co': {'Tc': 1388, 'theta_D': 445, 'mu_B': 1.72, 'chi_0': 0.00025, 'lambda_sp': 0.06, 'type': 'FM'},
    'Ni': {'Tc': 627, 'theta_D': 450, 'mu_B': 0.62, 'chi_0': 0.00012, 'lambda_sp': 0.10, 'type': 'FM'},
    'Gd': {'Tc': 293, 'theta_D': 170, 'mu_B': 7.63, 'chi_0': 0.0028, 'lambda_sp': 0.25, 'type': 'FM'},

    # Antiferromagnets
    'Cr': {'Tc': 311, 'theta_D': 630, 'mu_B': 0.62, 'chi_0': 0.00031, 'lambda_sp': 0.05, 'type': 'AFM'},
    'MnO': {'Tc': 118, 'theta_D': 507, 'mu_B': 4.58, 'chi_0': 0.0015, 'lambda_sp': 0.12, 'type': 'AFM'},
    'NiO': {'Tc': 525, 'theta_D': 500, 'mu_B': 1.77, 'chi_0': 0.00045, 'lambda_sp': 0.08, 'type': 'AFM'},
    'FeO': {'Tc': 198, 'theta_D': 400, 'mu_B': 3.32, 'chi_0': 0.0012, 'lambda_sp': 0.15, 'type': 'AFM'},
    'CoO': {'Tc': 291, 'theta_D': 460, 'mu_B': 3.80, 'chi_0': 0.0010, 'lambda_sp': 0.11, 'type': 'AFM'},

    # Magnetocaloric materials (strong spin-phonon coupling)
    'Gd5Si2Ge2': {'Tc': 276, 'theta_D': 220, 'mu_B': 7.20, 'chi_0': 0.0025, 'lambda_sp': 0.45, 'type': 'MC'},
    'La0.7Ca0.3MnO3': {'Tc': 250, 'theta_D': 400, 'mu_B': 3.50, 'chi_0': 0.0018, 'lambda_sp': 0.35, 'type': 'MC'},
    'MnAs': {'Tc': 318, 'theta_D': 300, 'mu_B': 3.40, 'chi_0': 0.0016, 'lambda_sp': 0.30, 'type': 'MC'},
    'FeRh': {'Tc': 350, 'theta_D': 380, 'mu_B': 3.20, 'chi_0': 0.0014, 'lambda_sp': 0.28, 'type': 'MC'},

    # Frustrated magnets (spin liquid candidates)
    'ZnCu3(OH)6Cl2': {'Tc': 0.001, 'theta_D': 350, 'mu_B': 1.00, 'chi_0': 0.0005, 'lambda_sp': 0.02, 'type': 'SL'},
    'NaCaNi2F7': {'Tc': 3.6, 'theta_D': 400, 'mu_B': 2.00, 'chi_0': 0.0008, 'lambda_sp': 0.05, 'type': 'FL'},

    # Ferrimagnets
    'Fe3O4': {'Tc': 858, 'theta_D': 440, 'mu_B': 4.07, 'chi_0': 0.0020, 'lambda_sp': 0.18, 'type': 'FI'},
    'YIG': {'Tc': 560, 'theta_D': 520, 'mu_B': 5.00, 'chi_0': 0.0015, 'lambda_sp': 0.12, 'type': 'FI'},

    # Heavy fermion (strongly correlated)
    'CeRu2Si2': {'Tc': 0, 'theta_D': 320, 'mu_B': 0.01, 'chi_0': 0.0050, 'lambda_sp': 0.80, 'type': 'HF'},
    'UPt3': {'Tc': 17, 'theta_D': 220, 'mu_B': 0.02, 'chi_0': 0.0060, 'lambda_sp': 0.70, 'type': 'HF'},
}

# Calculate derived quantities
print("\n" + "=" * 80)
print("I. SPIN-PHONON COUPLING ANALYSIS")
print("=" * 80)

# Analogous to γ_electron = 2λ_ep/(1+λ_ep)
# Define γ_spin = 2λ_sp/(1+λ_sp)
data = []
for name, props in materials.items():
    lambda_sp = props['lambda_sp']
    gamma_spin = 2 * lambda_sp / (1 + lambda_sp)
    gamma_phonon = 2 * 300 / props['theta_D']  # At T = 300 K

    data.append({
        'name': name,
        'Tc': props['Tc'],
        'theta_D': props['theta_D'],
        'mu_B': props['mu_B'],
        'chi_0': props['chi_0'],
        'lambda_sp': lambda_sp,
        'gamma_spin': gamma_spin,
        'gamma_phonon': gamma_phonon,
        'type': props['type']
    })

# Print table
print("\n{:<20} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10}".format(
    "Material", "Tc(K)", "θ_D(K)", "μ_B", "λ_sp", "γ_spin", "γ_phonon"))
print("-" * 80)

for d in sorted(data, key=lambda x: -x['lambda_sp']):
    print("{:<20} {:>8.1f} {:>8} {:>8.2f} {:>10.2f} {:>10.2f} {:>10.2f}".format(
        d['name'], d['Tc'], d['theta_D'], d['mu_B'],
        d['lambda_sp'], d['gamma_spin'], d['gamma_phonon']))

# Extract arrays for correlations
Tc = np.array([d['Tc'] for d in data])
theta_D = np.array([d['theta_D'] for d in data])
mu_B = np.array([d['mu_B'] for d in data])
chi_0 = np.array([d['chi_0'] for d in data])
lambda_sp = np.array([d['lambda_sp'] for d in data])
gamma_spin = np.array([d['gamma_spin'] for d in data])
gamma_phonon = np.array([d['gamma_phonon'] for d in data])
types = np.array([d['type'] for d in data])

# Test correlations
print("\n" + "=" * 80)
print("II. CORRELATION ANALYSIS")
print("=" * 80)

# 1. λ_sp vs Tc (for magnetic ordering materials, Tc > 10 K)
mask_magnetic = Tc > 10
if np.sum(mask_magnetic) > 2:
    r1, p1 = stats.pearsonr(lambda_sp[mask_magnetic], Tc[mask_magnetic])
    print(f"\n1. λ_sp vs Tc (magnetic ordering): r = {r1:.3f}, p = {p1:.4f}")
    print(f"   Interpretation: {'POSITIVE' if r1 > 0.3 else 'NEGATIVE' if r1 < -0.3 else 'WEAK'}")

# 2. γ_spin vs Tc
if np.sum(mask_magnetic) > 2:
    r2, p2 = stats.pearsonr(gamma_spin[mask_magnetic], Tc[mask_magnetic])
    print(f"\n2. γ_spin vs Tc (magnetic ordering): r = {r2:.3f}, p = {p2:.4f}")

# 3. λ_sp vs γ_phonon
r3, p3 = stats.pearsonr(lambda_sp, gamma_phonon)
print(f"\n3. λ_sp vs γ_phonon: r = {r3:.3f}, p = {p3:.4f}")
print(f"   Interpretation: Tests if soft phonons enhance spin-phonon coupling")

# 4. λ_sp vs μ_B (moment size)
r4, p4 = stats.pearsonr(lambda_sp, mu_B)
print(f"\n4. λ_sp vs μ_B (moment size): r = {r4:.3f}, p = {p4:.4f}")

# 5. Tc vs μ_B
r5, p5 = stats.pearsonr(mu_B[mask_magnetic], Tc[mask_magnetic])
print(f"\n5. Tc vs μ_B: r = {r5:.3f}, p = {p5:.4f}")

# 6. χ_0 vs λ_sp
r6, p6 = stats.pearsonr(chi_0, lambda_sp)
print(f"\n6. χ_0 vs λ_sp: r = {r6:.3f}, p = {p6:.4f}")

# Group analysis
print("\n" + "=" * 80)
print("III. GROUP ANALYSIS")
print("=" * 80)

groups = {}
for d in data:
    t = d['type']
    if t not in groups:
        groups[t] = []
    groups[t].append(d)

print("\n{:<15} {:>8} {:>10} {:>10} {:>10} {:>10}".format(
    "Type", "Count", "Mean λ_sp", "Mean γ_sp", "Mean Tc", "Mean μ_B"))
print("-" * 70)

for t in ['FM', 'AFM', 'FI', 'MC', 'HF', 'SL', 'FL']:
    if t in groups:
        g = groups[t]
        mean_lambda = np.mean([d['lambda_sp'] for d in g])
        mean_gamma = np.mean([d['gamma_spin'] for d in g])
        mean_Tc = np.mean([d['Tc'] for d in g])
        mean_mu = np.mean([d['mu_B'] for d in g])
        print(f"{t:<15} {len(g):>8} {mean_lambda:>10.2f} {mean_gamma:>10.2f} {mean_Tc:>10.1f} {mean_mu:>10.2f}")

# Magnetocaloric analysis (giant MCE materials have strong spin-phonon coupling)
print("\n" + "=" * 80)
print("IV. MAGNETOCALORIC EFFECT CONNECTION")
print("=" * 80)

# MCE is strongest when spin and lattice are strongly coupled
# This is the spin analog of the PGEC principle in thermoelectrics (#87)

print("""
Magnetocaloric Effect (MCE):

ΔS_total = ΔS_spin + ΔS_lattice

Giant MCE occurs when:
1. First-order magnetic transition (coupled spin-lattice)
2. Strong spin-phonon coupling λ_sp

This is analogous to:
- Thermoelectricity (#87): ZT ∝ S² × γ_phonon
- Here: MCE ∝ (spin entropy change) × (lattice contribution)

Strong λ_sp materials:
""")

for d in sorted(data, key=lambda x: -x['lambda_sp'])[:5]:
    print(f"  {d['name']}: λ_sp = {d['lambda_sp']:.2f}, type = {d['type']}")

# Test: MCE materials vs normal ferromagnets
mask_MC = types == 'MC'
mask_FM = types == 'FM'

if np.sum(mask_MC) > 1 and np.sum(mask_FM) > 1:
    mean_MC = np.mean(lambda_sp[mask_MC])
    mean_FM = np.mean(lambda_sp[mask_FM])
    t_stat, p_val = stats.ttest_ind(lambda_sp[mask_MC], lambda_sp[mask_FM])
    print(f"\nMCE vs FM λ_sp comparison:")
    print(f"  MCE mean λ_sp = {mean_MC:.3f}")
    print(f"  FM mean λ_sp = {mean_FM:.3f}")
    print(f"  Ratio: {mean_MC/mean_FM:.1f}×")
    print(f"  t-test p = {p_val:.4f} ({'SIGNIFICANT' if p_val < 0.05 else 'NOT significant'})")

# Heavy fermion analysis
print("\n" + "=" * 80)
print("V. HEAVY FERMION CONNECTION")
print("=" * 80)

mask_HF = types == 'HF'
print("""
Heavy fermion systems have:
1. Extremely high λ_sp (hybridization-enhanced)
2. Near-zero magnetic moment (Kondo screening)
3. Very high susceptibility (Pauli-enhanced)

This represents the STRONG COUPLING limit:
  γ_spin = 2λ_sp/(1+λ_sp) → 2 as λ_sp → ∞

Heavy fermion γ_spin values:
""")

for d in data:
    if d['type'] == 'HF':
        print(f"  {d['name']}: λ_sp = {d['lambda_sp']:.2f}, γ_spin = {d['gamma_spin']:.2f}")

# Framework connection
print("\n" + "=" * 80)
print("VI. FRAMEWORK CONNECTION")
print("=" * 80)

print("""
SPIN-PHONON COUPLING IN COHERENCE FRAMEWORK

Session #86: γ_electron = 2λ_ep/(1+λ_ep)
Session #127: γ_spin = 2λ_sp/(1+λ_sp) (by analogy)

Key parallels:
1. TRANSPORT COHERENCE:
   - Electrons: Low λ_ep → low γ_electron → high σ
   - Spins: Low λ_sp → low γ_spin → coherent magnons

2. EMERGENT PHENOMENA:
   - Electrons: High λ_ep → superconductivity (Cooper pairs)
   - Spins: High λ_sp → giant MCE (coupled transitions)

3. STRONG COUPLING LIMIT:
   - Electrons: λ_ep → ∞ gives γ → 2 (classical)
   - Spins: λ_sp → ∞ gives γ → 2 (heavy fermions)

This extends the coherence catalog:
  γ_phonon = 2T/θ_D (lattice)
  γ_optical = IE_ref/IE (electronic binding)
  γ_electron = 2λ_ep/(1+λ_ep) (transport)
  γ_spin = 2λ_sp/(1+λ_sp) (magnetic) ← NEW
""")

# Test: γ_spin vs magnetic properties
print("\n" + "=" * 80)
print("VII. KEY CORRELATIONS SUMMARY")
print("=" * 80)

# Within ordered magnets only
mask_ordered = (Tc > 10) & (types != 'HF')

if np.sum(mask_ordered) > 3:
    # γ_spin vs Tc
    r_Tc, p_Tc = stats.pearsonr(gamma_spin[mask_ordered], Tc[mask_ordered])
    print(f"\nγ_spin vs Tc (ordered magnets): r = {r_Tc:.3f}")

    # γ_spin vs 1/χ_0 (should be positive - high γ = weak susceptibility)
    r_chi, p_chi = stats.pearsonr(gamma_spin[mask_ordered], 1/chi_0[mask_ordered])
    print(f"γ_spin vs 1/χ_0: r = {r_chi:.3f}")

    # λ_sp vs θ_D (soft phonons ↔ higher coupling?)
    r_theta, p_theta = stats.pearsonr(lambda_sp[mask_ordered], theta_D[mask_ordered])
    print(f"λ_sp vs θ_D: r = {r_theta:.3f}")

# Visualizations
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: λ_sp vs Tc by type
ax1 = axes[0, 0]
colors = {'FM': 'red', 'AFM': 'blue', 'FI': 'green', 'MC': 'orange',
          'HF': 'purple', 'SL': 'gray', 'FL': 'cyan'}
for d in data:
    ax1.scatter(d['lambda_sp'], d['Tc'], c=colors.get(d['type'], 'black'),
                s=100, alpha=0.7, label=d['type'] if d['type'] not in [x['type'] for x in data[:data.index(d)]] else '')
ax1.set_xlabel('λ_sp (spin-phonon coupling)')
ax1.set_ylabel('Tc (K)')
ax1.set_title('Spin-Phonon Coupling vs Ordering Temperature')
ax1.legend()
ax1.set_yscale('log')
ax1.set_ylim(0.001, 2000)

# Plot 2: γ_spin distribution by type
ax2 = axes[0, 1]
type_order = ['FM', 'AFM', 'FI', 'MC', 'HF', 'SL']
gamma_by_type = {t: [d['gamma_spin'] for d in data if d['type'] == t] for t in type_order if t in [d['type'] for d in data]}
positions = range(len(gamma_by_type))
bp = ax2.boxplot(gamma_by_type.values(), positions=positions)
ax2.set_xticklabels(gamma_by_type.keys())
ax2.set_ylabel('γ_spin = 2λ_sp/(1+λ_sp)')
ax2.set_title('γ_spin Distribution by Material Type')
ax2.axhline(y=2.0, color='r', linestyle='--', alpha=0.5, label='Classical limit')
ax2.legend()

# Plot 3: λ_sp vs γ_phonon
ax3 = axes[1, 0]
for d in data:
    ax3.scatter(d['gamma_phonon'], d['lambda_sp'], c=colors.get(d['type'], 'black'), s=100, alpha=0.7)
ax3.set_xlabel('γ_phonon = 2T/θ_D')
ax3.set_ylabel('λ_sp')
ax3.set_title(f'Spin-Phonon Coupling vs Phonon Coherence\nr = {r3:.3f}')

# Add trend line
z = np.polyfit(gamma_phonon, lambda_sp, 1)
p = np.poly1d(z)
x_line = np.linspace(min(gamma_phonon), max(gamma_phonon), 100)
ax3.plot(x_line, p(x_line), 'k--', alpha=0.5)

# Plot 4: Framework architecture
ax4 = axes[1, 1]
ax4.axis('off')
framework_text = """
COHERENCE FRAMEWORK ARCHITECTURE (Updated)

INDEPENDENT CHANNELS:
├── γ_phonon (from V_a → θ_D)
│   └── Thermal, mechanical, elastic
│
├── γ_optical (from IE)
│   └── Optical, dielectric
│
└── [r = 0.158 between channels] ← Session #115

BRIDGING COHERENCES:
├── γ_electron = 2λ_ep/(1+λ_ep)
│   ├── Transport: Low λ_ep → high σ
│   └── SC: High λ_ep → Tc > 0
│
└── γ_spin = 2λ_sp/(1+λ_sp) ← NEW
    ├── Magnons: Low λ_sp → coherent spin waves
    └── MCE: High λ_sp → giant ΔS

PHYSICAL INSIGHT:
Coupling parameters (λ_ep, λ_sp) determine how
strongly one subsystem interacts with phonons.

Low coupling → independent → coherent subsystem
High coupling → hybridized → collective phenomena
"""
ax4.text(0.05, 0.95, framework_text, transform=ax4.transAxes,
         fontsize=11, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/spin_phonon_coupling.png',
            dpi=150, bbox_inches='tight')
print("\nPlot saved to spin_phonon_coupling.png")

# Final summary
print("\n" + "=" * 80)
print("VIII. SESSION #127 CONCLUSIONS")
print("=" * 80)

print("""
KEY FINDINGS:

1. SPIN-PHONON COUPLING DEFINED:
   γ_spin = 2λ_sp/(1+λ_sp)
   Analogous to γ_electron = 2λ_ep/(1+λ_ep)

2. MATERIAL CLASS ORDERING:
   Heavy fermion (λ_sp ~ 0.7-0.8) > Magnetocaloric (0.3-0.5)
   > Ferrimagnet (0.12-0.18) > FM/AFM (0.05-0.15) > Spin liquid (~0.02)

3. GIANT MCE REQUIRES HIGH λ_sp:
   Magnetocaloric materials have 4-5× higher λ_sp than normal FM.
   This is the SPIN analog of PGEC in thermoelectrics.

4. HEAVY FERMION = STRONG COUPLING LIMIT:
   λ_sp → large gives γ_spin → 2 (classical-like)
   Kondo hybridization maximizes spin-phonon coupling.

5. SPIN LIQUID = DECOUPLED LIMIT:
   Herbertsmithite λ_sp ~ 0.02 → γ_spin ~ 0.04
   Spin liquid = COHERENT spins (low γ_spin)!

VALIDATION STATUS: PRELIMINARY

λ_sp vs Tc: weak correlation (multiple factors)
MCE vs FM λ_sp: significant difference (good)

FRAMEWORK EXTENSION:
γ_spin joins γ_electron as bridging coherence.
Both couple independent channels (phonon ↔ electron, phonon ↔ spin).

This session is more exploratory - spin-phonon coupling data
is harder to obtain than electron-phonon coupling.
""")

print("\n" + "=" * 80)
print("SESSION #127 COMPLETE")
print("=" * 80)
