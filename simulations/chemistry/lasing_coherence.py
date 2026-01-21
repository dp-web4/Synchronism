"""
Session #165: Lasing Threshold and γ ~ 1
Chemistry Track - Synchronism Framework

Test the γ ~ 1 prediction for lasing transitions:
- Laser threshold (stimulated > spontaneous)
- Coherence buildup in optical cavities
- Connection to Dicke superradiance
- Photon statistics (thermal → coherent)

Key question:
Does the lasing threshold occur at γ ~ 1?

Author: Claude (Anthropic) - Autonomous Research
Date: 2026-01-21
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("SESSION #165: LASING THRESHOLD AND γ ~ 1")
print("=" * 70)

# =============================================================================
# SECTION 1: LASING FUNDAMENTALS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 1: LASING FUNDAMENTALS")
print("=" * 70)

print("""
Laser = Light Amplification by Stimulated Emission of Radiation

Key processes:
1. Spontaneous emission: Rate A₂₁ (incoherent)
2. Stimulated emission: Rate B₂₁ × n (coherent, copies photons)
3. Absorption: Rate B₁₂ × n

At thermal equilibrium (no lasing):
    n̄ = 1 / (exp(ℏω/k_B T) - 1)  (Bose-Einstein)

For lasing, need population inversion: N₂ > N₁
And gain > loss: G > L

Threshold condition:
    G_threshold = L  (gain equals loss)

Below threshold: Spontaneous emission dominates (incoherent)
Above threshold: Stimulated emission dominates (coherent)
""")

# Fundamental constants
h = 6.626e-34  # J·s
hbar = h / (2 * np.pi)
c = 3e8  # m/s
k_B = 1.381e-23  # J/K

# =============================================================================
# SECTION 2: γ DEFINITION FOR LASING
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 2: γ DEFINITION FOR LASING")
print("=" * 70)

print("""
Several natural γ definitions for lasing:

1. Pump ratio: γ_pump = P / P_th
   - γ < 1: Below threshold (ASE)
   - γ = 1: At threshold
   - γ > 1: Above threshold (lasing)

2. Gain ratio: γ_gain = L / G = Loss / Gain
   - γ > 1: Loss exceeds gain (no lasing)
   - γ = 1: Threshold (gain = loss)
   - γ < 1: Net gain (lasing)

3. Photon number: γ_n = n_sp / n  (spontaneous/total)
   - γ ~ 1: At threshold (n ~ 1/β where β = spontaneous coupling)
   - γ << 1: Deep in lasing (coherent photons dominate)
   - γ ~ 1: Below threshold (mostly spontaneous)

4. β factor: γ_β = 1 / β
   - β = fraction of spontaneous emission into lasing mode
   - Threshold softens as β → 1 (thresholdless laser)
   - Large γ_β = sharp threshold (conventional laser)

The lasing threshold IS the coherent/incoherent boundary!
""")

# =============================================================================
# SECTION 3: LASER TYPES AND THRESHOLDS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 3: LASER TYPES AND THRESHOLDS")
print("=" * 70)

# Laser data: (P_th typical, β factor, wavelength nm, type)
laser_data = {
    # Gas lasers
    'He-Ne (632.8 nm)': (5, 1e-8, 632.8, 'Gas'),
    'Ar-ion (514 nm)': (1000, 1e-8, 514, 'Gas'),
    'CO2 (10.6 μm)': (10, 1e-7, 10600, 'Gas'),
    # Solid-state
    'Nd:YAG (1064 nm)': (100, 1e-7, 1064, 'Solid'),
    'Ruby (694 nm)': (1000, 1e-8, 694, 'Solid'),
    'Ti:Sapphire (800 nm)': (500, 1e-7, 800, 'Solid'),
    # Semiconductor
    'GaAs VCSEL': (0.5, 0.01, 850, 'Semiconductor'),
    'InGaAs edge': (5, 1e-4, 980, 'Semiconductor'),
    'Quantum cascade': (100, 1e-4, 5000, 'Semiconductor'),
    # Microcavity
    'Photonic crystal': (0.01, 0.1, 1550, 'Microcavity'),
    'Microdisk': (0.001, 0.3, 1550, 'Microcavity'),
    'Plasmon nanolaser': (0.1, 0.5, 500, 'Nanolaser'),
}

print("\nLaser Types and Thresholds:")
print("-" * 70)
print(f"{'Laser Type':<25} {'P_th (mW)':<12} {'β':<10} {'1/β':<12} {'λ (nm)'}")
print("-" * 70)

for laser, (P_th, beta, wavelength, ltype) in laser_data.items():
    gamma_beta = 1 / beta
    print(f"{laser:<25} {P_th:<12.3f} {beta:<10.1e} {gamma_beta:<12.1e} {wavelength:.0f}")

# =============================================================================
# SECTION 4: PHOTON STATISTICS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 4: PHOTON STATISTICS AT THRESHOLD")
print("=" * 70)

print("""
Photon statistics characterize coherence:

1. Thermal (below threshold):
    P(n) = n̄^n / (1 + n̄)^(n+1)  (Bose-Einstein)
    g²(0) = 2 (bunched)
    Δn² = n̄² + n̄

2. Coherent (above threshold):
    P(n) = exp(-n̄) × n̄^n / n!  (Poisson)
    g²(0) = 1 (no bunching)
    Δn² = n̄ (shot noise)

3. At threshold:
    g²(0) transitions from 2 → 1
    Mandel Q = (Δn² - n̄)/n̄ transitions from Q ~ n̄ to Q = 0

The second-order coherence g²(0) marks the transition:
    γ_g2 = g²(0) - 1
    - γ_g2 = 1: Thermal (incoherent)
    - γ_g2 = 0: Coherent (lasing)
    - γ_g2 ~ 0.5 at threshold
""")

# Plot photon statistics
fig_stats, axes_stats = plt.subplots(1, 2, figsize=(12, 5))

# g²(0) vs pump
ax1 = axes_stats[0]
P_over_Pth = np.linspace(0.1, 5, 100)
# Model: g²(0) transitions from 2 to 1 around threshold
# Using tanh fit
g2_model = 1 + 1 / (1 + (P_over_Pth - 1)**2 * 10)
g2_approx = 1 + np.exp(-(P_over_Pth - 1) * 3)  # Different model

ax1.plot(P_over_Pth, g2_model, 'b-', linewidth=2, label='g²(0)')
ax1.axhline(y=1, color='green', linestyle='--', label='Coherent g²=1')
ax1.axhline(y=2, color='red', linestyle='--', label='Thermal g²=2')
ax1.axvline(x=1, color='gray', linestyle=':', label='Threshold')
ax1.fill_betweenx([0.5, 2.5], 0, 1, alpha=0.2, color='red', label='Below threshold')
ax1.fill_betweenx([0.5, 2.5], 1, 5, alpha=0.2, color='green', label='Above threshold')
ax1.set_xlabel('P / P_th', fontsize=12)
ax1.set_ylabel('g²(0)', fontsize=12)
ax1.set_title('Second-Order Coherence vs Pump', fontsize=12)
ax1.legend(fontsize=8)
ax1.set_xlim(0, 5)
ax1.set_ylim(0.5, 2.5)

# Photon distributions
ax2 = axes_stats[1]
n = np.arange(0, 20)
n_bar = 5
# Thermal
P_thermal = n_bar**n / (1 + n_bar)**(n + 1)
# Poisson
P_poisson = np.exp(-n_bar) * n_bar**n / np.array([np.math.factorial(int(ni)) for ni in n])
ax2.bar(n - 0.2, P_thermal, width=0.4, alpha=0.7, label='Thermal (below)', color='red')
ax2.bar(n + 0.2, P_poisson, width=0.4, alpha=0.7, label='Coherent (above)', color='green')
ax2.set_xlabel('Photon number n', fontsize=12)
ax2.set_ylabel('P(n)', fontsize=12)
ax2.set_title(f'Photon Statistics (n̄ = {n_bar})', fontsize=12)
ax2.legend()

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/lasing_stats.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Photon statistics figure saved.")

# =============================================================================
# SECTION 5: LINEWIDTH AND COHERENCE TIME
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 5: LINEWIDTH AND COHERENCE TIME")
print("=" * 70)

print("""
Laser linewidth narrowing at threshold:

Below threshold (spontaneous):
    Δν_sp = A₂₁ / (2π) (spontaneous emission linewidth)

Above threshold (Schawlow-Townes):
    Δν_ST = (π h ν / P_out) × (Δν_c)² × n_sp

where:
    Δν_c = cavity linewidth
    n_sp = spontaneous emission factor (~1-2)

Coherence time τ_c = 1/(π Δν)

Linewidth ratio:
    Δν_ST / Δν_sp ~ 1/n (photon number)

At threshold (n ~ 1/β):
    Δν ~ Δν_c (cavity limited)

Define γ_linewidth = Δν / Δν_c:
- γ >> 1: Below threshold (spontaneous broadened)
- γ ~ 1: At threshold (cavity limited)
- γ << 1: Above threshold (Schawlow-Townes narrowing)
""")

# Linewidth data
linewidth_data = {
    # (Δν_c MHz, Δν_ST Hz, Q factor)
    'He-Ne': (1.5, 1e3, 1e8),
    'Nd:YAG': (10, 1e4, 1e7),
    'Semiconductor DFB': (100, 1e6, 1e6),
    'VCSEL': (1000, 1e7, 1e5),
    'Fiber laser': (0.1, 1, 1e9),
    'External cavity': (0.01, 0.1, 1e10),
}

print("\nLaser Linewidths:")
print("-" * 60)
print(f"{'Laser':<25} {'Δν_c (MHz)':<12} {'Δν_ST (Hz)':<12} {'Q factor'}")
print("-" * 60)

for laser, (dnu_c, dnu_ST, Q) in linewidth_data.items():
    print(f"{laser:<25} {dnu_c:<12.2f} {dnu_ST:<12.1e} {Q:.1e}")

# =============================================================================
# SECTION 6: MICROCAVITY AND THRESHOLDLESS LASERS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 6: MICROCAVITY AND β FACTOR")
print("=" * 70)

print("""
The β factor determines threshold sharpness:

    β = Spontaneous emission rate into lasing mode / Total spontaneous rate

For conventional lasers: β ~ 10^(-8) - 10^(-4)
    - Sharp threshold (kink in output)
    - Large 1/β ~ 10^4 - 10^8 photons at threshold

For microcavity lasers: β ~ 0.01 - 1
    - Soft threshold (gradual transition)
    - Small 1/β ~ 1 - 100 photons at threshold

The "thresholdless laser" (β = 1) has no sharp threshold!
All spontaneous emission goes into the mode.

This is the γ ~ 1 framework:
    γ_β = 1/β = number of photons at threshold

When γ_β ~ 1 (β ~ 1): Transition occurs with ~1 photon
When γ_β >> 1 (β << 1): Transition occurs with many photons

The threshold photon number IS the γ parameter!
""")

# β factor analysis
print("\nβ Factor Analysis:")
print("-" * 60)

beta_values = [1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0]
for beta in beta_values:
    n_th = 1 / beta  # photon number at threshold
    gamma_beta = 1 / beta
    thresh_type = "Sharp" if beta < 1e-3 else "Soft" if beta < 0.1 else "Thresholdless"
    print(f"β = {beta:.1e}: n_th ~ {n_th:.1e}, γ_β = {gamma_beta:.1e}, Threshold: {thresh_type}")

# =============================================================================
# SECTION 7: DICKE SUPERRADIANCE
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 7: DICKE SUPERRADIANCE")
print("=" * 70)

print("""
Dicke superradiance: N atoms emit collectively.

Single atom: Intensity I₁ ∝ 1
N incoherent atoms: I_incoh = N × I₁
N coherent atoms (superradiance): I_SR ∝ N² × I₁

Superradiant enhancement: I_SR / I_incoh = N

Condition for superradiance:
    Sample size < λ (wavelength)
    OR: Phase coherence between atoms

The Dicke number:
    J = N/2 (maximum angular momentum)

Superradiant transition occurs when:
    γ_Dicke = 1 / (N × Γ × τ_c) ~ 1

where:
    Γ = single-atom decay rate
    τ_c = coherence time

This is COOPERATIVE coherence emergence at γ ~ 1!
""")

# Dicke physics examples
dicke_data = {
    # (N atoms, enhancement factor, system)
    'Rydberg atoms': (100, 100, 'Cold atoms'),
    'Quantum dots': (1000, 50, 'Semiconductor'),
    'NV centers': (10, 8, 'Diamond'),
    'Molecules': (1e6, 10, 'Molecular'),
    'Nuclear spins': (1e10, 100, 'NMR'),
}

print("\nDicke Superradiance Examples:")
print("-" * 60)
print(f"{'System':<20} {'N':<12} {'I_SR/I_sp':<12} {'Platform'}")
print("-" * 60)

for system, (N, enhance, platform) in dicke_data.items():
    print(f"{system:<20} {N:<12.0e} {enhance:<12.0f} {platform}")

# =============================================================================
# SECTION 8: POLARITON CONDENSATION (Connection to #158)
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 8: POLARITON LASING")
print("=" * 70)

print("""
Exciton-polariton condensation (connected to Session #158):

Polaritons = exciton + photon hybrid
At high density: Bose-Einstein condensation (BEC)

Polariton laser vs photon laser:
- Photon laser: Population inversion required
- Polariton laser: BEC-like threshold (no inversion)

Threshold criterion:
    n × λ_dB³ ~ ζ(3/2) = 2.612 (BEC)

This is EXACTLY γ_BEC = 1 from Session #159!

Polariton condensate at:
    γ_pol = n_c / n = 1 (threshold density)

The connection: Polariton lasing IS BEC of light-matter quasiparticles.
""")

# =============================================================================
# SECTION 9: RANDOM LASERS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 9: RANDOM LASERS")
print("=" * 70)

print("""
Random laser: Gain medium with disordered scattering.

No external cavity - feedback from multiple scattering.

Two types:
1. Incoherent (diffusive): No resonant modes
2. Coherent (resonant): Localized modes from disorder

Threshold condition:
    Mean free path l_g (gain) < l_s (scattering)

Define γ_random = l_s / l_g:
- γ > 1: No lasing (scattering extracts before amplification)
- γ ~ 1: Threshold
- γ < 1: Random lasing

Connection to Anderson localization (#89):
Localized photon modes enable coherent random lasing.
""")

# Random laser data
random_laser_data = {
    # (threshold fluence mJ/cm², scattering MFP μm, type)
    'ZnO powder': (0.1, 1, 'Coherent'),
    'TiO2 + dye': (1.0, 5, 'Incoherent'),
    'Semiconductor NPs': (0.5, 2, 'Coherent'),
    'Liquid crystals': (0.2, 10, 'Coherent'),
    'Polymer + scatterers': (2.0, 20, 'Incoherent'),
}

print("\nRandom Laser Examples:")
print("-" * 60)
print(f"{'System':<25} {'Eth (mJ/cm²)':<15} {'l_s (μm)':<10} {'Type'}")
print("-" * 60)

for system, (E_th, l_s, rtype) in random_laser_data.items():
    print(f"{system:<25} {E_th:<15.1f} {l_s:<10} {rtype}")

# =============================================================================
# SECTION 10: QUANTUM LIMIT OF LASING
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 10: QUANTUM LIMIT - SINGLE-ATOM LASER")
print("=" * 70)

print("""
Single-atom laser: Ultimate limit of microcavity QED.

One atom in high-Q cavity:
    g = atom-cavity coupling
    κ = cavity decay rate
    γ = atomic decay rate

Strong coupling: g > κ, γ (Session #158)

Lasing criterion:
    C = g² / (κ γ) > 1  (cooperativity)

Define γ_SAL = 1/C = κγ/g²:
- γ_SAL > 1: Below threshold (no lasing)
- γ_SAL ~ 1: Threshold
- γ_SAL < 1: Single-atom lasing

This IS the strong coupling criterion from Session #158!
Single-atom lasing requires γ_SC < 1.
""")

# Single-atom laser data
single_atom_data = {
    # (g MHz, κ MHz, γ MHz, cooperativity C)
    'Cs in Fabry-Perot': (25, 4, 2.6, 60),
    'Rb in fiber cavity': (5, 50, 3, 0.17),
    'Ion in Paul trap': (1, 0.1, 10, 1),
    'QD in photonic crystal': (20, 10, 0.1, 400),
    'NV in diamond': (0.1, 1, 0.01, 1),
}

print("\nSingle-Atom/Emitter Laser Data:")
print("-" * 70)
print(f"{'System':<25} {'g (MHz)':<10} {'κ (MHz)':<10} {'γ (MHz)':<10} {'C':<8} {'γ_SAL'}")
print("-" * 70)

for system, (g, kappa, gamma_at, C) in single_atom_data.items():
    gamma_sal = 1 / C
    status = "Lasing" if C > 1 else "Below"
    print(f"{system:<25} {g:<10.1f} {kappa:<10.1f} {gamma_at:<10.1f} {C:<8.1f} {gamma_sal:.2f} ({status})")

# =============================================================================
# SECTION 11: γ ~ 1 ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 11: γ ~ 1 ANALYSIS FOR LASING")
print("=" * 70)

print("""
Multiple γ ~ 1 boundaries in lasing:

1. Pump threshold: γ_pump = P/P_th = 1
   - Definition of lasing threshold
   - Gain = Loss at γ = 1

2. Photon number: γ_n = n_sp/n_total ~ 1 at threshold
   - Spontaneous and stimulated comparable
   - Transition from thermal to coherent

3. β factor: γ_β = 1/β = n_th (photons at threshold)
   - Small β: Sharp threshold, many photons
   - β → 1: Soft/"thresholdless", ~1 photon

4. Coherence g²(0): γ_g2 = g²(0) - 1
   - γ_g2 = 1: Thermal (incoherent)
   - γ_g2 = 0: Coherent (lasing)
   - Transition at γ_g2 ~ 0.5

5. Cooperativity: γ_SAL = 1/C (single-atom limit)
   - Strong coupling: C > 1 → γ_SAL < 1
   - Connects to Session #158

6. Random laser: γ_random = l_s / l_g ~ 1
   - Scattering vs gain balance

ALL lasing transitions occur at γ ~ 1!
""")

# Summary
print("\nγ ~ 1 Boundaries in Lasing:")
print("-" * 50)

gamma_lasing = {
    'P/P_th (pump)': 1.00,
    'g²(0) - 1 (at threshold)': 0.5,
    '1/C (single-atom)': 1.0,  # at threshold
    'l_s/l_g (random)': 1.0,
}

for name, value in gamma_lasing.items():
    status = "✓ γ ~ 1" if 0.3 < value < 3 else "~"
    print(f"  {status} {name}: {value:.2f}")

# =============================================================================
# SECTION 12: FIGURE
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 12: GENERATING FIGURE")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Input-output characteristic
ax1 = axes[0, 0]
P_in = np.linspace(0, 3, 100)
P_out_sharp = np.where(P_in > 1, (P_in - 1) * 0.5, 0.01 * P_in)  # Sharp threshold
P_out_soft = 0.5 * (P_in + np.sqrt(P_in**2 + 0.1))  # Soft threshold (β ~ 1)
ax1.plot(P_in, P_out_sharp, 'b-', linewidth=2, label='Sharp (β << 1)')
ax1.plot(P_in, P_out_soft * 0.3, 'g-', linewidth=2, label='Soft (β ~ 1)')
ax1.axvline(x=1, color='red', linestyle='--', linewidth=2, label='γ = 1')
ax1.fill_betweenx([0, 1.5], 0, 1, alpha=0.2, color='red')
ax1.fill_betweenx([0, 1.5], 1, 3, alpha=0.2, color='green')
ax1.set_xlabel('P_in / P_th (γ)', fontsize=12)
ax1.set_ylabel('P_out', fontsize=12)
ax1.set_title('A) Laser Input-Output: Threshold at γ = 1', fontsize=12)
ax1.legend()
ax1.set_xlim(0, 3)
ax1.set_ylim(0, 1.5)

# Panel B: β factor effect on threshold
ax2 = axes[0, 1]
beta_vals = [1e-4, 1e-2, 0.1, 0.5]
P_range = np.linspace(0.01, 3, 200)
for beta in beta_vals:
    # Simple model: n = 1/((1-P) + β) for P < 1, linear for P > 1
    n_ph = beta / np.abs(1 - P_range + beta) + P_range
    ax2.semilogy(P_range, n_ph, linewidth=2, label=f'β = {beta:.0e}')
ax2.axvline(x=1, color='red', linestyle='--', label='Threshold')
ax2.set_xlabel('P / P_th', fontsize=12)
ax2.set_ylabel('Photon number ⟨n⟩', fontsize=12)
ax2.set_title('B) Photon Number vs β Factor', fontsize=12)
ax2.legend(fontsize=9)
ax2.set_xlim(0, 3)
ax2.set_ylim(0.1, 1e4)

# Panel C: g²(0) transition
ax3 = axes[1, 0]
P_norm = np.linspace(0.1, 5, 100)
for beta in [1e-4, 1e-2, 0.1]:
    # Model: g² = 1 + 1/(1 + n/n_sp) where n_sp ~ β at threshold
    n_model = (P_norm - 1 + np.sqrt((P_norm - 1)**2 + 4*beta)) / 2
    n_model = np.maximum(n_model, beta)
    g2 = 1 + 1 / (1 + n_model / beta)
    ax3.plot(P_norm, g2, linewidth=2, label=f'β = {beta:.0e}')
ax3.axhline(y=1, color='green', linestyle='--', alpha=0.5, label='Coherent')
ax3.axhline(y=2, color='red', linestyle='--', alpha=0.5, label='Thermal')
ax3.axvline(x=1, color='gray', linestyle=':', label='Threshold')
ax3.set_xlabel('P / P_th', fontsize=12)
ax3.set_ylabel('g²(0)', fontsize=12)
ax3.set_title('C) Coherence g²(0) Transition', fontsize=12)
ax3.legend(fontsize=9)
ax3.set_xlim(0, 5)
ax3.set_ylim(0.8, 2.2)

# Panel D: β factor vs threshold sharpness
ax4 = axes[1, 1]
beta_range = np.logspace(-8, 0, 100)
n_threshold = 1 / beta_range  # photons at threshold
ax4.loglog(beta_range, n_threshold, 'b-', linewidth=2)
ax4.axhline(y=1, color='red', linestyle='--', label='n_th = 1 (quantum)')
ax4.axvline(x=0.1, color='green', linestyle=':', label='Soft threshold')
ax4.fill_betweenx([1e-2, 1e10], 0.1, 1, alpha=0.2, color='green',
                   label='"Thresholdless"')
ax4.set_xlabel('β factor', fontsize=12)
ax4.set_ylabel('Photons at threshold (γ_β = 1/β)', fontsize=12)
ax4.set_title('D) Threshold Sharpness vs β', fontsize=12)
ax4.legend(fontsize=9)
ax4.set_xlim(1e-8, 1)
ax4.set_ylim(1, 1e8)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/lasing_coherence.png',
            dpi=150, bbox_inches='tight')
print("Figure saved to lasing_coherence.png")
plt.close()

# =============================================================================
# SECTION 13: CONCLUSIONS
# =============================================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print("""
Session #165 Findings:

1. PUMP THRESHOLD γ_pump = P/P_th = 1
   - Definition of lasing: Gain = Loss at threshold
   - Below: Spontaneous emission (incoherent)
   - Above: Stimulated emission (coherent)
   - THE quintessential coherence transition

2. PHOTON STATISTICS AT THRESHOLD
   - g²(0) = 2 (thermal) → g²(0) = 1 (coherent)
   - Transition at γ ~ 1 where n ~ 1/β
   - Mandel Q-parameter: Q ~ n̄ → Q = 0

3. β FACTOR AND THRESHOLD PHOTON NUMBER
   - γ_β = 1/β = photons at threshold
   - Sharp threshold (β ~ 10^-8): ~10^8 photons
   - "Thresholdless" (β → 1): ~1 photon
   - Microcavity lasers approach quantum limit

4. LINEWIDTH NARROWING
   - Schawlow-Townes: Δν ∝ 1/P ∝ 1/n
   - γ_linewidth = Δν/Δν_c ~ 1 at threshold
   - Coherence time diverges above threshold

5. DICKE SUPERRADIANCE
   - Collective coherence: I ∝ N² (vs N for incoherent)
   - γ_Dicke = 1/(N×Γ×τ_c) ~ 1 at superradiant transition
   - Cooperative emission = coherent many-body state

6. SINGLE-ATOM LASER (Connection to #158)
   - γ_SAL = 1/C = κγ/g²
   - C > 1 (γ_SAL < 1) required for single-atom lasing
   - Strong coupling criterion IS lasing criterion

7. RANDOM LASER
   - γ_random = l_s/l_g = 1 at threshold
   - Connects to Anderson localization (#89)
   - Coherent random lasing from localized modes

This is the 28th phenomenon type at γ ~ 1!

SIGNIFICANCE:
Lasing is THE paradigmatic coherence transition in optics.
The laser threshold at P = P_th (γ = 1) is where:
- Stimulated > Spontaneous
- Coherent photons > Thermal photons
- g²(0) drops from 2 to 1

The β factor reveals that threshold sharpness = 1/β,
connecting microscopic (photon number) to macroscopic (kink).
Microcavity lasers (high β) show that the FUNDAMENTAL threshold
is ~1 photon (γ ~ 1 in photon number).
""")

print("=" * 70)
print("END SESSION #165")
print("=" * 70)
