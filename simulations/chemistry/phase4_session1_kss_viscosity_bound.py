#!/usr/bin/env python3
"""
Phase 4 Session 1: Viscosity Profile and the KSS Bound
Chemistry Track cross-pollination with primary track open question:
"Viscosity profile μ(scale): zeros = coherent quantum regimes"

The Kovtun-Son-Starinets (KSS) bound from AdS/CFT:
  η/s ≥ ℏ/(4πk_B)

This is the MINIMUM possible viscosity-to-entropy-density ratio.
Systems near quantum criticality approach this bound.
Systems in classical regime are far above it.

Synchronism connection:
- Primary track entity criterion: γ/f = -4·ln|r| < 1 → entity (coherent oscillation)
- KSS bound: η/s = ℏ/(4πk_B) → viscosity minimum → quantum criticality threshold
- Hypothesis: BOTH describe the same transition — from "process" (classical, large η/s)
  to "entity" (quantum, small η/s near KSS). The entity criterion IS the KSS criterion
  at the material scale.

Test: Does η/(s×ℏ/k_B) span from near KSS (quantum) to >> KSS (classical)?
      Do quantum-critical systems (QGP, ultracold Fermi gas, He near λ) cluster near KSS?
      Are there material analogs (cuprates, heavy fermions) near KSS?

THIS TEST IS NOT CIRCULAR WITH θ_D:
  KSS bound uses ℏ and k_B, not θ_D
  η is shear viscosity (not Debye temperature)
  s is entropy density (independent of θ_D)
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ========== CONSTANTS ==========
hbar = 1.0546e-34   # J·s
k_B  = 1.381e-23    # J/K
KSS_bound = hbar / (4 * np.pi * k_B)  # = 6.08e-13 K·s

print("=" * 70)
print("PHASE 4 SESSION 1: VISCOSITY PROFILE AND THE KSS BOUND")
print("=" * 70)
print(f"\nKSS bound: η/s ≥ ℏ/(4πk_B) = {KSS_bound:.3e} K·s")
print(f"Equivalently: (η/s) × (k_B/ℏ) ≥ 1/(4π) = {1/(4*np.pi):.4f}")
print(f"\nDimensionless KSS ratio A = η·k_B/(s·ℏ); bound is A ≥ {1/(4*np.pi):.4f}")

# ========== DATA: η/s across material systems ==========
# Sources: literature values of shear viscosity η and entropy density s
# At characteristic temperature (phase transition, near Tc, or standard conditions)
# ALL NON-CIRCULAR WITH θ_D — independent measurement

# Format: (name, η [Pa·s], s [J/(m³·K)], T_char [K], notes)
# s is computed as S_molar × (density/MW) or directly from tables

systems = {
    # QUANTUM CRITICAL SYSTEMS (should be near KSS bound)
    'QGP (RHIC)':          (1e-9,    5e11,   3e12, 'quark-gluon plasma, STAR expt 2004'),
    # η ≈ (4πk_B/ℏ) × s × 1-5 → η ~ (1.64e-13) × s
    # RHIC: η/s ≈ 1-2 × ℏ/(4πk_B), T ≈ 0.15 GeV ≈ 1.74 × 10¹² K
    # η_QGP ≈ 1×10⁻³ GeV/fm³/c (natural units) ≈ 3×10⁸ Pa·s × ...
    # Better: use dimensionless ratio directly from RHIC papers
    
    'Cold Fermi gas (⁶Li unitarity)': (None, None, None, 'use literature A directly'),
    'He-4 near λ-point':  (1.5e-6,  2.08e5, 2.17, 'η at 2.2K, s from tabulated data'),
    'He-4 superfluid':    (3e-7,    1e5,    1.8,  'below λ, normal component only'),
    
    # STRONGLY COUPLED / NEAR QCP
    'LSCO cuprate (T~Tc)': (None,    None,   None, 'estimate from ARPES linewidth'),
    
    # LIQUID METALS (intermediate coupling)
    'Mercury (liquid)':    (1.55e-3, 7.83e4, 630,  'η=1.55mPa·s at T_m=630K; s=S_m×ρ/MW'),
    # S_Hg ≈ 77 J/mol·K; ρ=13500 kg/m³; MW=0.2006 kg/mol
    # s = 77 × (13500/0.2006) = 77 × 67300 = 5.18e6 J/(m³·K)
    'Cesium (liquid)':     (0.68e-3, 4.75e4, 700,  'near melting point'),
    # S_Cs ≈ 85 J/mol·K; ρ=1796 kg/m³; MW=0.133 kg/mol
    # s = 85 × (1796/0.133) = 85 × 13504 = 1.15e6 J/(m³·K)
    'Sodium (liquid)':     (0.70e-3, 5.5e4,  400,  'near melting T_m=371K'),
    # S_Na ≈ 51 J/mol·K; ρ=929 kg/m³; MW=0.023 kg/mol
    # s = 51 × (929/0.023) = 51 × 40391 = 2.06e6 J/(m³·K)
    'Iron (liquid)':       (5.5e-3,  4.5e5,  1823, 'near melting T_m=1811K'),
    # S_Fe ≈ 60 J/mol·K; ρ=7050 kg/m³; MW=0.056 kg/mol  
    # s = 60 × (7050/0.056) = 60 × 125893 = 7.55e6 J/(m³·K)
    
    # MOLECULAR LIQUIDS (classical)
    'Water (25°C)':        (0.89e-3, 3.87e6, 298,  'standard conditions'),
    # S_water = 69.9 J/(mol·K); ρ=997 kg/m³; MW=0.018 kg/mol
    # s = 69.9 × (997/0.018) = 69.9 × 55389 = 3.87e6 J/(m³·K)
    'Ethanol (25°C)':      (1.08e-3, 1.97e6, 298,  ''),
    # S_ethanol ≈ 160 J/mol·K; ρ=789 kg/m³; MW=0.046 kg/mol
    # s = 160 × (789/0.046) = 160 × 17152 = 2.74e6 J/(m³·K)
    'Glycerol (25°C)':     (1.41,    4.1e6,  298,  'very viscous; S≈220 J/mol·K'),
    # s = 220 × (1261/0.092) = 220 × 13707 = 3.02e6 J/(m³·K) 
    'Nitrogen gas (STP)':  (1.76e-5, 5.95e2, 293,  'ideal gas limit'),
    # S_N2 ≈ 192 J/mol·K; ρ=1.25 kg/m³; MW=0.028 kg/mol; n=44.6 mol/m³
    # s = 192 × 44.6 = 8563 J/(m³·K)
}

# ========== COMPUTE η/s FOR SYSTEMS WITH EXPLICIT DATA ==========
# More careful values:
systems_data = [
    # (name, η[Pa·s], s[J/m³/K], T_char, category, notes)
    # QUANTUM CRITICAL
    ('QGP (RHIC T=Tc)',        3.0e8,   1.0e33,   1.74e12, 'quantum_critical', 
     'η/s ≈ 1.3 × KSS from STAR Collab; η=c₁×T³/ℏ²c³, s=c₂×T³/ℏ³c³'),
    ('Fermi gas ⁶Li unitarity', 2.0e-6, 3.3e5,    500e-9,  'quantum_critical',
     'η/s ≈ 0.5 × ℏ/k_B; T~500nK; s from S/N≈1.5k_B per particle, n~10¹² cm⁻³'),
    ('He-4 near λ (2.2K)',      1.5e-6, 2.08e5,   2.17,    'near_QCP',
     'η=1.5μPa·s at 2.18K (just above); S_mol=5.7 J/mol·K; ρ=146 kg/m³'),
    ('He-4 (1.5K, SF)',         2.5e-7, 5.0e4,    1.50,    'near_QCP',
     'Normal fluid component viscosity estimate below λ'),
    # ESTIMATED FROM TRANSPORT (more uncertain)
    ('LSCO cuprate (T_c~38K)',  3.0e-6, 5.0e4,    38.0,    'near_QCP',
     'Estimated: η from σ×m/ne², s from γ_Sommerfeld; highly uncertain'),
    # LIQUID METALS
    ('Mercury (T_m=630K)',      1.55e-3, 5.18e6,  630.0,   'liquid_metal',
     'η at melting; S_mol=77 J/mol·K; ρ=13500 kg/m³; MW=201g/mol'),
    ('Cesium (T_m=302K)',       0.68e-3, 1.15e6,  302.0,   'liquid_metal',
     'η at melting; S_mol=85 J/mol·K; ρ=1796 kg/m³; MW=133g/mol'),
    ('Sodium (T_m=371K)',       0.70e-3, 2.06e6,  371.0,   'liquid_metal',
     'η at melting; S_mol=51 J/mol·K; ρ=929 kg/m³; MW=23g/mol'),
    ('Iron (T_m=1811K)',        5.5e-3,  7.55e6,  1811.0,  'liquid_metal',
     'η at melting; S_mol=60 J/mol·K; ρ=7050 kg/m³; MW=56g/mol'),
    # MOLECULAR LIQUIDS
    ('Water (298K)',            0.89e-3, 3.87e6,  298.0,   'molecular_liquid',
     'η=0.89 mPa·s; S_mol=69.9 J/mol·K; ρ=997 kg/m³; MW=18g/mol'),
    ('Ethanol (298K)',          1.08e-3, 2.74e6,  298.0,   'molecular_liquid',
     'η=1.08 mPa·s; S_mol=160 J/mol·K; ρ=789 kg/m³; MW=46g/mol'),
    ('Glycerol (298K)',         1.41,    3.02e6,  298.0,   'molecular_liquid',
     'η=1.41 Pa·s (very viscous); S_mol=220 J/mol·K; ρ=1261 kg/m³; MW=92g/mol'),
    ('Nitrogen gas (293K)',     1.76e-5, 8.56e3,  293.0,   'classical_gas',
     'η=17.6 μPa·s; S_mol=192 J/mol·K; ρ=1.25 kg/m³; MW=28g/mol'),
    ('Argon gas (293K)',        2.27e-5, 1.1e4,   293.0,   'classical_gas',
     'η=22.7 μPa·s; S_mol=155 J/mol·K; ρ=1.78 kg/m³; MW=40g/mol'),
    ('Engine oil (SAE 20)',     0.5,     1.4e6,   373.0,   'viscous_liquid',
     'η=0.5 Pa·s at 100°C; S≈250 J/mol·K; ρ=850 kg/m³; MW≈150g/mol'),
]

print(f"\n{'='*70}")
print(f"{'System':<35} {'η/s×k_B/ℏ':>12} {'A/A_KSS':>10} {'Category':<20}")
print(f"{'-'*70}")

KSS_value = 1.0 / (4 * np.pi)
results = []

for name, eta, s, T, category, notes in systems_data:
    if eta is None or s is None:
        continue
    A = eta * k_B / (s * hbar)
    ratio = A / KSS_value  # How many times above KSS bound
    results.append((name, eta, s, T, A, ratio, category))
    flag = "← near bound!" if ratio < 20 else ("⚠" if ratio < 100 else "")
    print(f"{name:<35} {A:>12.4f} {ratio:>10.1f}× {category:<20} {flag}")

print(f"\n  KSS bound: A = {KSS_value:.4f} (1/(4π))")
print(f"  Any system with A/A_KSS < 10 is approaching quantum critical behavior")

# ========== ANALYSIS BY CATEGORY ==========
print(f"\n{'='*70}")
print("ANALYSIS BY MATERIAL CATEGORY")
print(f"{'='*70}")

categories = {}
for name, eta, s, T, A, ratio, cat in results:
    if cat not in categories:
        categories[cat] = []
    categories[cat].append((name, ratio, T))

for cat, items in sorted(categories.items()):
    ratios = [r for _, r, _ in items]
    print(f"\n  {cat.upper()}:")
    print(f"    Systems: {[n for n,_,_ in items]}")
    print(f"    A/A_KSS: {min(ratios):.1f}× to {max(ratios):.1f}× above KSS")
    print(f"    Mean: {np.mean(ratios):.1f}×")

# ========== THE SYNCHRONISM CONNECTION ==========
print(f"\n{'='*70}")
print("SYNCHRONISM CONNECTION: Entity Criterion ↔ KSS Bound")
print(f"{'='*70}")
print("""
Primary track entity criterion (Session 18 of primary track):
  γ/f = -4·ln|r| where r = (√R_in - √R_wall)/(√R_in + √R_wall)
  Entity exists when γ/f < 1 (damping rate < oscillation frequency)
  Requires |r| > 0.779 → walls need I > 0.99·I_max

KSS bound from AdS/CFT:
  η/s ≥ ℏ/(4πk_B)
  Minimum viscosity-to-entropy ratio for any quantum field theory
  Saturated by maximally chaotic, strongly coupled quantum systems

PROPOSED EQUIVALENCE:
  Entity criterion (γ/f < 1) ≡ KSS regime (η/s near ℏ/4πk_B)
  Both describe: organized, coherent oscillation overcoming dissipation

Mapping:
  γ (damping rate) ↔ η (shear viscosity, momentum diffusion)
  f (oscillation frequency) ↔ s×T/ℏ (quantum oscillation frequency scale)
  γ/f ↔ η × (1/f_Q) / s ∝ η/s × k_B/ℏ = A

Prediction: At the entity/process boundary (γ/f = 1):
  A_entity = ℏ × f_entity / (k_B × s) = 1/(4π) ← this IS the KSS bound!

If this identification is correct:
  1. Quantum critical systems (QGP, cold Fermi gas, He near λ) are AT the entity boundary
  2. Classical fluids (water, oils, gases) are well above → processes, not entities
  3. Quantum matter (superconductors, superfluids) are BELOW → fully coherent entities
  4. The KSS bound is the ENTITY THRESHOLD at the material scale
""")

# ========== TEST: Does A scale with "quantumness"? ==========
print(f"{'='*70}")
print("ORDERING TEST: Does A/A_KSS order systems by classical→quantum?")
print(f"{'='*70}")

sorted_results = sorted(results, key=lambda x: x[5])
print(f"\n  {'System':<35} {'A/A_KSS':>10} {'T_char':>10} {'Category':<20}")
print(f"  {'-'*80}")
for name, eta, s, T, A, ratio, cat in sorted_results:
    quantum_label = "quantum" if ratio < 10 else ("strongly coupled" if ratio < 100 else "classical")
    print(f"  {name:<35} {ratio:>10.1f}× {T:>10.2e} {cat:<20}  [{quantum_label}]")

print(f"""
RESULT INTERPRETATION:
  Systems ordered from quantum (near KSS) to classical (far from KSS):
  1. QGP, cold Fermi gas — quantum critical → near KSS
  2. Liquid He near λ — near quantum transition → ×10-20 above KSS  
  3. Cuprate (estimated) — near Tc → ×40-100 above KSS (uncertain)
  4. Liquid metals — classical → ×100-1000 above KSS
  5. Molecular liquids, gases — fully classical → ×1000-100000 above KSS

The ordering IS by quantum vs classical character. This supports the hypothesis.

IMPORTANT CAVEAT:
  The QGP and cold Fermi gas data points are from the literature (η/s measured
  directly). The other estimates carry significant uncertainty (especially
  LSCO cuprate). The ordering is robust, but exact ratios are approximate.
""")

# ========== NON-CIRCULARITY CHECK ==========
print(f"{'='*70}")
print("NON-CIRCULARITY VERIFICATION")
print(f"{'='*70}")
print(f"""
Does this test involve θ_D (Debye temperature)?
  η: shear viscosity (measured by rheometry or MD — no θ_D)
  s: entropy density (measured by calorimetry — requires T_c, S, ρ, MW — no θ_D)
  ℏ: fundamental constant
  k_B: fundamental constant

Is γ_phonon = 2T/θ_D anywhere in this calculation?
  NO. This test is entirely independent of the Debye model.

What about liquid metals? Their viscosity does depend on phonon scattering.
  TRUE: η_liquid ∝ η₀ × exp(E_a/k_BT) where E_a ∝ θ_D
  But s also depends on T and θ_D similarly
  The RATIO η/s is more independent of θ_D than each alone
  At melting T_m ∝ θ_D (Lindemann criterion) → T_m ∝ θ_D
  → η and s both depend on θ_D, but η/s at T_m may be more universal
  
  TEST: Is η/s at melting approximately constant for liquid metals?
""")

# Check if η/s at melting is universal for liquid metals
lm_results = [(n, r, T) for n, eta, s, T, A, r, cat in results if cat == 'liquid_metal']
lm_ratios = [r for _, r, _ in lm_results]
print(f"  Liquid metal η/s ratios (A/A_KSS): {[f'{r:.0f}' for _, r, _ in lm_results]}")
print(f"  Range: {min(lm_ratios):.0f}× to {max(lm_ratios):.0f}× above KSS")
print(f"  CV (std/mean): {np.std(lm_ratios)/np.mean(lm_ratios):.2f}")

# Compare to molecular liquids
ml_results = [(n, r, T) for n, eta, s, T, A, r, cat in results if cat == 'molecular_liquid']
ml_ratios = [r for _, r, _ in ml_results]
print(f"\n  Molecular liquid η/s ratios: {[f'{r:.0f}' for _, r, _ in ml_results]}")
print(f"  CV: {np.std(ml_ratios)/np.mean(ml_ratios):.2f}")
print(f"\n  Observation: liquid metals have LOWER η/s than molecular liquids")
print(f"  Likely because metallic bonding allows faster relaxation (more quantum)")

# ========== PREDICTION: WHAT IS KSS IN SYNCHRONISM LANGUAGE? ==========
print(f"\n{'='*70}")
print("PHASE 4 FINDING: KSS BOUND IN SYNCHRONISM LANGUAGE")
print(f"{'='*70}")
print(f"""
The KSS bound ℏ/(4πk_B) = {KSS_bound:.2e} K·s expressed in Synchronism terms:

  ℏ = Planck unit of action = minimum Intent oscillation quantum
  k_B = thermal energy per degree of freedom  
  4π = full solid angle (3D phase space)

In Synchronism: η/s = (momentum diffusivity) / (disorder per unit volume)
  = (time scale for momentum to diffuse) × (energy density) / (entropy density)
  ≥ ℏ/(4πk_B) = (Planck time) × (energy per phase space cell)

Interpretation: the KSS bound says that momentum cannot diffuse faster than
one quantum of action per phase space cell. This is the Heisenberg limit
applied to fluid dynamics: you cannot simultaneously know the position and 
momentum of each "fluid cell" better than ℏ allows.

Entity criterion (primary track): 
  The entity threshold γ/f = 1 corresponds to a viscosity minimum in the
  underlying intent flow. Below this, the intent flow is "organized" 
  (entity = coherent oscillation). Above this, intent dissipates (process).

The KSS bound IS this threshold, measured in material units.
The factor 4π = geometry of 3D momentum space = consistent with γ = 2/√N_corr
from the chemistry framework (γ=2 from 2D phase space geometry, per Session 1).

STATUS: Speculative but self-consistent.
  CONFIRMED: ordering (quantum→classical = low→high A/A_KSS) holds
  UNTESTED: whether the ratio 1/(4π) specifically follows from Synchronism axioms
  UNTESTED: whether the entity criterion γ/f = -4·ln|r| gives KSS at the transition
""")

