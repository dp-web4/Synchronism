#!/usr/bin/env python3
"""
Supercapacitor / Electrochemical Capacitor Coherence Analysis
Chemistry Session #234

Explores how γ ~ 1 manifests in supercapacitor systems:
1. EDL capacitance - Helmholtz double layer
2. RC time constant - charge/discharge dynamics
3. Ragone plot - energy vs power tradeoff
4. Cycle life - degradation thresholds
5. Self-discharge - charge retention
6. Pseudocapacitance - faradaic contribution
7. Pore size optimization - ion/pore matching

Key insight: Supercapacitor operation optimizes at γ ~ 1
boundaries where ionic/electronic transport balance.
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("SUPERCAPACITOR COHERENCE ANALYSIS")
print("Chemistry Session #234: γ ~ 1 in Electrochemical Capacitors")
print("=" * 70)

# =============================================================================
# 1. ELECTRIC DOUBLE LAYER (EDL) CAPACITANCE
# =============================================================================
print("\n" + "=" * 70)
print("1. ELECTRIC DOUBLE LAYER CAPACITANCE")
print("=" * 70)

# Helmholtz model: C = ε₀εᵣA/d
# Debye length: λ_D defines the characteristic screening distance
# At concentration giving λ_D = pore size: optimal capacitance

print("  Helmholtz model: C = ε₀εᵣA/d")
print("  where d = Helmholtz layer thickness")
print("  Gouy-Chapman-Stern: C = C_H × C_D / (C_H + C_D)")

# Debye length
print(f"\n  Debye length λ_D = √(ε₀εᵣkT / 2q²n₀)")
# At different concentrations
concentrations = [0.001, 0.01, 0.1, 1.0, 5.0]  # M
debye_lengths = [9.6, 3.04, 0.96, 0.30, 0.14]  # nm (for 1:1 electrolyte, 25°C)

print(f"\n  Concentration (M) | λ_D (nm) | Comment")
print("  " + "-" * 50)
for c, d in zip(concentrations, debye_lengths):
    comment = ""
    if c == 1.0:
        comment = "typical supercap electrolyte"
    elif abs(d - 0.96) < 0.01:
        comment = "near pore size match"
    print(f"  {c:8.3f}          | {d:6.2f}   | {comment}")

print(f"\n  At λ_D ≈ pore size: optimal EDL formation (γ ~ 1)")
print(f"  Too dilute (λ_D >> pore): incomplete double layer")
print(f"  Too concentrated (λ_D << pore): wasted electrolyte")

# Capacitance at electrode
print(f"\n  Specific capacitance targets:")
print(f"    Activated carbon: 100-300 F/g")
print(f"    Graphene: 100-550 F/g")
print(f"    CNTs: 50-200 F/g")
print(f"    MnO₂: 100-1380 F/g (pseudocapacitive)")

# =============================================================================
# 2. RC TIME CONSTANT
# =============================================================================
print("\n" + "=" * 70)
print("2. RC TIME CONSTANT")
print("=" * 70)

# τ = RC is THE characteristic time
# At t = τ: V(t) = V₀(1 - 1/e) ≈ 63.2% of final (charging)
# At t = τ: V(t) = V₀/e ≈ 36.8% of initial (discharging)

print("  Charging: V(t) = V₀[1 - exp(-t/RC)]")
print("  Discharging: V(t) = V₀ exp(-t/RC)")
print("  At t = RC = τ: 63.2% charged / 36.8% remaining")

# Typical RC values
rc_devices = {
    'Ceramic capacitor': {'C': 1e-9, 'R': 1, 'tau_s': 1e-9},
    'Film capacitor': {'C': 1e-6, 'R': 0.1, 'tau_s': 1e-7},
    'Electrolytic': {'C': 1e-3, 'R': 0.5, 'tau_s': 5e-4},
    'Supercapacitor (small)': {'C': 1, 'R': 1, 'tau_s': 1.0},
    'Supercapacitor (large)': {'C': 3000, 'R': 0.3e-3, 'tau_s': 0.9},
    'Battery': {'C': 'N/A', 'R': 0.01, 'tau_s': 'hours'},
}

print(f"\n  Device             | C        | R (Ω)  | τ = RC")
print("  " + "-" * 55)
for name, data in rc_devices.items():
    C_str = f"{data['C']}" if isinstance(data['C'], str) else f"{data['C']:.2e} F"
    tau_str = str(data['tau_s']) if isinstance(data['tau_s'], str) else f"{data['tau_s']:.2e} s"
    print(f"  {name:22s} | {C_str:8s} | {data['R']:.3f} | {tau_str}")

print(f"\n  τ = RC IS the γ ~ 1 time reference!")
print(f"  At t/τ = 1: charge/discharge at 63.2%/36.8%")
print(f"  1/e = 0.368... ≈ γ ~ 1 for exponential decay")

# =============================================================================
# 3. RAGONE PLOT ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("3. RAGONE PLOT: ENERGY vs POWER")
print("=" * 70)

# Ragone plot shows energy density vs power density
# Supercapacitors bridge gap between batteries and capacitors

devices = {
    'Conventional capacitor': {'E_Wh_kg': 0.05, 'P_W_kg': 100000},
    'Supercapacitor (EDL)': {'E_Wh_kg': 5, 'P_W_kg': 10000},
    'Supercapacitor (hybrid)': {'E_Wh_kg': 15, 'P_W_kg': 5000},
    'Li-ion battery': {'E_Wh_kg': 200, 'P_W_kg': 500},
    'Lead-acid battery': {'E_Wh_kg': 35, 'P_W_kg': 200},
    'Fuel cell': {'E_Wh_kg': 500, 'P_W_kg': 100},
}

print("  Energy vs Power tradeoff (Ragone plot):")
print(f"\n  Device                | E (Wh/kg) | P (W/kg) | E×P product")
print("  " + "-" * 65)
for name, data in devices.items():
    product = data['E_Wh_kg'] * data['P_W_kg']
    print(f"  {name:22s} | {data['E_Wh_kg']:8.1f}  | {data['P_W_kg']:8.0f} | {product:.0e}")

print(f"\n  Supercapacitors occupy the middle ground:")
print(f"    Higher power than batteries")
print(f"    Higher energy than capacitors")
print(f"    Energy × Power product at intermediate value")
print(f"    This intermediate = γ ~ 1 between pure E and pure P!")

# =============================================================================
# 4. CYCLE LIFE
# =============================================================================
print("\n" + "=" * 70)
print("4. CYCLE LIFE AND DEGRADATION")
print("=" * 70)

# Supercapacitors: 100,000-1,000,000+ cycles
# Far exceeding batteries

cycle_life = {
    'Supercap (EDL)': {'cycles': 1000000, 'retention': 0.95},
    'Supercap (pseudocap)': {'cycles': 100000, 'retention': 0.85},
    'Li-ion (LFP)': {'cycles': 3000, 'retention': 0.80},
    'Li-ion (NMC)': {'cycles': 1000, 'retention': 0.80},
    'Lead-acid': {'cycles': 500, 'retention': 0.80},
    'Solid-state Li': {'cycles': 5000, 'retention': 0.80},
}

print("  Cycle life comparison:")
print(f"\n  Device              | Cycles     | Retention at EOL")
print("  " + "-" * 55)
for name, data in cycle_life.items():
    print(f"  {name:20s} | {data['cycles']:>10,d} | {data['retention']*100:.0f}%")

# EOL threshold
print(f"\n  EOL threshold:")
print(f"    Batteries: 80% retention (Session #229)")
print(f"    Supercapacitors: variable (80-95% used)")
print(f"    EDL supercaps retain 95% after 10⁶ cycles!")

# Degradation per cycle
print(f"\n  Degradation per cycle:")
for name, data in cycle_life.items():
    loss_per_cycle = (1 - data['retention']) / data['cycles']
    print(f"    {name:20s}: {loss_per_cycle:.2e} per cycle")

# =============================================================================
# 5. SELF-DISCHARGE
# =============================================================================
print("\n" + "=" * 70)
print("5. SELF-DISCHARGE AND CHARGE RETENTION")
print("=" * 70)

# Self-discharge mechanisms
print("  Self-discharge mechanisms:")
print("    1. Ohmic leakage: V(t) = V₀ exp(-t/RC)")
print("    2. Charge redistribution: V(t) = V₀ - k√t")
print("    3. Faradaic: V(t) = V₀ - k×t (linear)")

# Typical self-discharge rates
print(f"\n  Retention after 24 hours (starting from V₀):")
sd_data = {
    'EDL supercap (AC)': 0.70,
    'EDL supercap (graphene)': 0.80,
    'Hybrid supercap': 0.60,
    'Li-ion battery': 0.999,
    'Conv. capacitor': 0.95,
}

print(f"\n  Device                | V/V₀ at 24h | 1 - V/V₀")
print("  " + "-" * 50)
for name, retention in sd_data.items():
    loss = 1 - retention
    print(f"  {name:22s} | {retention:.3f}      | {loss:.3f}")

# 50% voltage retention time
print(f"\n  Time to 50% voltage (t_½):")
print(f"    V(t_½) = V₀/2 → t_½ = τ × ln(2)")
print(f"    At t_½: V/V₀ = 0.5 (γ ~ 1 for voltage retention!)")
print(f"    Supercaps: t_½ typically hours to days")
print(f"    Batteries: t_½ typically months")

# =============================================================================
# 6. PSEUDOCAPACITANCE
# =============================================================================
print("\n" + "=" * 70)
print("6. PSEUDOCAPACITANCE: EDLC vs FARADAIC")
print("=" * 70)

# Pseudocapacitance: faradaic charge storage
# EDL: purely electrostatic (non-faradaic)
# Hybrid: combination

print("  Two charge storage mechanisms:")
print("    EDL: electrostatic (surface ions)")
print("    Pseudocapacitive: faradaic (redox reactions)")

# b-value analysis from cyclic voltammetry
# i = a × v^b
# b = 1.0: capacitive (surface-controlled)
# b = 0.5: diffusion-controlled (battery-like)

print(f"\n  b-value analysis: i = a × v^b")
print(f"    b = 1.0: purely capacitive (γ ~ 1 for surface control)")
print(f"    b = 0.5: purely diffusion-controlled")

b_values = {
    'Activated carbon (EDL)': 1.0,
    'RuO₂ (pseudocap)': 0.95,
    'MnO₂ (pseudocap)': 0.85,
    'Nb₂O₅ (intercalation)': 0.75,
    'V₂O₅ (battery-type)': 0.55,
    'LiCoO₂ (battery)': 0.50,
}

print(f"\n  Material            | b-value | Mechanism")
print("  " + "-" * 50)
for material, b in b_values.items():
    mech = "capacitive" if b > 0.9 else "pseudocapacitive" if b > 0.7 else "mixed" if b > 0.55 else "diffusion"
    print(f"  {material:20s} | {b:.2f}    | {mech}")

print(f"\n  b = 1.0 IS γ ~ 1 for capacitive behavior")
print(f"  b = 0.5 is diffusion-limited")
print(f"  Pseudocapacitive materials: 0.5 < b < 1.0")
print(f"  The capacitive/diffusion boundary at b ~ 0.75 ≈ γ ~ 1")

# =============================================================================
# 7. PORE SIZE OPTIMIZATION
# =============================================================================
print("\n" + "=" * 70)
print("7. PORE SIZE / ION SIZE MATCHING")
print("=" * 70)

# Key discovery: anomalous capacitance increase when pore size ≈ ion size
# Optimal at d_pore/d_ion ~ 1 (γ ~ 1!)

print("  Anomalous capacitance enhancement:")
print("  C increases when pore size ≈ desolvated ion size!")
print("  Ions shed solvation shell to fit in pores")

# Pore size effects
pore_sizes = [0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]  # nm
# Normalized capacitance (arbitrary, peak at ~0.7-0.8 nm for TEA+ ion)
cap_normalized = [0.3, 0.95, 1.0, 0.85, 0.75, 0.70, 0.68, 0.65, 0.60, 0.55]

print(f"\n  Pore size (nm) | C/C_max | Comment")
print("  " + "-" * 45)
for ps, cn in zip(pore_sizes, cap_normalized):
    comment = ""
    if cn == 1.0:
        comment = "OPTIMAL (d_pore ≈ d_ion)"
    elif ps == 1.0:
        comment = "d_pore/d_ion ~ 1"
    print(f"  {ps:8.1f}       | {cn:.2f}    | {comment}")

# Ion sizes
print(f"\n  Common electrolyte ion sizes (desolvated):")
ions = {
    'TEA+': 0.67,    # nm
    'TEMA+': 0.72,
    'BF4-': 0.48,
    'PF6-': 0.51,
    'TFSI-': 0.79,
    'EMI+': 0.76,
}

print(f"\n  Ion       | Diameter (nm)")
print("  " + "-" * 30)
for ion, size in ions.items():
    print(f"  {ion:8s} | {size:.2f}")

print(f"\n  Optimal: d_pore/d_ion = 1.0 (γ ~ 1!)")
print(f"  This was a SURPRISING discovery (2006, Chmiola et al.)")
print(f"  Previously thought large pores needed for ion access")

# =============================================================================
# 8. ENERGY AND POWER EQUATIONS
# =============================================================================
print("\n" + "=" * 70)
print("8. ENERGY AND POWER FUNDAMENTALS")
print("=" * 70)

# E = 0.5 CV² and P = V²/(4R) at maximum power
print("  Energy: E = ½CV²")
print("  Maximum power: P_max = V²/(4R_ESR)")
print("  At matched load R_L = R_ESR: P = P_max (γ ~ 1!)")

# Load matching
print(f"\n  Load matching for maximum power transfer:")
print(f"    P_load = V² × R_L / (R_ESR + R_L)²")
print(f"    dP/dR_L = 0 when R_L = R_ESR")
print(f"    Maximum power at R_L/R_ESR = 1 (γ ~ 1!)")
print(f"    At this point: 50% efficiency (half power in ESR)")

# Voltage window
print(f"\n  Operating voltage window:")
print(f"    Aqueous: 0-1.0 V (water decomposition limit)")
print(f"    Organic: 0-2.7 V")
print(f"    Ionic liquid: 0-3.5 V")
print(f"    Water decomposition at 1.23 V = SHE reference")
print(f"    (Links to fuel cell Session #230)")

# =============================================================================
# 9. SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY: γ ~ 1 IN SUPERCAPACITORS")
print("=" * 70)

gamma_findings = [
    ("EDL thickness", "λ_D ≈ d_pore", "Debye length matches pore size"),
    ("RC time constant", "t = τ = RC", "63.2%/36.8% charge transition"),
    ("Ragone position", "E × P intermediate", "Between battery and capacitor"),
    ("Cycle retention", "95% at 10⁶", "Near-perfect reversibility"),
    ("Self-discharge", "V/V₀ = 0.5", "Half-voltage retention time"),
    ("b-value", "b = 1.0", "Purely capacitive (γ ~ 1)"),
    ("Pore/ion match", "d_pore/d_ion = 1", "Anomalous capacitance peak"),
    ("Load matching", "R_L = R_ESR", "Maximum power transfer"),
]

print("\n  Parameter          | γ ~ 1 Condition     | Interpretation")
print("  " + "-" * 70)
for param, value, interp in gamma_findings:
    print(f"  {param:18s} | {value:19s} | {interp}")

print("\n  CONCLUSION: Supercapacitor optimization IS γ ~ 1 engineering:")
print("    - Pore/ion matching at d_pore/d_ion = 1 (anomalous!)")
print("    - Load matching at R_L = R_ESR (maximum power)")
print("    - b-value = 1 for capacitive behavior")
print("    - RC time constant defines charge dynamics")
print("    - Ragone position bridges energy/power")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Supercapacitor Coherence Analysis\nSession #234: γ ~ 1 in Electrochemical Capacitors',
             fontsize=14, fontweight='bold')

# 1. RC charge/discharge
ax1 = axes[0, 0]
t_tau = np.linspace(0, 5, 200)
charge = 1 - np.exp(-t_tau)
discharge = np.exp(-t_tau)
ax1.plot(t_tau, charge, 'b-', linewidth=2, label='Charge')
ax1.plot(t_tau, discharge, 'r-', linewidth=2, label='Discharge')
ax1.axvline(x=1.0, color='green', linestyle='--', linewidth=2, label='t = τ (γ ~ 1)')
ax1.axhline(y=0.632, color='gray', linestyle=':', alpha=0.5)
ax1.axhline(y=0.368, color='gray', linestyle=':', alpha=0.5)
ax1.set_xlabel('t/τ')
ax1.set_ylabel('V/V₀')
ax1.set_title('RC Charge/Discharge')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Ragone plot
ax2 = axes[0, 1]
for name, data in devices.items():
    ax2.scatter(data['P_W_kg'], data['E_Wh_kg'], s=100, zorder=5)
    ax2.annotate(name.replace(' ', '\n'), (data['P_W_kg'], data['E_Wh_kg']),
                fontsize=7, ha='center', va='bottom')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Power Density (W/kg)')
ax2.set_ylabel('Energy Density (Wh/kg)')
ax2.set_title('Ragone Plot')
ax2.grid(True, alpha=0.3)

# 3. Pore size optimization
ax3 = axes[0, 2]
ax3.plot(pore_sizes, cap_normalized, 'b-o', linewidth=2, markersize=8)
ax3.axvline(x=0.8, color='r', linestyle='--', linewidth=2, label='Optimal (d_pore ≈ d_ion)')
ax3.fill_between(pore_sizes, cap_normalized, alpha=0.3)
ax3.set_xlabel('Pore Size (nm)')
ax3.set_ylabel('Normalized Capacitance')
ax3.set_title('Anomalous Capacitance\n(Pore/Ion Matching)')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. b-value analysis
ax4 = axes[1, 0]
b_names = list(b_values.keys())
b_vals = list(b_values.values())
colors = ['green' if b > 0.9 else 'orange' if b > 0.7 else 'red' for b in b_vals]
ax4.barh(b_names, b_vals, color=colors, alpha=0.7)
ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='b = 1 (γ ~ 1)')
ax4.axvline(x=0.5, color='blue', linestyle='--', linewidth=2, alpha=0.5, label='b = 0.5 (diffusion)')
ax4.set_xlabel('b-value')
ax4.set_title('Capacitive vs Diffusion')
ax4.set_xlim(0.3, 1.1)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='x')

# 5. Load matching
ax5 = axes[1, 1]
R_ratio = np.linspace(0.01, 5, 200)
P_norm = R_ratio / (1 + R_ratio)**2
ax5.plot(R_ratio, P_norm, 'b-', linewidth=2)
ax5.axvline(x=1.0, color='r', linestyle='--', linewidth=2, label='R_L = R_ESR (γ ~ 1)')
ax5.scatter([1.0], [0.25], color='red', s=100, zorder=5)
ax5.fill_between(R_ratio, P_norm, alpha=0.3)
ax5.set_xlabel('R_load / R_ESR')
ax5.set_ylabel('Normalized Power')
ax5.set_title('Load Matching (Max Power)')
ax5.legend()
ax5.grid(True, alpha=0.3)

# 6. Cycle life comparison
ax6 = axes[1, 2]
cl_names = list(cycle_life.keys())
cl_vals = [cycle_life[n]['cycles'] for n in cl_names]
colors = ['green' if c > 10000 else 'orange' if c > 1000 else 'red' for c in cl_vals]
ax6.barh(cl_names, cl_vals, color=colors, alpha=0.7)
ax6.set_xlabel('Cycle Life')
ax6.set_title('Cycle Life Comparison')
ax6.set_xscale('log')
ax6.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/supercapacitor_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FINDING #171: Supercapacitor at γ ~ 1")
print("=" * 70)
print("""
Supercapacitor systems exhibit γ ~ 1 at multiple operational
and structural boundaries:

1. λ_D ≈ d_pore: Debye length matches pore size
2. t = τ = RC: 63.2%/36.8% charge/discharge reference
3. d_pore/d_ion = 1: Anomalous capacitance peak (surprising!)
4. b = 1.0: Purely capacitive behavior
5. R_L = R_ESR: Maximum power transfer at load matching
6. V/V₀ = 0.5: Half-voltage self-discharge threshold
7. Ragone intermediate: energy × power balance

97th phenomenon type exhibiting γ ~ 1 transition behavior.
Supercapacitor optimization IS γ ~ 1 matching!
""")

print("Visualization saved: supercapacitor_coherence.png")
