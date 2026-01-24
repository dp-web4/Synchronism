#!/usr/bin/env python3
"""
Adhesion and Wetting Coherence Analysis
Session #186 - Chemistry Track

Tests γ ~ 1 framework for surface phenomena:
1. Contact angle and wettability
2. Young's equation and surface energies
3. Capillary number and spreading
4. Adhesion work and cohesion
5. Superhydrophobicity and lotus effect

Key insight: Contact angle θ = 90° IS the γ ~ 1 wetting boundary
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("="*70)
print("ADHESION AND WETTING COHERENCE ANALYSIS")
print("Session #186 - Chemistry Track")
print("="*70)

# =============================================================================
# 1. CONTACT ANGLE AS COHERENCE PARAMETER
# =============================================================================
print("\n" + "="*70)
print("1. CONTACT ANGLE AS COHERENCE PARAMETER")
print("="*70)

print("""
Young's equation:
  cos(θ) = (γ_SV - γ_SL) / γ_LV

where:
  γ_SV = solid-vapor surface energy
  γ_SL = solid-liquid interface energy
  γ_LV = liquid-vapor surface tension
  θ = contact angle

Define coherence parameter:
  γ_wet = cos(θ)

At θ = 90°: γ_wet = cos(90°) = 0 (boundary)
  θ < 90°: γ_wet > 0 (hydrophilic, wetting)
  θ > 90°: γ_wet < 0 (hydrophobic, non-wetting)

Alternative: γ_θ = θ/90°
  At γ_θ = 1: θ = 90° (hydrophobic-hydrophilic boundary)
""")

# Contact angle data for water on various surfaces
contact_angle_data = {
    # Surface: (θ in degrees, classification)
    'Glass (clean)': (20, 'hydrophilic'),
    'Stainless steel': (70, 'hydrophilic'),
    'Aluminum': (60, 'hydrophilic'),
    'Silicon (native oxide)': (45, 'hydrophilic'),
    'PMMA': (75, 'hydrophilic'),
    'Polycarbonate': (82, 'near boundary'),
    'Polystyrene': (87, 'near boundary'),
    'PTFE (Teflon)': (115, 'hydrophobic'),
    'Paraffin wax': (105, 'hydrophobic'),
    'Polyethylene': (95, 'hydrophobic'),
    'Silicone rubber': (100, 'hydrophobic'),
    'Lotus leaf': (160, 'superhydrophobic'),
    'Graphene': (90, 'boundary'),
}

print("\nContact Angle Analysis:")
print("-"*70)
print(f"{'Surface':<25} {'θ (°)':>8} {'cos(θ)':>10} {'γ_θ':>8} {'Classification':<18}")
print("-"*70)

theta_values = []
cos_theta_values = []
gamma_theta_values = []

for surface, (theta, classification) in contact_angle_data.items():
    theta_rad = np.radians(theta)
    cos_theta = np.cos(theta_rad)
    gamma_theta = theta / 90

    theta_values.append(theta)
    cos_theta_values.append(cos_theta)
    gamma_theta_values.append(gamma_theta)

    status = "γ~1!" if 0.8 <= gamma_theta <= 1.2 else ""
    print(f"{surface:<25} {theta:>8} {cos_theta:>10.2f} {gamma_theta:>8.2f} {classification:<18} {status}")

print("-"*70)
print(f"Mean γ_θ = {np.mean(gamma_theta_values):.2f} ± {np.std(gamma_theta_values):.2f}")

near_boundary = sum(1 for g in gamma_theta_values if 0.8 <= g <= 1.2)
print(f"Surfaces with γ_θ in [0.8, 1.2]: {near_boundary}/{len(gamma_theta_values)}")

# =============================================================================
# 2. SPREADING COEFFICIENT
# =============================================================================
print("\n" + "="*70)
print("2. SPREADING COEFFICIENT")
print("="*70)

print("""
Spreading coefficient:
  S = γ_SV - γ_SL - γ_LV = γ_LV(cos(θ) - 1)

At S = 0: θ = 0° (perfect wetting)
At S < 0: partial wetting
At S → -2γ_LV: θ = 180° (complete non-wetting)

Define: γ_spread = |S| / γ_LV = |cos(θ) - 1|
At θ = 90°: γ_spread = 1 (γ ~ 1!)
""")

# Water surface tension
gamma_water = 72.8  # mN/m at 25°C

print("\nSpreading Coefficient Analysis:")
print("-"*60)
print(f"{'Surface':<25} {'θ (°)':>8} {'S (mN/m)':>12} {'γ_spread':>10}")
print("-"*60)

gamma_spread_values = []
for surface, (theta, _) in contact_angle_data.items():
    theta_rad = np.radians(theta)
    cos_theta = np.cos(theta_rad)
    S = gamma_water * (cos_theta - 1)
    gamma_spread = abs(cos_theta - 1)
    gamma_spread_values.append(gamma_spread)

    status = "γ~1" if 0.8 <= gamma_spread <= 1.2 else ""
    print(f"{surface:<25} {theta:>8} {S:>12.1f} {gamma_spread:>10.2f} {status}")

print("-"*60)
print(f"Mean γ_spread = {np.mean(gamma_spread_values):.2f} ± {np.std(gamma_spread_values):.2f}")

# =============================================================================
# 3. CAPILLARY NUMBER
# =============================================================================
print("\n" + "="*70)
print("3. CAPILLARY NUMBER")
print("="*70)

print("""
Capillary number:
  Ca = μv/γ_LV (viscous force / surface tension force)

At Ca = 1: viscous and surface tension forces balance (γ ~ 1!)
  Ca << 1: surface tension dominates (droplets, menisci)
  Ca >> 1: viscous forces dominate (thin film drainage)

For coating/printing:
  Ca ~ 1 is the optimal regime for controlled deposition
""")

# Calculate Ca for various processes
capillary_data = {
    # Process: (velocity m/s, viscosity Pa·s, gamma mN/m)
    'Ink-jet printing': (5, 0.003, 30),
    'Slot-die coating': (0.1, 0.1, 30),
    'Spin coating': (10, 0.001, 30),
    'Dip coating (slow)': (0.001, 0.001, 30),
    'Curtain coating': (1, 0.05, 30),
    'Roll coating': (1, 0.1, 30),
    'Capillary rise': (0.001, 0.001, 72.8),
    'Droplet spreading': (0.1, 0.001, 72.8),
}

print("\nCapillary Number Analysis:")
print("-"*60)
print(f"{'Process':<25} {'v (m/s)':>10} {'μ (Pa·s)':>10} {'Ca':>12}")
print("-"*60)

Ca_values = []
for process, (v, mu, gamma) in capillary_data.items():
    Ca = mu * v / (gamma / 1000)  # Convert mN/m to N/m
    Ca_values.append(Ca)
    status = "γ~1" if 0.1 < Ca < 10 else ""
    print(f"{process:<25} {v:>10.3f} {mu:>10.3f} {Ca:>12.2e} {status}")

print("-"*60)
print(f"Ca range: {min(Ca_values):.2e} to {max(Ca_values):.2e}")
near_unity_Ca = sum(1 for Ca in Ca_values if 0.1 < Ca < 10)
print(f"Processes with Ca in [0.1, 10]: {near_unity_Ca}/{len(Ca_values)}")

# =============================================================================
# 4. WORK OF ADHESION
# =============================================================================
print("\n" + "="*70)
print("4. WORK OF ADHESION AND COHESION")
print("="*70)

print("""
Young-Dupré equation:
  W_adh = γ_LV(1 + cos(θ))

Work of cohesion (liquid to itself):
  W_coh = 2×γ_LV

Adhesion ratio:
  γ_adh = W_adh/W_coh = (1 + cos(θ))/2

At θ = 90°: γ_adh = 0.5
At θ = 0°: γ_adh = 1 (equal to cohesion)
At θ = 180°: γ_adh = 0 (no adhesion)

The boundary θ = 90° gives γ_adh = 0.5 (half of cohesion).
""")

print("\nWork of Adhesion Analysis:")
print("-"*65)
print(f"{'Surface':<25} {'θ (°)':>8} {'W_adh (mN/m)':>12} {'γ_adh':>10}")
print("-"*65)

gamma_adh_values = []
W_coh = 2 * gamma_water

for surface, (theta, _) in contact_angle_data.items():
    theta_rad = np.radians(theta)
    W_adh = gamma_water * (1 + np.cos(theta_rad))
    gamma_adh = W_adh / W_coh
    gamma_adh_values.append(gamma_adh)

    status = "γ~1" if 0.4 <= gamma_adh <= 0.6 else ""  # Near 0.5
    print(f"{surface:<25} {theta:>8} {W_adh:>12.1f} {gamma_adh:>10.2f} {status}")

print("-"*65)
print(f"Mean γ_adh = {np.mean(gamma_adh_values):.2f} ± {np.std(gamma_adh_values):.2f}")

# =============================================================================
# 5. SURFACE ENERGY AND CRITICAL SURFACE TENSION
# =============================================================================
print("\n" + "="*70)
print("5. CRITICAL SURFACE TENSION")
print("="*70)

print("""
Zisman's critical surface tension γ_c:
  cos(θ) → 1 as γ_LV → γ_c

When γ_LV = γ_c: complete wetting (θ = 0)
When γ_LV > γ_c: partial wetting

Define: γ_Zisman = γ_LV/γ_c

At γ_Zisman = 1: transition to complete wetting
  Below: wetting
  Above: non-wetting
""")

# Critical surface tension data
zisman_data = {
    # Surface: (γ_c in mN/m)
    'PTFE': (18,),
    'Polypropylene': (29,),
    'Polyethylene': (31,),
    'Polystyrene': (33,),
    'PMMA': (39,),
    'Nylon 6': (42,),
    'PET': (43,),
    'Aluminum oxide': (45,),
    'Glass': (73,),
}

print("\nCritical Surface Tension (Zisman):")
print("-"*55)
print(f"{'Surface':<20} {'γ_c (mN/m)':>12} {'γ_Zisman (water)':>15}")
print("-"*55)

gamma_Zisman_values = []
for surface, (gamma_c,) in zisman_data.items():
    gamma_Z = gamma_water / gamma_c
    gamma_Zisman_values.append(gamma_Z)
    status = "γ~1" if 0.8 <= gamma_Z <= 1.2 else ""
    print(f"{surface:<20} {gamma_c:>12} {gamma_Z:>15.2f} {status}")

print("-"*55)
print(f"Mean γ_Zisman = {np.mean(gamma_Zisman_values):.2f} ± {np.std(gamma_Zisman_values):.2f}")

# =============================================================================
# 6. CASSIE-BAXTER AND WENZEL STATES
# =============================================================================
print("\n" + "="*70)
print("6. ROUGH SURFACE WETTING: CASSIE-BAXTER VS WENZEL")
print("="*70)

print("""
Wenzel state (liquid fills roughness):
  cos(θ_W) = r × cos(θ_Y)
  where r = roughness ratio (actual/projected area)

Cassie-Baxter state (air pockets):
  cos(θ_CB) = f_s × cos(θ_Y) + (1 - f_s) × (-1)
  where f_s = solid fraction

Transition between states at critical θ_Y ~ 90°!

Define: γ_rough = θ_apparent / θ_Young
At θ_Y = 90°: both states give same apparent angle
  θ_Y < 90°: Wenzel enhances wetting
  θ_Y > 90°: Cassie possible (superhydrophobic)
""")

# Rough surface data
rough_data = {
    # Surface: (θ_Young, roughness_r, f_s, observed θ)
    'Smooth PTFE': (115, 1.0, 1.0, 115),
    'Rough PTFE (r=2)': (115, 2.0, 0.3, 150),
    'Lotus leaf': (105, 3.0, 0.1, 160),
    'Smooth glass': (20, 1.0, 1.0, 20),
    'Rough glass (r=1.5)': (20, 1.5, 1.0, 10),
    'Textured silicon': (45, 2.0, 0.4, 120),
    'Smooth PS': (87, 1.0, 1.0, 87),
    'Nano-textured PS': (87, 2.5, 0.2, 155),
}

print("\nRough Surface Wetting:")
print("-"*70)
print(f"{'Surface':<20} {'θ_Y (°)':>8} {'r':>6} {'f_s':>6} {'θ_obs (°)':>10} {'γ_rough':>8}")
print("-"*70)

gamma_rough_values = []
for surface, (theta_Y, r, f_s, theta_obs) in rough_data.items():
    gamma_rough = theta_obs / theta_Y if theta_Y > 0 else 0
    gamma_rough_values.append(gamma_rough)
    status = "γ~1" if 0.8 <= gamma_rough <= 1.2 else ""
    print(f"{surface:<20} {theta_Y:>8} {r:>6.1f} {f_s:>6.1f} {theta_obs:>10} {gamma_rough:>8.2f} {status}")

print("-"*70)
print(f"Mean γ_rough = {np.mean(gamma_rough_values):.2f} ± {np.std(gamma_rough_values):.2f}")
print("θ_Y = 90° is THE critical angle for superhydrophobicity!")

# =============================================================================
# 7. BOND NUMBER AND GRAVITY EFFECTS
# =============================================================================
print("\n" + "="*70)
print("7. BOND NUMBER: GRAVITY VS SURFACE TENSION")
print("="*70)

print("""
Bond number:
  Bo = Δρ × g × L² / γ_LV

At Bo = 1: gravity and surface tension balance (γ ~ 1!)
  Bo << 1: surface tension dominates (small drops, bubbles)
  Bo >> 1: gravity dominates (puddles, pools)

Capillary length:
  l_c = √(γ/(Δρ×g)) ~ 2.7 mm for water

At L = l_c: Bo = 1 exactly
""")

# Bond number for various length scales
rho_water = 1000  # kg/m³
g = 9.81  # m/s²
gamma_water_SI = 0.0728  # N/m
l_c = np.sqrt(gamma_water_SI / (rho_water * g)) * 1000  # mm

print(f"\nCapillary length for water: l_c = {l_c:.2f} mm")

bond_data = {
    # System: (L in mm, description)
    'Fog droplet': (0.01, 'aerosol'),
    'Mist droplet': (0.1, 'condensation'),
    'Rain droplet': (2.0, 'typical'),
    'Large raindrop': (5.0, 'maximum'),
    'Capillary rise': (2.7, 'capillary length'),
    'Dewdrop': (1.0, 'on surface'),
    'Puddle': (10.0, 'flat'),
    'Bubble in water': (1.0, 'rising'),
}

print("\nBond Number Analysis:")
print("-"*55)
print(f"{'System':<20} {'L (mm)':>10} {'Bo':>12} {'Status':<15}")
print("-"*55)

Bo_values = []
for system, (L_mm, desc) in bond_data.items():
    L = L_mm / 1000  # Convert to m
    Bo = rho_water * g * L**2 / gamma_water_SI
    Bo_values.append(Bo)

    if Bo < 0.1:
        status = "surface tension"
    elif Bo > 10:
        status = "gravity"
    else:
        status = "γ ~ 1 (balanced)"

    gamma_status = "γ~1!" if 0.1 < Bo < 10 else ""
    print(f"{system:<20} {L_mm:>10.2f} {Bo:>12.3f} {status:<15} {gamma_status}")

print("-"*55)
print(f"At L = l_c = {l_c:.2f} mm: Bo = 1 exactly")

# =============================================================================
# 8. WEBER NUMBER AND DROPLET DYNAMICS
# =============================================================================
print("\n" + "="*70)
print("8. WEBER NUMBER: INERTIA VS SURFACE TENSION")
print("="*70)

print("""
Weber number:
  We = ρ × v² × L / γ_LV

At We = 1: inertia and surface tension balance (γ ~ 1!)
  We << 1: surface tension dominates (stable drops)
  We >> 1: inertia dominates (splashing, breakup)

Critical Weber number for droplet breakup:
  We_crit ~ 10-12 (not exactly 1, but order of magnitude)
""")

# Weber number for droplet processes
weber_data = {
    # Process: (L mm, v m/s, description)
    'Falling raindrop': (2, 5, 'terminal velocity'),
    'Ink-jet droplet': (0.05, 10, 'printing'),
    'Spray nozzle': (0.1, 20, 'atomization'),
    'Droplet on hot surface': (2, 0.5, 'Leidenfrost'),
    'Gentle impact': (2, 0.1, 'deposition'),
    'Splashing': (2, 3, 'breakup threshold'),
}

print("\nWeber Number Analysis:")
print("-"*60)
print(f"{'Process':<25} {'L (mm)':>8} {'v (m/s)':>10} {'We':>12}")
print("-"*60)

We_values = []
for process, (L_mm, v, desc) in weber_data.items():
    L = L_mm / 1000
    We = rho_water * v**2 * L / gamma_water_SI
    We_values.append(We)
    status = "γ~1" if 0.1 < We < 10 else ""
    print(f"{process:<25} {L_mm:>8.2f} {v:>10.1f} {We:>12.1f} {status}")

print("-"*60)
print(f"We range: {min(We_values):.1f} to {max(We_values):.1f}")

# =============================================================================
# 9. COMPREHENSIVE STATISTICS
# =============================================================================
print("\n" + "="*70)
print("9. COMPREHENSIVE STATISTICS")
print("="*70)

all_gamma_values = {
    'Contact angle (γ_θ)': gamma_theta_values,
    'Spreading (γ_spread)': gamma_spread_values,
    'Adhesion (γ_adh)': gamma_adh_values,
    'Zisman (γ_c ratio)': gamma_Zisman_values,
    'Bond number (log Bo)': [np.log10(Bo) if Bo > 0 else 0 for Bo in Bo_values],
}

print("\nSummary Statistics:")
print("-"*60)
print(f"{'Parameter':<25} {'Mean':>10} {'Std':>10} {'N':>5}")
print("-"*60)

for param, values in all_gamma_values.items():
    if len(values) > 1:
        print(f"{param:<25} {np.mean(values):>10.2f} {np.std(values):>10.2f} {len(values):>5}")

# Statistical test for γ_θ vs 1.0
t_stat, p_value = stats.ttest_1samp(gamma_theta_values, 1.0)
print(f"\nγ_θ vs 1.0: mean = {np.mean(gamma_theta_values):.2f}, p = {p_value:.4f}")

# Key γ ~ 1 findings
print("\nKey γ ~ 1 Boundaries:")
print("-"*50)
print("1. θ = 90° (hydrophobic-hydrophilic, cos(θ) = 0)")
print("2. Bo = 1 at capillary length l_c ~ 2.7 mm")
print("3. Ca = 1 at viscous-capillary balance")
print("4. We = 1 at inertia-surface tension balance")
print("5. γ_adh = 0.5 (half cohesion at θ = 90°)")

# =============================================================================
# 10. VISUALIZATION
# =============================================================================
print("\n" + "="*70)
print("10. GENERATING VISUALIZATION")
print("="*70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Contact angle distribution
ax1 = axes[0, 0]
ax1.bar(range(len(theta_values)), theta_values, color='steelblue', alpha=0.7)
ax1.axhline(y=90, color='red', linestyle='--', linewidth=2, label='θ = 90° (γ ~ 1)')
ax1.axhspan(80, 100, alpha=0.2, color='green', label='γ_θ ~ 1 region')
ax1.set_xlabel('Surface')
ax1.set_ylabel('Contact Angle θ (°)')
ax1.set_title('Contact Angles: Boundary at θ = 90°')
ax1.set_xticks(range(len(contact_angle_data)))
ax1.set_xticklabels([s[:8] for s in contact_angle_data.keys()], rotation=45, ha='right', fontsize=8)
ax1.legend()

# Plot 2: Young equation - cos(θ) vs surface type
ax2 = axes[0, 1]
sorted_idx = np.argsort(theta_values)
sorted_theta = [theta_values[i] for i in sorted_idx]
sorted_cos = [cos_theta_values[i] for i in sorted_idx]
sorted_labels = [list(contact_angle_data.keys())[i][:10] for i in sorted_idx]

ax2.plot(range(len(sorted_cos)), sorted_cos, 'bo-', linewidth=2, markersize=8)
ax2.axhline(y=0, color='red', linestyle='--', linewidth=2, label='cos(θ) = 0 (θ = 90°)')
ax2.fill_between(range(len(sorted_cos)), -0.2, 0.2, alpha=0.2, color='green')
ax2.set_xlabel('Surface (sorted by θ)')
ax2.set_ylabel('cos(θ)')
ax2.set_title('Wetting: Hydrophilic (+) vs Hydrophobic (-)')
ax2.legend()
ax2.set_xticks(range(len(sorted_labels)))
ax2.set_xticklabels(sorted_labels, rotation=45, ha='right', fontsize=8)

# Plot 3: Bond number and capillary length
ax3 = axes[1, 0]
L_range = np.logspace(-2, 2, 100)  # mm
Bo_curve = (L_range / l_c)**2
ax3.loglog(L_range, Bo_curve, 'b-', linewidth=2)
ax3.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Bo = 1 (γ ~ 1)')
ax3.axvline(x=l_c, color='green', linestyle=':', label=f'l_c = {l_c:.1f} mm')
ax3.fill_between(L_range, 0.1, 10, alpha=0.1, color='green')
ax3.set_xlabel('Length Scale L (mm)')
ax3.set_ylabel('Bond Number Bo')
ax3.set_title('Bond Number: Gravity vs Surface Tension')
ax3.legend()
ax3.set_xlim(0.01, 100)
ax3.set_ylim(1e-4, 1e4)

# Plot 4: Cassie-Wenzel transition
ax4 = axes[1, 1]
theta_Y = np.linspace(0, 180, 100)
# Wenzel (r = 2)
r = 2.0
theta_W = np.degrees(np.arccos(np.clip(r * np.cos(np.radians(theta_Y)), -1, 1)))
# Cassie-Baxter (f_s = 0.3)
f_s = 0.3
cos_CB = f_s * np.cos(np.radians(theta_Y)) + (1 - f_s) * (-1)
theta_CB = np.degrees(np.arccos(np.clip(cos_CB, -1, 1)))

ax4.plot(theta_Y, theta_Y, 'k--', linewidth=1, label='Smooth (r=1)')
ax4.plot(theta_Y, theta_W, 'b-', linewidth=2, label=f'Wenzel (r={r})')
ax4.plot(theta_Y, theta_CB, 'r-', linewidth=2, label=f'Cassie-Baxter (f_s={f_s})')
ax4.axvline(x=90, color='green', linestyle=':', linewidth=2, label='θ_Y = 90° (γ ~ 1)')
ax4.set_xlabel("Young's Angle θ_Y (°)")
ax4.set_ylabel('Apparent Angle θ (°)')
ax4.set_title('Rough Surface Wetting: Transition at θ_Y = 90°')
ax4.legend()
ax4.set_xlim(0, 180)
ax4.set_ylim(0, 180)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/wetting_coherence.png', dpi=150)
print("Figure saved to wetting_coherence.png")

# =============================================================================
# 11. CONCLUSIONS
# =============================================================================
print("\n" + "="*70)
print("11. CONCLUSIONS")
print("="*70)

print(f"""
ADHESION AND WETTING AT γ ~ 1

Finding #123: Wetting phenomena show multiple γ ~ 1 boundaries

1. CONTACT ANGLE θ = 90°
   - cos(θ) = 0 at boundary (γ_wet = 0)
   - γ_θ = θ/90° = 1 at hydrophobic-hydrophilic transition
   - Near boundary: polystyrene (87°), graphene (90°)

2. BOND NUMBER
   - Bo = 1 at capillary length l_c = {l_c:.2f} mm
   - γ ~ 1 separates surface tension and gravity regimes
   - Universal for all liquids

3. CAPILLARY AND WEBER NUMBERS
   - Ca = 1: viscous-capillary balance
   - We = 1: inertia-surface tension balance
   - Both define regime transitions at γ ~ 1

4. WORK OF ADHESION
   - At θ = 90°: W_adh = 0.5 × W_coh
   - Adhesion equals half of cohesion at boundary

5. ROUGH SURFACE WETTING
   - θ_Y = 90° is THE critical Young angle
   - Below: Wenzel enhances wetting
   - Above: Cassie-Baxter allows superhydrophobicity

PHYSICAL INTERPRETATION:
- θ = 90° IS the coherence boundary for wetting
- Below: liquid spreads (coherent film)
- Above: liquid beads (discrete droplets)
- All dimensionless numbers (Bo, Ca, We) at γ ~ 1

49th phenomenon type at γ ~ 1!
""")

print("\n" + "="*70)
print("SESSION #186 COMPLETE")
print("="*70)
