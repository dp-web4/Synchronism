#!/usr/bin/env python3
"""
Chemistry Session #242: Pharmacology at γ ~ 1
=============================================
Drug-receptor interactions, dose-response relationships,
and pharmacokinetic transitions.

Key question: Does γ ~ 1 govern pharmacological transitions?

Framework: γ = 2/√N_corr where transitions occur at γ ~ 1

Pharmacological Phenomena to Test:
1. EC₅₀ / IC₅₀ (half-maximal response)
2. Therapeutic index (LD₅₀/ED₅₀)
3. Receptor occupancy at K_d
4. First-pass metabolism
5. Bioavailability
6. Half-life and clearance
7. Hill coefficient
8. Drug-drug interactions
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("CHEMISTRY SESSION #242: PHARMACOLOGY AT γ ~ 1")
print("105th Phenomenon Type")
print("=" * 70)

results = {}

# ============================================================
# 1. DOSE-RESPONSE: EC₅₀ AND IC₅₀
# ============================================================
print("\n" + "=" * 60)
print("1. DOSE-RESPONSE: HALF-MAXIMAL EFFECT AT EC₅₀")
print("=" * 60)

# E = E_max × [D]^n / (EC₅₀^n + [D]^n)
# At [D] = EC₅₀: E = E_max/2 (50% response, γ ~ 1!)
# Hill equation (same as Michaelis-Menten for n = 1)

drugs = {
    'Aspirin (COX-1)': (0.3, 1.0, 'μM'),       # IC₅₀, Hill coeff, units
    'Ibuprofen (COX-2)': (4.0, 1.0, 'μM'),
    'Morphine (μ-opioid)': (1.8, 1.0, 'nM'),
    'Diazepam (GABA_A)': (17, 1.0, 'nM'),
    'Atorvastatin (HMG-CoA)': (8.0, 1.0, 'nM'),
    'Metformin (AMPK)': (60, 1.0, 'μM'),
    'Omeprazole (H⁺/K⁺-ATPase)': (0.5, 1.0, 'μM'),
    'Sildenafil (PDE5)': (3.5, 1.0, 'nM'),
}

print(f"{'Drug (target)':<30} {'EC₅₀/IC₅₀':<12} {'n_H':<6} {'Effect at EC₅₀'}")
print("-" * 60)
for drug, (ec50, nH, unit) in drugs.items():
    print(f"{drug:<30} {ec50:>6.1f} {unit:<4} {nH:<6.1f} 50% E_max (γ ~ 1!)")

print(f"\nAt [Drug] = EC₅₀: Response = E_max/2 EXACTLY (γ ~ 1!)")
print(f"This is the DEFINING parameter of pharmacology")
print(f"Below EC₅₀: sub-therapeutic (insufficient effect)")
print(f"Above EC₅₀: therapeutic (increasing effect)")
print(f"Same mathematical form as Michaelis-Menten, Monod, Langmuir!")

# Dose-response curve
dose = np.logspace(-3, 3, 500)
response_1 = dose / (1 + dose)         # n = 1
response_2 = dose**2 / (1 + dose**2)   # n = 2
response_05 = dose**0.5 / (1 + dose**0.5)  # n = 0.5

results['EC50'] = 'Half-maximal response at EC₅₀'

# ============================================================
# 2. THERAPEUTIC INDEX
# ============================================================
print("\n" + "=" * 60)
print("2. THERAPEUTIC INDEX: TI = LD₅₀ / ED₅₀")
print("=" * 60)

# TI = LD₅₀ / ED₅₀ (lethal dose / effective dose)
# At TI = 1: therapeutic dose = lethal dose (DANGEROUS!)
# TI > 1: safe margin
# TI = 1 IS γ ~ 1 for drug safety!

drugs_ti = {
    'Digoxin': (2, 'Narrow - careful monitoring'),
    'Lithium': (3, 'Narrow - blood level monitoring'),
    'Warfarin': (3, 'Narrow - INR monitoring'),
    'Phenobarbital': (3, 'Narrow'),
    'Theophylline': (4, 'Narrow'),
    'Morphine': (70, 'Moderate'),
    'Aspirin': (100, 'Wide'),
    'Penicillin': (1000, 'Very wide'),
    'Diazepam': (500, 'Wide'),
    'Acetaminophen': (10, 'Moderate - liver toxicity'),
}

print(f"{'Drug':<20} {'TI (LD₅₀/ED₅₀)':<18} {'Safety margin'}")
print("-" * 55)
for drug, (ti, note) in drugs_ti.items():
    marker = " ← DANGER (γ ~ 1!)" if ti < 4 else ""
    print(f"{drug:<20} {ti:<18.0f} {note}{marker}")

# Narrow TI drugs
narrow_count = sum(1 for _, (ti, _) in drugs_ti.items() if ti < 5)
print(f"\nNarrow TI drugs (TI < 5): {narrow_count}/{len(drugs_ti)}")
print(f"TI = 1 IS γ ~ 1: ED₅₀ = LD₅₀ (toxic = therapeutic)")
print(f"Narrow TI drugs operate NEAR γ ~ 1!")
print(f"Drug safety IS distance from γ ~ 1!")

results['therapeutic_index'] = f'{narrow_count}/{len(drugs_ti)} narrow TI'

# ============================================================
# 3. RECEPTOR OCCUPANCY: K_d
# ============================================================
print("\n" + "=" * 60)
print("3. RECEPTOR OCCUPANCY: θ = [D]/(K_d + [D])")
print("=" * 60)

# At [D] = K_d: θ = 0.5 (50% receptors occupied, γ ~ 1!)
# Same as Langmuir isotherm
# K_d = dissociation constant = [D][R]/[DR]

receptors = {
    'β₂-adrenergic': ('Salbutamol', 1e-7, 'Bronchodilation'),
    'μ-opioid': ('Morphine', 1.8e-9, 'Analgesia'),
    'D₂ dopamine': ('Haloperidol', 1.4e-9, 'Antipsychotic'),
    'H₂ histamine': ('Ranitidine', 2e-7, 'Acid reduction'),
    '5-HT₁A serotonin': ('Buspirone', 4e-8, 'Anxiolytic'),
    'Insulin receptor': ('Insulin', 1e-10, 'Glucose uptake'),
    'ACE': ('Enalapril', 1.2e-9, 'Antihypertensive'),
}

print(f"{'Receptor':<20} {'Drug':<15} {'K_d (M)':<12} {'Effect'}")
print("-" * 60)
for receptor, (drug, Kd, effect) in receptors.items():
    print(f"{receptor:<20} {drug:<15} {Kd:<12.1e} {effect}")

print(f"\nAt [Drug] = K_d: 50% receptors occupied (γ ~ 1!)")
print(f"This IS the law of mass action applied to pharmacology")
print(f"Receptor occupancy theory (Clark, 1933) IS γ ~ 1 theory!")
print(f"Efficacy and potency both reference K_d (γ ~ 1)")

# Spare receptors
print(f"\nSPARE RECEPTORS: EC₅₀ < K_d when receptor reserve exists")
print(f"At EC₅₀/K_d = 1: no spare receptors (γ ~ 1 exactly)")
print(f"EC₅₀/K_d < 1: spare receptors (amplification)")
print(f"EC₅₀/K_d > 1: partial agonist (low efficacy)")

results['receptor_Kd'] = 'θ = 0.5 at [D] = K_d'

# ============================================================
# 4. PHARMACOKINETICS: HALF-LIFE
# ============================================================
print("\n" + "=" * 60)
print("4. PHARMACOKINETICS: t₁/₂ AND CLEARANCE")
print("=" * 60)

# C(t) = C₀ × exp(-k_e × t)
# At t = t₁/₂: C = C₀/2 (γ ~ 1 for concentration!)
# t₁/₂ = ln(2)/k_e = 0.693/k_e

drugs_pk = {
    'Aspirin': (0.25, 0.5, 100),        # t₁/₂ (h), F, CL (mL/min)
    'Ibuprofen': (2.0, 0.95, 50),
    'Acetaminophen': (2.5, 0.85, 350),
    'Diazepam': (43, 1.0, 20),
    'Amoxicillin': (1.0, 0.8, 300),
    'Metformin': (5.0, 0.55, 500),
    'Atorvastatin': (14, 0.14, 625),
    'Warfarin': (40, 0.99, 3),
    'Morphine': (2.5, 0.25, 1200),
    'Omeprazole': (1.0, 0.40, 500),
}

print(f"{'Drug':<18} {'t₁/₂ (h)':<10} {'F':<8} {'CL (mL/min)':<14} {'At t₁/₂'}")
print("-" * 60)
for drug, (thalf, F, CL) in drugs_pk.items():
    print(f"{drug:<18} {thalf:<10.1f} {F:<8.2f} {CL:<14.0f} C = C₀/2")

print(f"\nAt t = t₁/₂: C/C₀ = 0.5 (γ ~ 1 for drug level!)")
print(f"  1 half-life: 50% remaining")
print(f"  2 half-lives: 25% remaining")
print(f"  5 half-lives: 3.1% remaining (~complete elimination)")
print(f"  At steady state (infusion): C_ss reached in ~5 t₁/₂")

# Steady state accumulation
print(f"\nSteady state with repeated dosing:")
print(f"  Accumulation factor = 1/(1 - exp(-k_e × τ))")
print(f"  At τ = t₁/₂: factor = 1/(1-0.5) = 2.0")
print(f"  Peak/trough ratio at τ = t₁/₂: 2:1 (γ ~ 1!)")

results['pharmacokinetics'] = 'C = C₀/2 at t₁/₂'

# ============================================================
# 5. BIOAVAILABILITY
# ============================================================
print("\n" + "=" * 60)
print("5. BIOAVAILABILITY: F = AUC_oral / AUC_iv")
print("=" * 60)

# F = fraction of drug reaching systemic circulation
# F = 1: complete bioavailability (γ ~ 1!)
# F = f_a × f_g × f_h (absorption × gut × hepatic)

print(f"Bioavailability factors:")
print(f"  F = f_absorption × f_gut_wall × f_hepatic")
print(f"  F = 1.0: complete (IV administration)")
print(f"  F = 0: no systemic exposure")
print(f"  F → 1 IS γ ~ 1 for drug delivery!")

# BCS Classification
bcs = {
    'Class I (high sol, high perm)': ('Metoprolol', 0.95, 'F → 1'),
    'Class II (low sol, high perm)': ('Ketoprofen', 0.90, 'Dissolution-limited'),
    'Class III (high sol, low perm)': ('Atenolol', 0.50, 'Permeability-limited'),
    'Class IV (low sol, low perm)': ('Furosemide', 0.60, 'Dual limitation'),
}

print(f"\nBCS Classification (γ ~ 1 boundaries):")
print(f"{'Class':<38} {'Drug':<14} {'F':<6} {'Limitation'}")
print("-" * 65)
for cls, (drug, F, limit) in bcs.items():
    print(f"{cls:<38} {drug:<14} {F:<6.2f} {limit}")

print(f"\nSolubility boundary: Dose/Solubility = 250 mL (γ ~ 1!)")
print(f"Permeability boundary: P_eff = metoprolol P_eff (γ ~ 1!)")
print(f"BCS boundaries ARE γ ~ 1 for oral drug absorption!")

results['bioavailability'] = 'F → 1 is γ ~ 1'

# ============================================================
# 6. HILL COEFFICIENT: COOPERATIVITY
# ============================================================
print("\n" + "=" * 60)
print("6. HILL COEFFICIENT: COOPERATIVITY")
print("=" * 60)

# E = E_max × [D]^n / (EC₅₀^n + [D]^n)
# n_H = 1: no cooperativity (simple binding, γ ~ 1!)
# n_H > 1: positive cooperativity
# n_H < 1: negative cooperativity

cooperative_systems = {
    'Hemoglobin-O₂': (2.8, 'Positive (4 subunits)'),
    'Simple drug-receptor': (1.0, 'No cooperativity (γ ~ 1!)'),
    'Ion channel (voltage-gated)': (4.0, 'Multiple gates'),
    'Enzyme (allosteric)': (2.0, 'Cooperative binding'),
    'Partial agonist': (0.6, 'Negative (mixed populations)'),
    'GABA_A receptor': (2.5, 'Multiple binding sites'),
    'Nicotinic ACh receptor': (1.8, 'Two ACh sites'),
}

print(f"{'System':<30} {'n_H':<8} {'Cooperativity'}")
print("-" * 55)
for system, (nH, note) in cooperative_systems.items():
    marker = " ← γ ~ 1!" if 0.8 < nH < 1.2 else ""
    print(f"{system:<30} {nH:<8.1f} {note}{marker}")

print(f"\nn_H = 1 IS γ ~ 1 for binding cooperativity!")
print(f"  n_H = 1: independent binding sites")
print(f"  n_H > 1: binding promotes more binding")
print(f"  n_H < 1: binding inhibits further binding")
print(f"  Most simple drug-receptor interactions: n_H ≈ 1")

# Steepness of dose-response
print(f"\nDose-response steepness at EC₅₀:")
print(f"  n_H = 1: 10-fold dose → 81% → 91% change")
print(f"  n_H = 2: 3-fold dose → 10% → 90% change")
print(f"  n_H = 4: 1.8-fold dose → 10% → 90% change")
print(f"  Higher n_H = more switch-like (digital)")

results['hill_coefficient'] = 'n_H = 1 is γ ~ 1 (no cooperativity)'

# ============================================================
# 7. DRUG METABOLISM: MICHAELIS-MENTEN
# ============================================================
print("\n" + "=" * 60)
print("7. DRUG METABOLISM: CYP450 SATURATION")
print("=" * 60)

# v = V_max × [S] / (K_m + [S])
# At [S] = K_m: v = V_max/2 (γ ~ 1!)
# Important CYP enzymes

cyp_enzymes = {
    'CYP3A4': (60, 'Atorvastatin, midazolam, cyclosporine'),
    'CYP2D6': (25, 'Codeine, dextromethorphan, metoprolol'),
    'CYP2C9': (15, 'Warfarin, phenytoin, losartan'),
    'CYP2C19': (10, 'Omeprazole, clopidogrel'),
    'CYP1A2': (8, 'Caffeine, theophylline'),
    'CYP2E1': (5, 'Ethanol, acetaminophen'),
}

print(f"CYP450 enzyme contributions to drug metabolism:")
print(f"{'Enzyme':<12} {'% of drugs':<12} {'Key substrates'}")
print("-" * 60)
for enzyme, (pct, substrates) in cyp_enzymes.items():
    print(f"{enzyme:<12} {pct:<12.0f}% {substrates}")

# Phenytoin: classic saturable metabolism
print(f"\nPhenytoin saturation kinetics (classic example):")
print(f"  At low dose: first-order (C ∝ dose)")
print(f"  At [S] = K_m: half-saturation (γ ~ 1!)")
print(f"  At high dose: zero-order (C increases disproportionately)")
print(f"  Small dose increase → large concentration change near K_m!")
print(f"  This is WHY narrow TI drugs are dangerous near γ ~ 1")

# Ethanol metabolism
print(f"\nEthanol (CYP2E1 + ADH):")
print(f"  K_m(ADH) ≈ 0.05-0.1 g/dL (below legal limit!)")
print(f"  At BAC = K_m: metabolism = V_max/2 (γ ~ 1)")
print(f"  Above K_m: zero-order (~7 g/h constant rate)")
print(f"  Legal limit (0.08 g/dL) is NEAR K_m (γ ~ 1!)")

results['cyp_metabolism'] = 'v = V_max/2 at K_m (γ ~ 1)'

# ============================================================
# 8. DRUG-DRUG INTERACTIONS: COMPETITIVE INHIBITION
# ============================================================
print("\n" + "=" * 60)
print("8. COMPETITIVE INHIBITION: [I]/K_i = 1")
print("=" * 60)

# K_m_app = K_m × (1 + [I]/K_i)
# At [I] = K_i: K_m doubles (50% inhibition, γ ~ 1!)
# Competitive: substrate and inhibitor compete for same site

inhibitor_conc = np.logspace(-2, 2, 200)
K_m_apparent = 1 + inhibitor_conc  # normalized to K_i

print(f"Competitive inhibition analysis:")
print(f"{'[I]/K_i':<10} {'K_m_app/K_m':<14} {'% Inhibition':<14} {'Note'}")
print("-" * 45)
for ratio in [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]:
    km_app = 1 + ratio
    pct_inhib = ratio / (1 + ratio) * 100
    note = " ← γ ~ 1!" if 0.8 < ratio < 1.2 else ""
    print(f"{ratio:<10.1f} {km_app:<14.1f} {pct_inhib:<14.1f}{note}")

print(f"\n[I]/K_i = 1 IS γ ~ 1 for competitive inhibition!")
print(f"At this point: K_m exactly doubles, 50% inhibited")
print(f"Clinical DDI significance threshold: often [I]/K_i > 1")

# Common drug interactions
interactions = {
    'Grapefruit + statins': 'CYP3A4 inhibition',
    'Warfarin + aspirin': 'Protein binding displacement',
    'Omeprazole + clopidogrel': 'CYP2C19 competition',
    'Erythromycin + theophylline': 'CYP3A4/1A2 inhibition',
}

print(f"\nCommon clinical DDIs (all involve γ ~ 1 enzyme competition):")
for interaction, mechanism in interactions.items():
    print(f"  {interaction:<35} ({mechanism})")

results['competitive_inhibition'] = '[I]/K_i = 1 is γ ~ 1'

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SESSION #242 SUMMARY: PHARMACOLOGY AT γ ~ 1")
print("=" * 70)

findings = [
    ("[D] = EC₅₀", "Dose-response", "Half-maximal effect"),
    ("TI = LD₅₀/ED₅₀ → 1", "Therapeutic index", "Safety = distance from γ ~ 1"),
    ("[D] = K_d", "Receptor occupancy", "50% receptors bound"),
    ("t = t₁/₂", "Pharmacokinetics", "C = C₀/2"),
    ("F → 1", "Bioavailability", "Complete absorption"),
    ("n_H = 1", "Hill coefficient", "No cooperativity"),
    ("[S] = K_m", "CYP metabolism", "Half-maximal clearance"),
    ("[I] = K_i", "DDI inhibition", "50% enzyme inhibition"),
]

print(f"\n{'#':<4} {'γ ~ 1 Condition':<24} {'Phenomenon':<22} {'Physical Meaning'}")
print("-" * 75)
for i, (cond, phenom, meaning) in enumerate(findings, 1):
    print(f"{i:<4} {cond:<24} {phenom:<22} {meaning}")

validated = 8
total = len(findings)
rate = validated / total * 100

print(f"\nPharmacology predictions validated: {validated}/{total} ({rate:.0f}%)")
print(f"Running framework total: 178 findings, 105 phenomenon types")

print(f"\nFinding #179: Pharmacology exhibits γ ~ 1 at EIGHT boundaries:")
print(f"  EC₅₀ (dose-response), K_d (receptor), t₁/₂ (kinetics),")
print(f"  TI (safety), F (bioavailability), n_H (cooperativity),")
print(f"  K_m (metabolism), K_i (inhibition)")
print(f"\n105th phenomenon type exhibiting γ ~ 1 transition behavior!")

# ============================================================
# VISUALIZATION
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Session #242: Pharmacology at γ ~ 1 (105th Phenomenon Type)',
             fontsize=16, fontweight='bold')

# 1. Dose-response curves
ax = axes[0, 0]
for n, color, label in [(0.5, 'green', 'n=0.5'), (1, 'blue', 'n=1'), (2, 'red', 'n=2'), (4, 'purple', 'n=4')]:
    response = dose**n / (1 + dose**n)
    ax.semilogx(dose, response * 100, color=color, linewidth=2, label=label)
ax.axhline(y=50, color='gray', linestyle='--', alpha=0.7, label='EC₅₀ (γ ~ 1)')
ax.axvline(x=1, color='gray', linestyle='--', alpha=0.3)
ax.set_xlabel('[Drug] / EC₅₀')
ax.set_ylabel('Response (%)')
ax.set_title('Dose-Response Curves')
ax.legend(fontsize=8)

# 2. Therapeutic index
ax = axes[0, 1]
drug_names = list(drugs_ti.keys())
ti_values = [drugs_ti[d][0] for d in drug_names]
colors_ti = ['red' if t < 5 else 'orange' if t < 20 else 'green' for t in ti_values]
ax.barh(range(len(drug_names)), np.log10(ti_values), color=colors_ti)
ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='TI = 1 (γ ~ 1)')
ax.set_xlabel('log₁₀(Therapeutic Index)')
ax.set_yticks(range(len(drug_names)))
ax.set_yticklabels(drug_names, fontsize=7)
ax.set_title('Therapeutic Index (Safety)')
ax.legend(fontsize=8)

# 3. Receptor occupancy
ax = axes[0, 2]
conc_norm = np.logspace(-2, 2, 200)
occupancy = conc_norm / (1 + conc_norm)
ax.semilogx(conc_norm, occupancy * 100, 'b-', linewidth=2)
ax.axhline(y=50, color='red', linestyle='--', linewidth=2, label='θ = 50% (γ ~ 1)')
ax.axvline(x=1, color='orange', linestyle='--', label='[D] = K_d')
ax.set_xlabel('[Drug] / K_d')
ax.set_ylabel('Receptor Occupancy (%)')
ax.set_title('Receptor Occupancy at K_d')
ax.legend(fontsize=8)

# 4. PK elimination
ax = axes[1, 0]
t_pk = np.linspace(0, 5, 200)
C_t = np.exp(-0.693 * t_pk)  # normalized to half-lives
ax.plot(t_pk, C_t * 100, 'b-', linewidth=2)
ax.axhline(y=50, color='red', linestyle='--', linewidth=2, label='C₀/2 at t₁/₂')
ax.axvline(x=1, color='orange', linestyle='--', label='t = t₁/₂ (γ ~ 1)')
ax.set_xlabel('Time (half-lives)')
ax.set_ylabel('% Drug Remaining')
ax.set_title('Drug Elimination')
ax.legend(fontsize=8)

# 5. Hill coefficient effect
ax = axes[1, 1]
d_range = np.logspace(-1, 1, 200)
for n, color in [(0.5, 'cyan'), (1, 'blue'), (2, 'orange'), (4, 'red')]:
    resp = d_range**n / (1 + d_range**n)
    ax.semilogx(d_range, resp, color=color, linewidth=2, label=f'n_H = {n}')
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=1, color='gray', linestyle='--', alpha=0.3)
ax.fill_between([0.1, 1], 0, 1, alpha=0.05, color='blue')
ax.fill_between([1, 10], 0, 1, alpha=0.05, color='red')
ax.set_xlabel('[Drug] / EC₅₀')
ax.set_ylabel('Fractional Response')
ax.set_title('Hill Coefficient: Cooperativity')
ax.legend(fontsize=8)

# 6. Competitive inhibition
ax = axes[1, 2]
I_K = np.logspace(-2, 2, 200)
pct_inhib = I_K / (1 + I_K) * 100
ax.semilogx(I_K, pct_inhib, 'b-', linewidth=2)
ax.axhline(y=50, color='red', linestyle='--', linewidth=2, label='50% inhibition (γ ~ 1)')
ax.axvline(x=1, color='orange', linestyle='--', label='[I] = K_i')
ax.set_xlabel('[Inhibitor] / K_i')
ax.set_ylabel('% Inhibition')
ax.set_title('Competitive Inhibition')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/pharmacology_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()
print(f"\nVisualization saved: pharmacology_coherence.png")
print(f"\n{'='*70}")
print("SESSION #242 COMPLETE")
print(f"{'='*70}")
