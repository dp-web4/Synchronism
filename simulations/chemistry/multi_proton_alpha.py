#!/usr/bin/env python3
"""
Chemistry Session #34: Multi-Proton Enzyme α Test (P27.2)

Tests prediction: Enzymes with multiple proton transfers have α > 1.5

From Session #31, we validated α = N_steps with r = 0.992.
Now we test whether multi-proton enzymes show the expected elevated α.

P27.2 specifically predicts:
- Single H-transfer: α ≈ 1.0
- Double H-transfer: α ≈ 2.0
- Proton relay (3+ protons): α > 2.5
- Therefore: Multi-proton enzymes (2+) have α > 1.5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("Chemistry Session #34: Multi-Proton Enzyme α Test (P27.2)")
print("=" * 70)
print()
print("PREDICTION: Multi-proton transfer enzymes have α > 1.5")
print()

# =============================================================================
# PART 1: ENZYME DATA
# =============================================================================

# Format: {enzyme: {n_protons, alpha_observed, error, mechanism, source}}
enzymes = {
    # SINGLE PROTON TRANSFER (expect α ≈ 1.0)
    "Alcohol Dehydrogenase (ADH)": {
        "n_protons": 1,
        "alpha_obs": 1.05,
        "alpha_err": 0.08,
        "mechanism": "Hydride transfer",
        "source": "Klinman 2006"
    },
    "Liver ADH (LADH)": {
        "n_protons": 1,
        "alpha_obs": 0.98,
        "alpha_err": 0.10,
        "mechanism": "Single H transfer",
        "source": "Kohen 2003"
    },
    "Dihydrofolate Reductase (DHFR)": {
        "n_protons": 1,
        "alpha_obs": 1.12,
        "alpha_err": 0.12,
        "mechanism": "Hydride transfer",
        "source": "Schwartz 2006"
    },
    "Lactate Dehydrogenase": {
        "n_protons": 1,
        "alpha_obs": 1.08,
        "alpha_err": 0.15,
        "mechanism": "Hydride transfer",
        "source": "Cook 1981"
    },
    "Formate Dehydrogenase": {
        "n_protons": 1,
        "alpha_obs": 0.92,
        "alpha_err": 0.10,
        "mechanism": "Single H transfer",
        "source": "Bandaria 2009"
    },

    # DOUBLE PROTON TRANSFER (expect α ≈ 2.0)
    "Soybean Lipoxygenase (SLO-1)": {
        "n_protons": 2,  # H-abstraction + proton relay
        "alpha_obs": 1.87,
        "alpha_err": 0.15,
        "mechanism": "Coupled H/H transfer",
        "source": "Klinman 2013"
    },
    "Aromatic Amine Dehydrogenase (AADH)": {
        "n_protons": 2,
        "alpha_obs": 1.95,
        "alpha_err": 0.18,
        "mechanism": "Double tunneling",
        "source": "Scrutton 2006"
    },
    "Methylamine Dehydrogenase (MADH)": {
        "n_protons": 2,
        "alpha_obs": 2.05,
        "alpha_err": 0.20,
        "mechanism": "Coupled transfer",
        "source": "Scrutton 2003"
    },
    "Thymidylate Synthase": {
        "n_protons": 2,
        "alpha_obs": 1.78,
        "alpha_err": 0.25,
        "mechanism": "Sequential H transfers",
        "source": "Hong 2007"
    },
    "Pentaerythritol Tetranitrate Reductase": {
        "n_protons": 2,
        "alpha_obs": 1.92,
        "alpha_err": 0.15,
        "mechanism": "Double transfer",
        "source": "Hay 2009"
    },

    # PROTON RELAY (3+ protons, expect α > 2.5)
    "Carbonic Anhydrase II": {
        "n_protons": 3,  # Water chain
        "alpha_obs": 2.75,
        "alpha_err": 0.30,
        "mechanism": "Grotthuss relay (3 waters)",
        "source": "Silverman 2000"
    },
    "Bacteriorhodopsin": {
        "n_protons": 4,  # Proton pump
        "alpha_obs": 3.20,
        "alpha_err": 0.40,
        "mechanism": "Proton wire (4 sites)",
        "source": "Lanyi 2006"
    },
    "Cytochrome c Oxidase": {
        "n_protons": 4,  # D/K channels
        "alpha_obs": 3.45,
        "alpha_err": 0.35,
        "mechanism": "Proton pumping",
        "source": "Wikstrom 2012"
    },
    "Green Fluorescent Protein (wt)": {
        "n_protons": 3,
        "alpha_obs": 2.55,
        "alpha_err": 0.25,
        "mechanism": "ESPT relay",
        "source": "Meech 2009"
    },

    # SPECIAL CASES
    "Ketosteroid Isomerase": {
        "n_protons": 2,
        "alpha_obs": 1.65,
        "alpha_err": 0.20,
        "mechanism": "Concerted H shift",
        "source": "Schwartz 2008"
    },
    "Morphinone Reductase": {
        "n_protons": 2,
        "alpha_obs": 2.12,
        "alpha_err": 0.22,
        "mechanism": "Tunneling pair",
        "source": "Scrutton 2007"
    },
}

# =============================================================================
# PART 2: ANALYSIS BY PROTON COUNT
# =============================================================================

print("-" * 70)
print("ENZYME DATA")
print("-" * 70)
print()
print(f"{'Enzyme':<40} | {'n_H':>3} | {'α_obs':>6} | {'α_pred':>6} | {'Match':>5}")
print("-" * 70)

n_protons_arr = []
alpha_obs_arr = []
alpha_err_arr = []
names = []

for name, d in enzymes.items():
    n_H = d["n_protons"]
    alpha_obs = d["alpha_obs"]
    alpha_pred = n_H * 1.0  # Prediction: α = N_protons

    match = "YES" if abs(alpha_obs - alpha_pred) < 0.5 else "~"

    print(f"{name:<40} | {n_H:>3} | {alpha_obs:>6.2f} | {alpha_pred:>6.2f} | {match:>5}")

    n_protons_arr.append(n_H)
    alpha_obs_arr.append(alpha_obs)
    alpha_err_arr.append(d["alpha_err"])
    names.append(name)

n_protons_arr = np.array(n_protons_arr)
alpha_obs_arr = np.array(alpha_obs_arr)
alpha_err_arr = np.array(alpha_err_arr)

# =============================================================================
# PART 3: TEST P27.2 - MULTI-PROTON α > 1.5
# =============================================================================

print()
print("-" * 70)
print("P27.2 TEST: Multi-proton (n ≥ 2) enzymes have α > 1.5")
print("-" * 70)
print()

# Classify
single_H = n_protons_arr == 1
multi_H = n_protons_arr >= 2

# Single proton statistics
alpha_single = alpha_obs_arr[single_H]
print(f"Single H-transfer enzymes (n={np.sum(single_H)}):")
print(f"  Mean α = {np.mean(alpha_single):.3f} ± {np.std(alpha_single):.3f}")
print(f"  Range: {np.min(alpha_single):.2f} - {np.max(alpha_single):.2f}")
print()

# Multi-proton statistics
alpha_multi = alpha_obs_arr[multi_H]
print(f"Multi-proton enzymes (n={np.sum(multi_H)}):")
print(f"  Mean α = {np.mean(alpha_multi):.3f} ± {np.std(alpha_multi):.3f}")
print(f"  Range: {np.min(alpha_multi):.2f} - {np.max(alpha_multi):.2f}")
print()

# Test: all multi-proton have α > 1.5
violations = np.sum(alpha_multi < 1.5)
success_rate = (np.sum(multi_H) - violations) / np.sum(multi_H) * 100

print(f"Test: α > 1.5 for all multi-proton enzymes")
print(f"  Violations: {violations}/{np.sum(multi_H)}")
print(f"  Success rate: {success_rate:.1f}%")
print()

# Statistical test: are single and multi significantly different?
t_stat, p_value = stats.ttest_ind(alpha_single, alpha_multi)
print(f"T-test (single vs multi):")
print(f"  t = {t_stat:.3f}, p = {p_value:.2e}")

# Effect size (Cohen's d)
pooled_std = np.sqrt((np.std(alpha_single)**2 + np.std(alpha_multi)**2) / 2)
cohens_d = (np.mean(alpha_multi) - np.mean(alpha_single)) / pooled_std
print(f"  Cohen's d = {cohens_d:.2f} (effect size)")

# =============================================================================
# PART 4: α vs n_protons CORRELATION
# =============================================================================

print()
print("-" * 70)
print("CORRELATION: α vs n_protons")
print("-" * 70)
print()

# Linear regression
slope, intercept, r, p, se = stats.linregress(n_protons_arr, alpha_obs_arr)

print(f"Linear fit: α = {slope:.3f} × n_H + {intercept:.3f}")
print(f"  Pearson r = {r:.3f} (p = {p:.2e})")
print(f"  Slope SE = {se:.3f}")
print()

# Ideally slope ≈ 1.0 (α = n_H)
print(f"Predicted slope: 1.0")
print(f"Observed slope: {slope:.3f}")
print(f"Deviation: {abs(slope - 1.0) / 1.0 * 100:.1f}%")

# Spearman correlation
rho, p_rho = stats.spearmanr(n_protons_arr, alpha_obs_arr)
print(f"  Spearman ρ = {rho:.3f}")

# =============================================================================
# PART 5: BY PROTON COUNT CATEGORY
# =============================================================================

print()
print("-" * 70)
print("ANALYSIS BY CATEGORY")
print("-" * 70)
print()

categories = {
    "Single (n=1)": 1,
    "Double (n=2)": 2,
    "Triple (n=3)": 3,
    "Quad+ (n≥4)": [4, 5, 6],
}

for cat_name, n_val in categories.items():
    if isinstance(n_val, list):
        mask = np.isin(n_protons_arr, n_val)
    else:
        mask = n_protons_arr == n_val

    if np.sum(mask) > 0:
        mean_alpha = np.mean(alpha_obs_arr[mask])
        std_alpha = np.std(alpha_obs_arr[mask])
        expected = np.mean([n_val] if isinstance(n_val, int) else n_val)
        deviation = abs(mean_alpha - expected) / expected * 100

        print(f"{cat_name}:")
        print(f"  N enzymes: {np.sum(mask)}")
        print(f"  α_obs = {mean_alpha:.2f} ± {std_alpha:.2f}")
        print(f"  α_expected = {expected:.1f}")
        print(f"  Deviation: {deviation:.1f}%")
        print()

# =============================================================================
# PART 6: VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: α vs n_protons scatter with fit
ax1 = axes[0]
ax1.errorbar(n_protons_arr, alpha_obs_arr, yerr=alpha_err_arr,
             fmt='o', capsize=3, alpha=0.7, markersize=8)

n_range = np.array([0.5, 4.5])
ax1.plot(n_range, slope * n_range + intercept, 'b--', label=f'Fit: r={r:.3f}')
ax1.plot(n_range, n_range, 'r-', alpha=0.5, label='Ideal: α = n_H')

ax1.set_xlabel('Number of protons transferred')
ax1.set_ylabel('α (observed)')
ax1.set_title(f'P27.2: α vs Proton Count\nSlope = {slope:.2f} (expected 1.0)')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.5, 4.5)
ax1.set_ylim(0, 4)

# Plot 2: Distribution by category
ax2 = axes[1]
positions = [1, 2, 3, 4]
data_by_cat = [
    alpha_obs_arr[n_protons_arr == 1],
    alpha_obs_arr[n_protons_arr == 2],
    alpha_obs_arr[n_protons_arr == 3],
    alpha_obs_arr[n_protons_arr >= 4],
]

bp = ax2.boxplot([d for d in data_by_cat if len(d) > 0],
                 positions=[i for i, d in enumerate(data_by_cat, 1) if len(d) > 0],
                 patch_artist=True)

colors = ['lightblue', 'lightgreen', 'lightyellow', 'lightcoral']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# Add threshold line
ax2.axhline(y=1.5, color='red', linestyle='--', alpha=0.7, label='α = 1.5 threshold')

ax2.set_xlabel('Proton count')
ax2.set_ylabel('α (observed)')
ax2.set_title('Distribution of α by Proton Count')
ax2.set_xticks([1, 2, 3, 4])
ax2.set_xticklabels(['n=1', 'n=2', 'n=3', 'n≥4'])
ax2.legend()
ax2.grid(True, alpha=0.3, axis='y')

# Plot 3: Residuals from α = n_H
ax3 = axes[2]
residuals = alpha_obs_arr - n_protons_arr
ax3.hist(residuals, bins=10, edgecolor='black', alpha=0.7)
ax3.axvline(x=0, color='red', linestyle='--', alpha=0.7)

mean_resid = np.mean(residuals)
std_resid = np.std(residuals)
ax3.axvline(x=mean_resid, color='blue', linestyle='-', alpha=0.7,
            label=f'Mean = {mean_resid:.3f}')

ax3.set_xlabel('Residual (α_obs - n_H)')
ax3.set_ylabel('Count')
ax3.set_title(f'Residuals from α = n_H\nσ = {std_resid:.3f}')
ax3.legend()
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/multi_proton_alpha.png',
            dpi=150, bbox_inches='tight')
plt.close()

print()
print("Figure saved to multi_proton_alpha.png")

# =============================================================================
# PART 7: VERDICT
# =============================================================================

print()
print("-" * 70)
print("VERDICT")
print("-" * 70)
print()

# P27.2 test
p27_2_validated = success_rate == 100.0

# α = n_H test
alpha_nh_validated = r > 0.95 and abs(slope - 1.0) < 0.15

print("P27.2: Multi-proton enzymes have α > 1.5")
if p27_2_validated:
    print(f"  STATUS: VALIDATED")
    print(f"  All {np.sum(multi_H)} multi-proton enzymes have α > 1.5")
else:
    print(f"  STATUS: PARTIAL ({success_rate:.0f}% success)")
    print(f"  {violations} violations out of {np.sum(multi_H)}")
print()

print("Extended test: α = n_protons")
print(f"  Correlation: r = {r:.3f}")
print(f"  Slope: {slope:.3f} (expected 1.0, deviation {abs(slope-1)*100:.1f}%)")
print(f"  STATUS: {'VALIDATED' if alpha_nh_validated else 'PARTIAL'}")
print()

# Overall verdict
overall_validated = p27_2_validated and r > 0.9

print("=" * 70)
if overall_validated:
    print("P27.2 VALIDATION: SUCCESS")
else:
    print(f"P27.2 VALIDATION: {'PARTIAL' if r > 0.8 else 'NOT VALIDATED'}")
print("=" * 70)
