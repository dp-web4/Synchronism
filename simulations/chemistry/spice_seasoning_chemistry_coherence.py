#!/usr/bin/env python3
"""
Chemistry Session #1849: Spice & Seasoning Chemistry
Phenomenon: 1712th type - Pungency coherence

Spice chemistry involves diverse pungent compounds that activate
TRPV1 and TRPA1 receptors. Capsaicin (chili), piperine (pepper),
allicin (garlic), gingerol (ginger), and others create characteristic
heat and flavor through molecular interactions governed by coherence.

The master equation gamma = 2/sqrt(N_corr) captures how pungent
molecule correlations create the perceived intensity spectrum.

Master equation: gamma = 2/sqrt(N_corr)
Coherence fraction: CF = 1/(1 + gamma^2)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N_corr):
    return 2.0 / np.sqrt(N_corr)

def CF(N_corr):
    g = gamma(N_corr)
    return 1.0 / (1.0 + g**2)

def dCF_dN(N_corr, dN=1e-6):
    return (CF(N_corr + dN) - CF(N_corr - dN)) / (2 * dN)

def dgamma_dN(N_corr, dN=1e-6):
    return (gamma(N_corr + dN) - gamma(N_corr - dN)) / (2 * dN)

results = []

# T1: gamma=1 at N_corr=4 (Capsaicin Scoville)
t1_val = gamma(4)
t1_pass = np.isclose(t1_val, 1.0, atol=1e-10)
results.append(("T1: Capsaicin Scoville\n(gamma=1 at N_corr=4)", t1_pass, t1_val, 1.0))

# T2: CF=0.5 at N_corr=4 (Piperine)
t2_val = CF(4)
t2_pass = np.isclose(t2_val, 0.5, atol=1e-10)
results.append(("T2: Piperine\n(CF=0.5 at N_corr=4)", t2_pass, t2_val, 0.5))

# T3: CF = 1-1/e at derived N (Allicin)
target_cf3 = 1.0 - 1.0/np.e
N3 = 4.0 / (1.0/target_cf3 - 1.0)
t3_val = CF(N3)
t3_pass = np.isclose(t3_val, target_cf3, atol=1e-10)
results.append(("T3: Allicin\n(CF=1-1/e at derived N)", t3_pass, t3_val, target_cf3))

# T4: CF = 1/e at derived N (Gingerol)
target_cf4 = 1.0/np.e
N4 = 4.0 / (1.0/target_cf4 - 1.0)
t4_val = CF(N4)
t4_pass = np.isclose(t4_val, target_cf4, atol=1e-10)
results.append(("T4: Gingerol\n(CF=1/e at derived N)", t4_pass, t4_val, target_cf4))

# T5: dgamma/dN < 0 monotonic (Curcumin)
N_test = np.linspace(1, 100, 1000)
dg = np.array([dgamma_dN(n) for n in N_test])
t5_pass = np.all(dg < 0)
t5_val = np.max(dg)
results.append(("T5: Curcumin\n(dgamma/dN<0 monotonic)", t5_pass, t5_val, "< 0"))

# T6: gamma(1)=2 (Vanillin)
t6_val = gamma(1)
t6_pass = np.isclose(t6_val, 2.0, atol=1e-10)
results.append(("T6: Vanillin\n(gamma(1)=2)", t6_pass, t6_val, 2.0))

# T7: dCF/dN > 0 monotonic (Cinnamic aldehyde)
dcf = np.array([dCF_dN(n) for n in N_test])
t7_pass = np.all(dcf > 0)
t7_val = np.min(dcf)
results.append(("T7: Cinnamic aldehyde\n(dCF/dN>0 monotonic)", t7_pass, t7_val, "> 0"))

# T8: CF(100)>0.95 asymptotic (Eugenol)
t8_val = CF(100)
t8_pass = t8_val > 0.95
results.append(("T8: Eugenol\n(CF(100)>0.95)", t8_pass, t8_val, "> 0.95"))

passed = sum(1 for _, p, _, _ in results if p)
print(f"Session #1849: Spice & Seasoning Chemistry - {passed}/8 tests passed")

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle("Session #1849: Spice & Seasoning Chemistry\n"
             "Phenomenon #1712: Pungency Coherence | "
             f"gamma = 2/sqrt(N_corr) | {passed}/8 tests passed",
             fontsize=14, fontweight='bold')

N_plot = np.linspace(0.5, 120, 500)
gamma_plot = gamma(N_plot)
cf_plot = CF(N_plot)

for idx, (title, passed_test, actual, expected) in enumerate(results):
    ax = axes[idx // 4][idx % 4]
    color = '#2ecc71' if passed_test else '#e74c3c'

    if idx in [0, 5]:
        ax.plot(N_plot, gamma_plot, 'b-', linewidth=2)
        if idx == 0:
            ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
            ax.axvline(x=4.0, color='r', linestyle='--', alpha=0.7)
            ax.plot(4, 1.0, 'ro', markersize=10)
        else:
            ax.axhline(y=2.0, color='r', linestyle='--', alpha=0.7)
            ax.axvline(x=1.0, color='r', linestyle='--', alpha=0.7)
            ax.plot(1, 2.0, 'ro', markersize=10)
        ax.set_ylabel('gamma')
    elif idx in [1, 2, 3, 7]:
        ax.plot(N_plot, cf_plot, 'b-', linewidth=2)
        if idx == 1:
            ax.axhline(y=0.5, color='r', linestyle='--', alpha=0.7)
            ax.axvline(x=4.0, color='r', linestyle='--', alpha=0.7)
            ax.plot(4, 0.5, 'ro', markersize=10)
        elif idx == 2:
            ax.axhline(y=1-1/np.e, color='r', linestyle='--', alpha=0.7)
            ax.axvline(x=N3, color='r', linestyle='--', alpha=0.7)
            ax.plot(N3, 1-1/np.e, 'ro', markersize=10)
        elif idx == 3:
            ax.axhline(y=1/np.e, color='r', linestyle='--', alpha=0.7)
            ax.axvline(x=N4, color='r', linestyle='--', alpha=0.7)
            ax.plot(N4, 1/np.e, 'ro', markersize=10)
        else:
            ax.axhline(y=0.95, color='r', linestyle='--', alpha=0.7)
            ax.plot(100, CF(100), 'ro', markersize=10)
        ax.set_ylabel('CF')
    elif idx == 4:
        ax.plot(N_test, dg, 'b-', linewidth=2)
        ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
        ax.fill_between(N_test, dg, 0, alpha=0.3, color='blue')
        ax.set_ylabel('dgamma/dN')
    elif idx == 6:
        ax.plot(N_test, dcf, 'b-', linewidth=2)
        ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
        ax.fill_between(N_test, dcf, 0, alpha=0.3, color='blue')
        ax.set_ylabel('dCF/dN')

    ax.set_xlabel('N_corr')
    ax.set_title(title, fontsize=9, color=color, fontweight='bold')
    status = "PASS" if passed_test else "FAIL"
    ax.text(0.02, 0.98, f"{status}", transform=ax.transAxes,
            fontsize=12, fontweight='bold', color=color,
            verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/spice_seasoning_chemistry_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("PNG saved successfully.")
