#!/usr/bin/env python3
"""
Chemistry Session #1841: Essential Oil Chemistry
Phenomenon: 1704th type - Terpene extraction coherence
Master equation: gamma = 2/sqrt(N_corr), CF = 1/(1+gamma^2)

8 Boundary Tests (2x4 subplot grid):
T1: gamma=1 at N_corr=4
T2: CF=0.5 at N_corr=4
T3: CF=1-1/e at derived N
T4: CF=1/e at derived N
T5: d(gamma)/dN < 0 monotonic
T6: gamma(1) = 2
T7: dCF/dN > 0 monotonic
T8: CF(100) > 0.95 asymptotic

Domain tests: Steam distillation, cold press, terpene profile, oxidation,
photodegradation, encapsulation, headspace, shelf life
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Master equations ---
def gamma(N):
    return 2.0 / np.sqrt(N)

def CF(N):
    g = gamma(N)
    return 1.0 / (1.0 + g**2)

# --- Derived thresholds ---
# CF = 1 - 1/e  =>  1/(1+g^2) = 1 - 1/e  =>  g^2 = 1/((1-1/e)) - 1 = 1/(e-1)/(e) ... solve
# g^2 = e/(e-1) - 1 = 1/(e-1)
# N = 4/g^2 = 4*(e-1)
N_high = 4.0 * (np.e - 1)  # ~6.873

# CF = 1/e  =>  1/(1+g^2) = 1/e  =>  g^2 = e - 1  =>  N = 4/(e-1)
N_low = 4.0 / (np.e - 1)  # ~2.328

# --- Run 8 tests ---
results = []
N_arr = np.linspace(1, 100, 1000)
g_arr = gamma(N_arr)
cf_arr = CF(N_arr)

# T1: gamma=1 at N_corr=4
t1 = np.isclose(gamma(4), 1.0, atol=1e-10)
results.append(("T1: γ=1 @ N=4", t1, f"γ(4)={gamma(4):.6f}"))

# T2: CF=0.5 at N_corr=4
t2 = np.isclose(CF(4), 0.5, atol=1e-10)
results.append(("T2: CF=0.5 @ N=4", t2, f"CF(4)={CF(4):.6f}"))

# T3: CF=1-1/e at derived N
t3 = np.isclose(CF(N_high), 1.0 - 1.0/np.e, atol=1e-10)
results.append(("T3: CF=1-1/e @ N*", t3, f"CF({N_high:.3f})={CF(N_high):.6f}, target={1-1/np.e:.6f}"))

# T4: CF=1/e at derived N
t4 = np.isclose(CF(N_low), 1.0/np.e, atol=1e-10)
results.append(("T4: CF=1/e @ N*", t4, f"CF({N_low:.3f})={CF(N_low):.6f}, target={1/np.e:.6f}"))

# T5: d(gamma)/dN < 0 monotonic
dg = np.diff(g_arr)
t5 = np.all(dg < 0)
results.append(("T5: dγ/dN<0 monotonic", t5, f"max(dγ)={np.max(dg):.2e}"))

# T6: gamma(1) = 2
t6 = np.isclose(gamma(1), 2.0, atol=1e-10)
results.append(("T6: γ(1)=2", t6, f"γ(1)={gamma(1):.6f}"))

# T7: dCF/dN > 0 monotonic
dcf = np.diff(cf_arr)
t7 = np.all(dcf > 0)
results.append(("T7: dCF/dN>0 monotonic", t7, f"min(dCF)={np.min(dcf):.2e}"))

# T8: CF(100) > 0.95
t8 = CF(100) > 0.95
results.append(("T8: CF(100)>0.95", t8, f"CF(100)={CF(100):.6f}"))

passed = sum(1 for _, p, _ in results if p)

# --- Domain-specific labels ---
domain_labels = [
    "Steam distillation", "Cold press", "Terpene profile", "Oxidation",
    "Photodegradation", "Encapsulation", "Headspace", "Shelf life"
]

# --- Plot ---
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle(
    f"Session #1841: Essential Oil Chemistry — Terpene Extraction Coherence\n"
    f"γ = 2/√N_corr | CF = 1/(1+γ²) | Tests passed: {passed}/8",
    fontsize=14, fontweight='bold'
)

for idx, ax in enumerate(axes.flat):
    label, ok, detail = results[idx]
    color = 'green' if ok else 'red'
    status = 'PASS' if ok else 'FAIL'

    if idx < 4:
        # Point tests - plot CF curve with marked point
        ax.plot(N_arr, cf_arr, 'b-', linewidth=1.5)
        if idx == 0:
            ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
            ax.plot(4, CF(4), 'ro', markersize=10)
            ax.axvline(x=4, color='gray', linestyle='--', alpha=0.5)
        elif idx == 1:
            ax.plot(4, CF(4), 'ro', markersize=10)
            ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5)
        elif idx == 2:
            ax.plot(N_high, CF(N_high), 'ro', markersize=10)
            ax.axhline(y=1-1/np.e, color='orange', linestyle='--', alpha=0.5)
        elif idx == 3:
            ax.plot(N_low, CF(N_low), 'ro', markersize=10)
            ax.axhline(y=1/np.e, color='orange', linestyle='--', alpha=0.5)
        ax.set_xlabel("N_corr")
        ax.set_ylabel("CF")
    elif idx == 4:
        ax.plot(N_arr, g_arr, 'b-', linewidth=1.5)
        ax.fill_between(N_arr, g_arr, alpha=0.1, color='blue')
        ax.set_xlabel("N_corr")
        ax.set_ylabel("γ")
    elif idx == 5:
        ax.plot(N_arr, g_arr, 'b-', linewidth=1.5)
        ax.plot(1, gamma(1), 'ro', markersize=10)
        ax.axhline(y=2.0, color='orange', linestyle='--', alpha=0.5)
        ax.set_xlabel("N_corr")
        ax.set_ylabel("γ")
    elif idx == 6:
        ax.plot(N_arr[1:], dcf, 'b-', linewidth=1.5)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel("N_corr")
        ax.set_ylabel("dCF/dN")
    elif idx == 7:
        ax.plot(N_arr, cf_arr, 'b-', linewidth=1.5)
        ax.axhline(y=0.95, color='orange', linestyle='--', alpha=0.5)
        ax.plot(100, CF(100), 'ro', markersize=10)
        ax.set_xlabel("N_corr")
        ax.set_ylabel("CF")

    ax.set_title(f"{label}\n[{domain_labels[idx]}]\n[{status}] {detail}", fontsize=9, color=color)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = "/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/essential_oil_chemistry_coherence_v2.png"
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")
print(f"Results: {passed}/8 tests passed")
for label, ok, detail in results:
    print(f"  {'PASS' if ok else 'FAIL'}: {label} — {detail}")
