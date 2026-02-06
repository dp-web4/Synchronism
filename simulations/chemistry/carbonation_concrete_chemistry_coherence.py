#!/usr/bin/env python3
"""
Chemistry Session #1866 - Carbonation Chemistry Coherence
1729th Phenomenon in Synchronism Framework

CO2 diffusion coherence in concrete carbonation:
  Ca(OH)2 + CO2 -> CaCO3 + H2O

The carbonation front advances as CO2 diffuses through the pore network.
Coherence emerges when CO2 molecular transport couples with the reaction
zone geometry, creating a phase-locked diffusion-reaction front.

Coherence metrics:
  gamma = 2 / sqrt(N_corr)    -- coherence parameter
  CF = 1 / (1 + gamma^2)      -- coherence fraction

where N_corr = number of correlated CO2-pore interaction domains.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N):
    return 2.0 / np.sqrt(N)

def CF(N):
    g = gamma(N)
    return 1.0 / (1.0 + g**2)

def dGamma_dN(N, dN=1e-6):
    return (gamma(N + dN) - gamma(N - dN)) / (2 * dN)

def dCF_dN(N, dN=1e-6):
    return (CF(N + dN) - CF(N - dN)) / (2 * dN)

# --- 8 Boundary Tests ---
tests = []

# T1: gamma = 1 at N = 4
g4 = gamma(4)
tests.append(("T1: γ=1 @ N=4", np.isclose(g4, 1.0, atol=1e-10), f"γ(4) = {g4:.10f}"))

# T2: CF = 0.5 at N = 4
cf4 = CF(4)
tests.append(("T2: CF=0.5 @ N=4", np.isclose(cf4, 0.5, atol=1e-10), f"CF(4) = {cf4:.10f}"))

# T3: CF = 1 - 1/e threshold
target3 = 1 - 1/np.e
# Solve: 1/(1+4/N) = 1-1/e => N = 4/((1/(1-1/e))-1) = 4*(e-1)/(1) = 4*(e-1)/1
# Actually: 1/(1+4/N)=1-1/e => 1+4/N = e/(e-1) => 4/N = 1/(e-1) => N=4(e-1)
N_t3 = 4*(np.e - 1)
cf_t3 = CF(N_t3)
tests.append(("T3: CF=1-1/e threshold", np.isclose(cf_t3, target3, atol=1e-10), f"CF({N_t3:.4f}) = {cf_t3:.10f}, target = {target3:.10f}"))

# T4: CF = 1/e threshold
target4 = 1/np.e
# Solve: 1/(1+4/N)=1/e => 1+4/N=e => 4/N=e-1 => N=4/(e-1)
N_t4 = 4/(np.e - 1)
cf_t4 = CF(N_t4)
tests.append(("T4: CF=1/e threshold", np.isclose(cf_t4, target4, atol=1e-10), f"CF({N_t4:.4f}) = {cf_t4:.10f}, target = {target4:.10f}"))

# T5: dγ/dN < 0 (gamma decreases with N)
dg = dGamma_dN(10)
tests.append(("T5: dγ/dN < 0", dg < 0, f"dγ/dN(10) = {dg:.10f}"))

# T6: gamma(1) = 2
g1 = gamma(1)
tests.append(("T6: γ(1)=2", np.isclose(g1, 2.0, atol=1e-10), f"γ(1) = {g1:.10f}"))

# T7: dCF/dN > 0 (CF increases with N)
dcf = dCF_dN(10)
tests.append(("T7: dCF/dN > 0", dcf > 0, f"dCF/dN(10) = {dcf:.10f}"))

# T8: CF(100) > 0.95
cf100 = CF(100)
tests.append(("T8: CF(100)>0.95", cf100 > 0.95, f"CF(100) = {cf100:.10f}"))

# --- Plotting ---
N_range = np.linspace(1, 120, 1000)

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle("Chemistry Session #1866 — Carbonation Concrete Chemistry Coherence\n"
             "1729th Phenomenon · CO₂ Diffusion Coherence · γ = 2/√N_corr, CF = 1/(1+γ²)",
             fontsize=13, fontweight='bold')

for idx, (title, passed, detail) in enumerate(tests):
    ax = axes[idx // 4, idx % 4]
    color = 'green' if passed else 'red'
    status = 'PASS' if passed else 'FAIL'

    if idx == 0:  # T1
        ax.plot(N_range, gamma(N_range), 'b-', lw=2)
        ax.axhline(1.0, color='gray', ls='--', alpha=0.5)
        ax.axvline(4.0, color='gray', ls='--', alpha=0.5)
        ax.plot(4, 1.0, 'o', color=color, ms=12, zorder=5)
        ax.set_ylabel('γ')
        ax.set_xlabel('N_corr')
    elif idx == 1:  # T2
        ax.plot(N_range, CF(N_range), 'b-', lw=2)
        ax.axhline(0.5, color='gray', ls='--', alpha=0.5)
        ax.axvline(4.0, color='gray', ls='--', alpha=0.5)
        ax.plot(4, 0.5, 'o', color=color, ms=12, zorder=5)
        ax.set_ylabel('CF')
        ax.set_xlabel('N_corr')
    elif idx == 2:  # T3
        ax.plot(N_range, CF(N_range), 'b-', lw=2)
        ax.axhline(target3, color='orange', ls='--', alpha=0.7, label=f'1-1/e={target3:.4f}')
        ax.axvline(N_t3, color='gray', ls='--', alpha=0.5)
        ax.plot(N_t3, target3, 'o', color=color, ms=12, zorder=5)
        ax.legend(fontsize=8)
        ax.set_ylabel('CF')
        ax.set_xlabel('N_corr')
    elif idx == 3:  # T4
        ax.plot(N_range, CF(N_range), 'b-', lw=2)
        ax.axhline(target4, color='orange', ls='--', alpha=0.7, label=f'1/e={target4:.4f}')
        ax.axvline(N_t4, color='gray', ls='--', alpha=0.5)
        ax.plot(N_t4, target4, 'o', color=color, ms=12, zorder=5)
        ax.legend(fontsize=8)
        ax.set_ylabel('CF')
        ax.set_xlabel('N_corr')
    elif idx == 4:  # T5
        N_d = np.linspace(1, 50, 500)
        dg_vals = np.array([dGamma_dN(n) for n in N_d])
        ax.plot(N_d, dg_vals, 'b-', lw=2)
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
        ax.fill_between(N_d, dg_vals, 0, where=(dg_vals < 0), alpha=0.2, color=color)
        ax.set_ylabel('dγ/dN')
        ax.set_xlabel('N_corr')
    elif idx == 5:  # T6
        ax.plot(N_range, gamma(N_range), 'b-', lw=2)
        ax.axhline(2.0, color='gray', ls='--', alpha=0.5)
        ax.axvline(1.0, color='gray', ls='--', alpha=0.5)
        ax.plot(1, 2.0, 'o', color=color, ms=12, zorder=5)
        ax.set_ylabel('γ')
        ax.set_xlabel('N_corr')
    elif idx == 6:  # T7
        N_d = np.linspace(1, 50, 500)
        dcf_vals = np.array([dCF_dN(n) for n in N_d])
        ax.plot(N_d, dcf_vals, 'b-', lw=2)
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
        ax.fill_between(N_d, dcf_vals, 0, where=(dcf_vals > 0), alpha=0.2, color=color)
        ax.set_ylabel('dCF/dN')
        ax.set_xlabel('N_corr')
    elif idx == 7:  # T8
        ax.plot(N_range, CF(N_range), 'b-', lw=2)
        ax.axhline(0.95, color='gray', ls='--', alpha=0.5)
        ax.axvline(100, color='gray', ls='--', alpha=0.5)
        ax.plot(100, cf100, 'o', color=color, ms=12, zorder=5)
        ax.set_ylabel('CF')
        ax.set_xlabel('N_corr')

    ax.set_title(f"{title}\n[{status}] {detail}", fontsize=9, color=color)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/carbonation_concrete_chemistry_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("=" * 70)
print("Chemistry Session #1866 — Carbonation Concrete Chemistry Coherence")
print("1729th Phenomenon — CO2 Diffusion Coherence")
print("=" * 70)
for title, passed, detail in tests:
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] {title}: {detail}")
print("=" * 70)
all_pass = all(p for _, p, _ in tests)
print(f"Overall: {'ALL 8 TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
print(f"Plot saved: carbonation_concrete_chemistry_coherence.png")
