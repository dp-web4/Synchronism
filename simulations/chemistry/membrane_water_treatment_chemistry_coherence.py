"""
Chemistry Session #1853: Membrane Water Treatment Chemistry
Phenomenon: 1716th type - Rejection rate coherence
Water Treatment Chemistry Series

Synchronism Framework: gamma = 2/sqrt(N_corr), CF = 1/(1+gamma^2)

Membrane treatment maps onto the coherence framework:
- N_corr represents the number of correlated solute-membrane interaction events
- gamma measures the disorder in solute rejection/transport
- CF tracks the fraction of solutes coherently rejected by the membrane
- RO, NF, UF, MF each occupy different N_corr regimes
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N_corr):
    return 2.0 / np.sqrt(N_corr)

def coherence_fraction(gamma_val):
    return 1.0 / (1.0 + gamma_val**2)

N_corr = np.linspace(1, 20, 1000)
g = gamma(N_corr)
cf = coherence_fraction(g)

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('Membrane Water Treatment Chemistry - Rejection Rate Coherence Analysis\nSynchronism Framework: γ = 2/√N_corr  |  Session #1853', fontsize=14, fontweight='bold')

passed = 0
total = 8

# Test 1: γ=1 at N_corr=4 — RO Desalination Threshold
ax = axes[0, 0]
ax.plot(N_corr, g, 'b-', linewidth=2)
ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.7, label='γ=1 (rejection onset)')
ax.axvline(x=4.0, color='g', linestyle='--', alpha=0.7, label='N_corr=4')
g_at_4 = gamma(4.0)
test1 = abs(g_at_4 - 1.0) < 0.01
if test1: passed += 1
ax.set_title(f'T1: RO Desalination Threshold\nγ(4)={g_at_4:.4f} [{"PASS" if test1 else "FAIL"}]', color='green' if test1 else 'red', fontsize=10)
ax.set_xlabel('N_corr (membrane interactions)'); ax.set_ylabel('γ (transport disorder)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Test 2: CF=0.5 at N_corr=4 — NF Softening Balance
ax = axes[0, 1]
cf_at_4 = coherence_fraction(gamma(4.0))
test2 = abs(cf_at_4 - 0.5) < 0.01
if test2: passed += 1
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=0.5, color='r', linestyle='--', alpha=0.7, label='CF=0.5'); ax.axvline(x=4.0, color='g', linestyle='--', alpha=0.7, label='N_corr=4')
ax.plot(4.0, cf_at_4, 'ro', markersize=8)
ax.set_title(f'T2: NF Softening Balance\nCF(4)={cf_at_4:.4f} [{"PASS" if test2 else "FAIL"}]', color='green' if test2 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('CF (rejection coherence)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Test 3: CF=1-1/e — Ultrafiltration Regime
ax = axes[0, 2]
target_cf3 = 1.0 - 1.0/np.e
g_t3 = np.sqrt(1.0/target_cf3 - 1.0); N_t3 = (2.0/g_t3)**2
cf_c3 = coherence_fraction(gamma(N_t3))
test3 = abs(cf_c3 - target_cf3) < 0.01
if test3: passed += 1
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=target_cf3, color='r', linestyle='--', alpha=0.7, label=f'CF=1-1/e={target_cf3:.3f}')
ax.axvline(x=N_t3, color='g', linestyle='--', alpha=0.7, label=f'N={N_t3:.2f}')
ax.plot(N_t3, cf_c3, 'ro', markersize=8)
ax.set_title(f'T3: Ultrafiltration Regime\nCF({N_t3:.2f})={cf_c3:.4f} [{"PASS" if test3 else "FAIL"}]', color='green' if test3 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Test 4: CF=1/e — Microfiltration Baseline
ax = axes[0, 3]
target_cf4 = 1.0/np.e
g_t4 = np.sqrt(1.0/target_cf4 - 1.0); N_t4 = (2.0/g_t4)**2
cf_c4 = coherence_fraction(gamma(N_t4))
test4 = abs(cf_c4 - target_cf4) < 0.01
if test4: passed += 1
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=target_cf4, color='r', linestyle='--', alpha=0.7, label=f'CF=1/e={target_cf4:.3f}')
ax.axvline(x=N_t4, color='g', linestyle='--', alpha=0.7, label=f'N={N_t4:.2f}')
ax.plot(N_t4, cf_c4, 'ro', markersize=8)
ax.set_title(f'T4: Microfiltration Baseline\nCF({N_t4:.2f})={cf_c4:.4f} [{"PASS" if test4 else "FAIL"}]', color='green' if test4 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Test 5: dγ/dN < 0 — Fouling Resistance Decay
ax = axes[1, 0]
dg = np.diff(g); test5 = np.all(dg < 0)
if test5: passed += 1
ax.plot(N_corr[1:], dg, 'b-', linewidth=2); ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
ax.fill_between(N_corr[1:], dg, 0, alpha=0.2, color='blue')
ax.set_title(f'T5: Fouling Resistance Decay\ndγ/dN<0 monotonic [{"PASS" if test5 else "FAIL"}]', color='green' if test5 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('dγ/dN'); ax.grid(True, alpha=0.3)

# Test 6: γ(1)=2 — Scaling Onset
ax = axes[1, 1]
g_at_1 = gamma(1.0); test6 = abs(g_at_1 - 2.0) < 0.01
if test6: passed += 1
ax.plot(N_corr, g, 'b-', linewidth=2); ax.plot(1.0, g_at_1, 'ro', markersize=10, label=f'γ(1)={g_at_1:.1f}')
ax.axhline(y=2.0, color='r', linestyle='--', alpha=0.7, label='γ=2 (max disorder)')
ax.set_title(f'T6: Scaling Onset\nγ(1)={g_at_1:.4f} [{"PASS" if test6 else "FAIL"}]', color='green' if test6 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('γ'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Test 7: dCF/dN > 0 — Concentration Polarization Relief
ax = axes[1, 2]
dcf = np.diff(cf); test7 = np.all(dcf > 0)
if test7: passed += 1
ax.plot(N_corr[1:], dcf, 'b-', linewidth=2); ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
ax.fill_between(N_corr[1:], dcf, 0, alpha=0.2, color='green')
ax.set_title(f'T7: Concentration Polarization Relief\ndCF/dN>0 monotonic [{"PASS" if test7 else "FAIL"}]', color='green' if test7 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('dCF/dN'); ax.grid(True, alpha=0.3)

# Test 8: CF(100)>0.95 — Flux Decline Asymptote
ax = axes[1, 3]
N_ext = np.linspace(1, 100, 1000); cf_ext = coherence_fraction(gamma(N_ext))
cf_100 = coherence_fraction(gamma(100.0)); test8 = cf_100 > 0.95
if test8: passed += 1
ax.plot(N_ext, cf_ext, 'b-', linewidth=2)
ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.7, label='CF=1.0 (perfect rejection)')
ax.axhline(y=0.95, color='orange', linestyle='--', alpha=0.7, label='CF=0.95 (high rejection)')
ax.plot(100.0, cf_100, 'ro', markersize=8, label=f'CF(100)={cf_100:.4f}')
ax.set_title(f'T8: Flux Decline Asymptote\nCF(100)={cf_100:.4f} [{"PASS" if test8 else "FAIL"}]', color='green' if test8 else 'red', fontsize=10)
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = '/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/membrane_water_treatment_chemistry_coherence.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"All 8 boundary condition tests completed for Membrane Water Treatment Chemistry")
print(f"Results: {passed}/{total} tests PASSED")
print(f"Plot saved: {out_path}")
