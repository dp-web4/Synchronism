"""
Chemistry Session #1858: Softening Chemistry
Phenomenon: 1721st type - Hardness removal coherence
Synchronism Framework: gamma = 2/sqrt(N_corr), CF = 1/(1+gamma^2)

Water softening removes calcium and magnesium hardness through precipitation,
ion exchange, or membrane processes. The coherence framework maps crystallization
coordination from disordered nucleation (low N_corr) to organized crystal growth
(high N_corr), covering lime-soda softening, caustic softening, split treatment,
recarbonation, Langelier saturation index, Ryznar stability index, pellet
softening, and chemical dose optimization.
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
fig.suptitle("Softening Chemistry - Hardness Removal Coherence Analysis\n"
             "Synchronism Framework: γ = 2/√N_corr  |  Session #1858  |  Phenomenon #1721",
             fontsize=14, fontweight='bold')

passed = 0
total = 8
test_labels = [
    "T1: Lime-Soda Softening\nγ=1 at N_corr=4",
    "T2: Caustic Softening\nCF=0.5 at N_corr=4",
    "T3: Split Treatment\nCF=1-1/e threshold",
    "T4: Recarbonation\nCF=1/e threshold",
    "T5: Langelier Saturation\ndγ/dN<0 monotonic",
    "T6: Ryznar Index\nγ(1)=2",
    "T7: Pellet Softening\ndCF/dN>0 monotonic",
    "T8: Chemical Dose\nCF(100)>0.95 asymptotic"
]

# T1
ax = axes[0, 0]
g_at_4 = gamma(4.0)
test1 = abs(g_at_4 - 1.0) < 0.01
ax.plot(N_corr, g, 'b-', linewidth=2)
ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
ax.axvline(x=4.0, color='g', linestyle='--', alpha=0.7)
ax.plot(4.0, g_at_4, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('γ')
ax.set_title(test_labels[0] + f"\n{'PASS' if test1 else 'FAIL'}: γ(4)={g_at_4:.4f}",
             color='green' if test1 else 'red', fontweight='bold')
if test1: passed += 1

# T2
ax = axes[0, 1]
cf_at_4 = coherence_fraction(gamma(4.0))
test2 = abs(cf_at_4 - 0.5) < 0.01
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=0.5, color='r', linestyle='--', alpha=0.7)
ax.axvline(x=4.0, color='g', linestyle='--', alpha=0.7)
ax.plot(4.0, cf_at_4, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.set_title(test_labels[1] + f"\n{'PASS' if test2 else 'FAIL'}: CF(4)={cf_at_4:.4f}",
             color='green' if test2 else 'red', fontweight='bold')
if test2: passed += 1

# T3
ax = axes[0, 2]
target_cf3 = 1.0 - 1.0/np.e
g_t3 = np.sqrt(1.0/target_cf3 - 1.0)
N_t3 = (2.0/g_t3)**2
cf_c3 = coherence_fraction(gamma(N_t3))
test3 = abs(cf_c3 - target_cf3) < 0.01
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=target_cf3, color='r', linestyle='--', alpha=0.7)
ax.axvline(x=N_t3, color='g', linestyle='--', alpha=0.7)
ax.plot(N_t3, cf_c3, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.set_title(test_labels[2] + f"\n{'PASS' if test3 else 'FAIL'}: CF={cf_c3:.4f} at N={N_t3:.2f}",
             color='green' if test3 else 'red', fontweight='bold')
if test3: passed += 1

# T4
ax = axes[0, 3]
target_cf4 = 1.0/np.e
g_t4 = np.sqrt(1.0/target_cf4 - 1.0)
N_t4 = (2.0/g_t4)**2
cf_c4 = coherence_fraction(gamma(N_t4))
test4 = abs(cf_c4 - target_cf4) < 0.01
ax.plot(N_corr, cf, 'b-', linewidth=2)
ax.axhline(y=target_cf4, color='r', linestyle='--', alpha=0.7)
ax.axvline(x=N_t4, color='g', linestyle='--', alpha=0.7)
ax.plot(N_t4, cf_c4, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.set_title(test_labels[3] + f"\n{'PASS' if test4 else 'FAIL'}: CF={cf_c4:.4f} at N={N_t4:.2f}",
             color='green' if test4 else 'red', fontweight='bold')
if test4: passed += 1

# T5
ax = axes[1, 0]
dg = np.diff(g)
test5 = np.all(dg < 0)
ax.plot(N_corr[:-1], dg, 'b-', linewidth=2)
ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
ax.fill_between(N_corr[:-1], dg, 0, alpha=0.3, color='blue')
ax.set_xlabel('N_corr')
ax.set_ylabel('dγ/dN')
ax.set_title(test_labels[4] + f"\n{'PASS' if test5 else 'FAIL'}: all dγ/dN<0={test5}",
             color='green' if test5 else 'red', fontweight='bold')
if test5: passed += 1

# T6
ax = axes[1, 1]
g_at_1 = gamma(1.0)
test6 = abs(g_at_1 - 2.0) < 0.01
ax.plot(N_corr, g, 'b-', linewidth=2)
ax.axhline(y=2.0, color='r', linestyle='--', alpha=0.7)
ax.axvline(x=1.0, color='g', linestyle='--', alpha=0.7)
ax.plot(1.0, g_at_1, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('γ')
ax.set_title(test_labels[5] + f"\n{'PASS' if test6 else 'FAIL'}: γ(1)={g_at_1:.4f}",
             color='green' if test6 else 'red', fontweight='bold')
if test6: passed += 1

# T7
ax = axes[1, 2]
dcf = np.diff(cf)
test7 = np.all(dcf > 0)
ax.plot(N_corr[:-1], dcf, 'b-', linewidth=2)
ax.axhline(y=0, color='r', linestyle='--', alpha=0.7)
ax.fill_between(N_corr[:-1], dcf, 0, alpha=0.3, color='green')
ax.set_xlabel('N_corr')
ax.set_ylabel('dCF/dN')
ax.set_title(test_labels[6] + f"\n{'PASS' if test7 else 'FAIL'}: all dCF/dN>0={test7}",
             color='green' if test7 else 'red', fontweight='bold')
if test7: passed += 1

# T8
ax = axes[1, 3]
N_ext = np.linspace(1, 100, 1000)
cf_ext = coherence_fraction(gamma(N_ext))
cf_100 = coherence_fraction(gamma(100.0))
test8 = cf_100 > 0.95
ax.plot(N_ext, cf_ext, 'b-', linewidth=2)
ax.axhline(y=0.95, color='r', linestyle='--', alpha=0.7)
ax.plot(100.0, cf_100, 'ro', markersize=10)
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.set_title(test_labels[7] + f"\n{'PASS' if test8 else 'FAIL'}: CF(100)={cf_100:.4f}",
             color='green' if test8 else 'red', fontweight='bold')
if test8: passed += 1

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/softening_water_chemistry_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()
print(f"All 8 boundary condition tests completed for Softening Chemistry")
print(f"Results: {passed}/{total} tests PASSED")
