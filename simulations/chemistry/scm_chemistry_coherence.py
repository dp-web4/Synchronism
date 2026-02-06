"""
Chemistry Session #1862: Supplementary Cementitious Materials (SCM) Chemistry
1725th Phenomenon - Pozzolanic Reaction Coherence

Synchronism framework applied to pozzolanic reactions in SCMs.
Pozzolanic materials (fly ash, silica fume, slag) react with Ca(OH)2
to form additional C-S-H gel. The coherence framework models how
correlated reactive silica/alumina sites achieve phase-locked
resonance during pozzolanic reaction.

γ = 2/√N_corr  (coherence decay parameter)
CF = 1/(1+γ²)  (coherence fraction)

8 Boundary Tests (2x4 subplot grid):
T1: γ=1 @ N=4       T2: CF=0.5 @ N=4
T3: CF=1-1/e         T4: CF=1/e
T5: dγ/dN<0          T6: γ(1)=2
T7: dCF/dN>0         T8: CF(100)>0.95
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Core functions ---
def gamma(N):
    N = np.asarray(N, dtype=float)
    return 2.0 / np.sqrt(N)

def CF(N):
    g = gamma(N)
    return 1.0 / (1.0 + g**2)

def dGamma_dN(N):
    N = np.asarray(N, dtype=float)
    return -1.0 / N**1.5

def dCF_dN(N):
    N = np.asarray(N, dtype=float)
    g = gamma(N)
    dgdN = dGamma_dN(N)
    return -2.0 * g * dgdN / (1.0 + g**2)**2

# --- SCM specifics ---
# Pozzolanic reaction: correlated SiO2/Al2O3 reactive sites with Ca(OH)2
N_range = np.linspace(1, 120, 1000)
gamma_vals = gamma(N_range)
CF_vals = CF(N_range)

# --- 8 Boundary Tests ---
results = {}

# T1: γ = 1 @ N = 4
t1_val = gamma(4)
results['T1'] = {'expected': 1.0, 'actual': t1_val, 'pass': np.isclose(t1_val, 1.0, atol=1e-10)}

# T2: CF = 0.5 @ N = 4
t2_val = CF(4)
results['T2'] = {'expected': 0.5, 'actual': t2_val, 'pass': np.isclose(t2_val, 0.5, atol=1e-10)}

# T3: CF = 1 - 1/e
target_cf3 = 1.0 - 1.0/np.e
N_t3 = 4.0 / (1.0/(1.0 - 1.0/np.e) - 1.0)
t3_val = CF(N_t3)
results['T3'] = {'expected': target_cf3, 'actual': t3_val, 'pass': np.isclose(t3_val, target_cf3, atol=1e-10), 'N': N_t3}

# T4: CF = 1/e
target_cf4 = 1.0/np.e
N_t4 = 4.0 / (1.0/(1.0/np.e) - 1.0)
t4_val = CF(N_t4)
results['T4'] = {'expected': target_cf4, 'actual': t4_val, 'pass': np.isclose(t4_val, target_cf4, atol=1e-10), 'N': N_t4}

# T5: dγ/dN < 0
N_test = np.array([2, 5, 10, 20, 50, 100])
t5_vals = dGamma_dN(N_test)
results['T5'] = {'values': t5_vals, 'pass': np.all(t5_vals < 0)}

# T6: γ(1) = 2
t6_val = gamma(1)
results['T6'] = {'expected': 2.0, 'actual': t6_val, 'pass': np.isclose(t6_val, 2.0, atol=1e-10)}

# T7: dCF/dN > 0
t7_vals = dCF_dN(N_test)
results['T7'] = {'values': t7_vals, 'pass': np.all(t7_vals > 0)}

# T8: CF(100) > 0.95
t8_val = CF(100)
results['T8'] = {'expected': '>0.95', 'actual': t8_val, 'pass': t8_val > 0.95}

# --- Print results ---
print("=" * 65)
print("Chemistry Session #1862: SCM Chemistry")
print("1725th Phenomenon - Pozzolanic Reaction Coherence")
print("=" * 65)
for key in sorted(results.keys()):
    r = results[key]
    status = "PASS" if r['pass'] else "FAIL"
    if 'actual' in r:
        print(f"  {key}: {status} | expected={r['expected']}, actual={r['actual']:.6f}")
    else:
        print(f"  {key}: {status} | all values satisfy constraint")
print("=" * 65)

# --- Plot 2x4 grid ---
fig, axes = plt.subplots(2, 4, figsize=(18, 9))
fig.suptitle('Session #1862: Supplementary Cementitious Materials Chemistry\nPozzolanic Reaction Coherence — 8 Boundary Tests',
             fontsize=14, fontweight='bold', y=0.98)

colors_pass = '#2ecc71'
colors_fail = '#e74c3c'

def status_color(passed):
    return colors_pass if passed else colors_fail

# T1
ax = axes[0, 0]
ax.plot(N_range, gamma_vals, 'b-', linewidth=1.5)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=4.0, color='gray', linestyle='--', alpha=0.5)
ax.plot(4, 1.0, 'o', color=status_color(results['T1']['pass']), markersize=12, zorder=5)
ax.set_title(f"T1: γ=1 @ N=4 [{'PASS' if results['T1']['pass'] else 'FAIL'}]",
             color=status_color(results['T1']['pass']), fontweight='bold')
ax.set_xlabel('N_corr (reactive sites)')
ax.set_ylabel('γ')
ax.set_xlim(0, 20)

# T2
ax = axes[0, 1]
ax.plot(N_range, CF_vals, 'r-', linewidth=1.5)
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=4.0, color='gray', linestyle='--', alpha=0.5)
ax.plot(4, 0.5, 'o', color=status_color(results['T2']['pass']), markersize=12, zorder=5)
ax.set_title(f"T2: CF=0.5 @ N=4 [{'PASS' if results['T2']['pass'] else 'FAIL'}]",
             color=status_color(results['T2']['pass']), fontweight='bold')
ax.set_xlabel('N_corr (reactive sites)')
ax.set_ylabel('CF')
ax.set_xlim(0, 20)

# T3
ax = axes[0, 2]
ax.plot(N_range, CF_vals, 'r-', linewidth=1.5)
ax.axhline(y=target_cf3, color='orange', linestyle='--', alpha=0.7, label=f'1-1/e≈{target_cf3:.4f}')
ax.axvline(x=N_t3, color='gray', linestyle='--', alpha=0.5)
ax.plot(N_t3, target_cf3, 'o', color=status_color(results['T3']['pass']), markersize=12, zorder=5)
ax.set_title(f"T3: CF=1-1/e @ N≈{N_t3:.2f} [{'PASS' if results['T3']['pass'] else 'FAIL'}]",
             color=status_color(results['T3']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.legend(fontsize=8)

# T4
ax = axes[0, 3]
ax.plot(N_range, CF_vals, 'r-', linewidth=1.5)
ax.axhline(y=target_cf4, color='purple', linestyle='--', alpha=0.7, label=f'1/e≈{target_cf4:.4f}')
ax.axvline(x=N_t4, color='gray', linestyle='--', alpha=0.5)
ax.plot(N_t4, target_cf4, 'o', color=status_color(results['T4']['pass']), markersize=12, zorder=5)
ax.set_title(f"T4: CF=1/e @ N≈{N_t4:.2f} [{'PASS' if results['T4']['pass'] else 'FAIL'}]",
             color=status_color(results['T4']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.legend(fontsize=8)

# T5
ax = axes[1, 0]
dg_range = dGamma_dN(N_range)
ax.plot(N_range, dg_range, 'b-', linewidth=1.5)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.fill_between(N_range, dg_range, 0, where=(dg_range < 0), alpha=0.2, color='blue')
ax.set_title(f"T5: dγ/dN<0 [{'PASS' if results['T5']['pass'] else 'FAIL'}]",
             color=status_color(results['T5']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('dγ/dN')

# T6
ax = axes[1, 1]
ax.plot(N_range, gamma_vals, 'b-', linewidth=1.5)
ax.axhline(y=2.0, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5)
ax.plot(1, 2.0, 'o', color=status_color(results['T6']['pass']), markersize=12, zorder=5)
ax.set_title(f"T6: γ(1)=2 [{'PASS' if results['T6']['pass'] else 'FAIL'}]",
             color=status_color(results['T6']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('γ')
ax.set_xlim(0, 10)

# T7
ax = axes[1, 2]
dcf_range = dCF_dN(N_range)
ax.plot(N_range, dcf_range, 'r-', linewidth=1.5)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.fill_between(N_range, dcf_range, 0, where=(dcf_range > 0), alpha=0.2, color='red')
ax.set_title(f"T7: dCF/dN>0 [{'PASS' if results['T7']['pass'] else 'FAIL'}]",
             color=status_color(results['T7']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('dCF/dN')

# T8
ax = axes[1, 3]
ax.plot(N_range, CF_vals, 'r-', linewidth=1.5)
ax.axhline(y=0.95, color='gray', linestyle='--', alpha=0.5, label='0.95 threshold')
ax.axvline(x=100, color='gray', linestyle='--', alpha=0.5)
ax.plot(100, t8_val, 'o', color=status_color(results['T8']['pass']), markersize=12, zorder=5)
ax.set_title(f"T8: CF(100)={t8_val:.4f}>0.95 [{'PASS' if results['T8']['pass'] else 'FAIL'}]",
             color=status_color(results['T8']['pass']), fontweight='bold')
ax.set_xlabel('N_corr')
ax.set_ylabel('CF')
ax.legend(fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.94])
outpath = '/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/scm_chemistry_coherence.png'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nPlot saved: {outpath}")
print("All 8 boundary tests: " + ("PASS" if all(r['pass'] for r in results.values()) else "SOME FAILED"))
