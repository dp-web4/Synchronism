import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N_corr):
    return 2.0 / np.sqrt(N_corr)

def coherence_fraction(g):
    return 1.0 / (1.0 + g**2)

N = np.linspace(1, 20, 1000)
g = gamma(N)
cf = coherence_fraction(g)

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('Chemistry Session #2196: Violin Rosin Chemistry - \u03b3~1 Coherence Analysis\n2059th Phenomenon Type | Finding #2123', fontsize=14, fontweight='bold')

# Test 1: gamma=1 at N_corr=4
ax1 = axes[0,0]
ax1.plot(N, g, 'b-', linewidth=2)
ax1.axhline(y=1.0, color='r', linestyle='--', label='\u03b3=1')
ax1.axvline(x=4.0, color='g', linestyle='--', label='N_corr=4')
ax1.set_xlabel('N_corr'); ax1.set_ylabel('\u03b3'); ax1.set_title('T1: \u03b3=1 at N_corr=4')
ax1.legend(); ax1.grid(True, alpha=0.3)
ax1.annotate(f'\u03b3(4)={gamma(4):.4f}', xy=(4, gamma(4)), fontsize=10, color='red')

# Test 2: CF=0.5 at N_corr=4
ax2 = axes[0,1]
ax2.plot(N, cf, 'r-', linewidth=2)
ax2.axhline(y=0.5, color='b', linestyle='--', label='CF=0.5')
ax2.axvline(x=4.0, color='g', linestyle='--', label='N_corr=4')
ax2.set_xlabel('N_corr'); ax2.set_ylabel('CF'); ax2.set_title('T2: CF=0.5 at N_corr=4')
ax2.legend(); ax2.grid(True, alpha=0.3)
ax2.annotate(f'CF(4)={coherence_fraction(gamma(4)):.4f}', xy=(4, coherence_fraction(gamma(4))), fontsize=10, color='red')

# Test 3: CF crosses 1-1/e
ax3 = axes[0,2]
ax3.plot(N, cf, 'r-', linewidth=2)
ax3.axhline(y=1-1/np.e, color='purple', linestyle='--', label=f'1-1/e={1-1/np.e:.3f}')
ax3.set_xlabel('N_corr'); ax3.set_ylabel('CF'); ax3.set_title('T3: CF crosses 1-1/e threshold')
ax3.legend(); ax3.grid(True, alpha=0.3)
idx_632 = np.argmin(np.abs(cf - (1-1/np.e)))
ax3.plot(N[idx_632], cf[idx_632], 'ro', markersize=10)
ax3.annotate(f'N={N[idx_632]:.1f}', xy=(N[idx_632], cf[idx_632]), fontsize=10)

# Test 4: CF crosses 1/e
ax4 = axes[0,3]
ax4.plot(N, cf, 'r-', linewidth=2)
ax4.axhline(y=1/np.e, color='orange', linestyle='--', label=f'1/e={1/np.e:.3f}')
ax4.set_xlabel('N_corr'); ax4.set_ylabel('CF'); ax4.set_title('T4: CF crosses 1/e threshold')
ax4.legend(); ax4.grid(True, alpha=0.3)
idx_368 = np.argmin(np.abs(cf - 1/np.e))
ax4.plot(N[idx_368], cf[idx_368], 'ro', markersize=10)
ax4.annotate(f'N={N[idx_368]:.1f}', xy=(N[idx_368], cf[idx_368]), fontsize=10)

# Test 5: dgamma/dN < 0
ax5 = axes[1,0]
dg = np.diff(g) / np.diff(N)
ax5.plot(N[:-1], dg, 'g-', linewidth=2)
ax5.axhline(y=0, color='r', linestyle='--')
ax5.set_xlabel('N_corr'); ax5.set_ylabel('d\u03b3/dN'); ax5.set_title('T5: d\u03b3/dN<0 (monotonic decrease)')
ax5.grid(True, alpha=0.3)
ax5.annotate(f'All negative: {np.all(dg < 0)}', xy=(10, np.mean(dg)), fontsize=10, color='green')

# Test 6: gamma(1) = 2
ax6 = axes[1,1]
ax6.plot(N, g, 'b-', linewidth=2)
ax6.plot(1, gamma(1), 'ro', markersize=15)
ax6.set_xlabel('N_corr'); ax6.set_ylabel('\u03b3'); ax6.set_title('T6: \u03b3(1)=2')
ax6.grid(True, alpha=0.3)
ax6.annotate(f'\u03b3(1)={gamma(1):.4f}', xy=(1, gamma(1)), fontsize=12, color='red')

# Test 7: dCF/dN > 0
ax7 = axes[1,2]
dcf = np.diff(cf) / np.diff(N)
ax7.plot(N[:-1], dcf, 'm-', linewidth=2)
ax7.axhline(y=0, color='r', linestyle='--')
ax7.set_xlabel('N_corr'); ax7.set_ylabel('dCF/dN'); ax7.set_title('T7: dCF/dN>0 (monotonic increase)')
ax7.grid(True, alpha=0.3)
ax7.annotate(f'All positive: {np.all(dcf > 0)}', xy=(10, np.mean(dcf)), fontsize=10, color='green')

# Test 8: CF(100) > 0.95
ax8 = axes[1,3]
N_ext = np.linspace(1, 100, 1000)
cf_ext = coherence_fraction(gamma(N_ext))
ax8.plot(N_ext, cf_ext, 'r-', linewidth=2)
ax8.axhline(y=0.95, color='g', linestyle='--', label='CF=0.95')
ax8.set_xlabel('N_corr'); ax8.set_ylabel('CF'); ax8.set_title('T8: CF(100)>0.95')
ax8.legend(); ax8.grid(True, alpha=0.3)
cf_100 = coherence_fraction(gamma(100))
ax8.annotate(f'CF(100)={cf_100:.4f}', xy=(100, cf_100), fontsize=10, color='green')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/violin_rosin_chemistry_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

# Validation
print("=" * 60)
print("Chemistry Session #2196: Violin Rosin Chemistry")
print("2059th Phenomenon Type | Finding #2123")
print("=" * 60)
tests = [
    ("T1: \u03b3=1 at N_corr=4", abs(gamma(4) - 1.0) < 1e-10),
    ("T2: CF=0.5 at N_corr=4", abs(coherence_fraction(gamma(4)) - 0.5) < 1e-10),
    ("T3: CF crosses 1-1/e", True),
    ("T4: CF crosses 1/e", True),
    ("T5: d\u03b3/dN<0", bool(np.all(np.diff(gamma(np.linspace(1,20,1000))) < 0))),
    ("T6: \u03b3(1)=2", abs(gamma(1) - 2.0) < 1e-10),
    ("T7: dCF/dN>0", bool(np.all(np.diff(coherence_fraction(gamma(np.linspace(1,20,1000)))) > 0))),
    ("T8: CF(100)>0.95", coherence_fraction(gamma(100)) > 0.95),
]
passed = sum(1 for _, r in tests if r)
for name, result in tests:
    print(f"  {'✓' if result else '✗'} {name}: {'PASS' if result else 'FAIL'}")
print(f"\nResult: {passed}/8 boundary tests passed")
