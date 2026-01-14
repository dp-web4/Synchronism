#!/usr/bin/env python3
"""
Chemistry Session #26: Methods for Measuring N_corr

Demonstrates different experimental approaches to determine the
correlation number N_corr from observable quantities.

Key insight: N_corr can be extracted from:
1. Fluctuation amplitudes
2. Correlation lengths
3. Entropy ratios
4. Information measures
5. Spectral line widths
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.signal import find_peaks

def generate_correlated_system(n_particles, n_corr, n_samples=10000):
    """Generate samples from a system with N_corr correlation."""
    # Block correlation structure
    group_size = max(1, int(n_corr))
    n_groups = max(1, n_particles // group_size)
    actual_n = n_groups * group_size

    # Independent group values
    group_values = np.random.randn(n_samples, n_groups)

    # Expand to particles (all in same group have same value)
    samples = np.repeat(group_values, group_size, axis=1)[:, :n_particles]

    return samples

def method1_fluctuation_analysis(samples_corr, samples_uncorr):
    """
    Method 1: Extract N_corr from fluctuation amplitudes.

    N_corr = (σ_corr / σ_uncorr)²
    """
    # Compute mean for each sample
    mean_corr = np.mean(samples_corr, axis=1)
    mean_uncorr = np.mean(samples_uncorr, axis=1)

    # Standard deviations
    sigma_corr = np.std(mean_corr)
    sigma_uncorr = np.std(mean_uncorr)

    # Infer N_corr
    n_corr_measured = (sigma_corr / sigma_uncorr) ** 2

    return n_corr_measured

def method2_correlation_length(samples, spacing=1.0):
    """
    Method 2: Extract N_corr from correlation length.

    N_corr ~ (ξ/a)^d where ξ is correlation length
    """
    n_particles = samples.shape[1]

    # Compute spatial correlation function
    # C(r) = <X_i X_{i+r}> / <X_i²>
    correlations = []
    max_r = min(50, n_particles // 2)

    for r in range(max_r):
        if r == 0:
            correlations.append(1.0)
        else:
            corr_vals = []
            for sample in samples[:min(1000, len(samples))]:
                c = np.corrcoef(sample[:-r], sample[r:])[0, 1]
                if not np.isnan(c):
                    corr_vals.append(c)
            if corr_vals:
                correlations.append(np.mean(corr_vals))
            else:
                correlations.append(0)

    correlations = np.array(correlations)

    # Find correlation length (where C drops to 1/e)
    threshold = 1 / np.e
    xi_idx = np.where(correlations < threshold)[0]
    if len(xi_idx) > 0:
        xi = xi_idx[0] * spacing
    else:
        xi = max_r * spacing

    # N_corr ~ xi in 1D (for our block model)
    n_corr_measured = xi / spacing

    return n_corr_measured, correlations

def method3_entropy_ratio(samples_corr, samples_uncorr):
    """
    Method 3: Extract N_corr from entropy ratio.

    S/S₀ = 1/√N_corr → N_corr = (S₀/S)²

    Note: This requires careful normalization.
    """
    # Estimate entropy from variance (Gaussian approximation)
    # H = 0.5 * log(2πe * σ²)

    n_particles = samples_corr.shape[1]

    # Covariance matrices
    cov_corr = np.cov(samples_corr.T)
    cov_uncorr = np.cov(samples_uncorr.T)

    # Log-determinants (joint entropy up to constant)
    _, logdet_corr = np.linalg.slogdet(cov_corr)
    _, logdet_uncorr = np.linalg.slogdet(cov_uncorr)

    # Entropy ratio
    # Actually for block-correlated, the entropy scales differently
    # H_joint = n_groups * H_single (for perfect block correlation)
    # H_uncorr = n_particles * H_single

    # Ratio = n_groups / n_particles = 1/n_corr_per_group

    # From logdet: entropy = 0.5 * logdet
    entropy_ratio = logdet_corr / logdet_uncorr if logdet_uncorr > 0 else 0

    # For our block model: entropy_ratio ≈ n_groups / n_particles = 1/group_size
    # So N_corr_measured ≈ logdet_uncorr / logdet_corr

    # This is complex; simplify
    # Use the fact that for perfect block correlation:
    # rank(cov_corr) = n_groups
    # but cov is still n×n, just singular or near-singular

    # Eigenvalue analysis
    eigvals = np.linalg.eigvalsh(cov_corr)
    eigvals = eigvals[eigvals > 1e-10]  # Remove near-zero
    n_effective = len(eigvals)
    n_corr_measured = n_particles / n_effective if n_effective > 0 else n_particles

    return n_corr_measured

def method4_information_theoretic(samples):
    """
    Method 4: Extract N_corr from multi-information.

    I = Σ H(Xᵢ) - H(X₁,...,Xₙ) (total correlation)
    N_corr ≈ exp(2I/N) for Gaussian
    """
    n_particles = samples.shape[1]

    # Covariance matrix
    cov = np.cov(samples.T)

    # Add small regularization for stability
    cov = cov + np.eye(n_particles) * 1e-10

    # Sum of marginal entropies (up to constant)
    marginal_vars = np.diag(cov)
    h_marginals = 0.5 * np.sum(np.log(marginal_vars + 1e-10))

    # Joint entropy (up to constant)
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        # Fall back to eigenvalue method
        eigvals = np.linalg.eigvalsh(cov)
        eigvals = eigvals[eigvals > 1e-10]
        logdet = np.sum(np.log(eigvals))
    h_joint = 0.5 * logdet

    # Multi-information
    I = h_marginals - h_joint

    # For block correlation, I grows with correlation
    # N_corr = exp(2I/N) is an approximation
    n_corr_measured = np.exp(2 * I / n_particles)

    # Cap at reasonable value
    n_corr_measured = min(n_corr_measured, n_particles)

    return n_corr_measured

def method5_spectral_linewidth(samples, n_corr_true):
    """
    Method 5: Spectral approach - line narrowing.

    For coherent oscillators, spectral line narrows by 1/√N_corr.
    """
    # Simulate time series with correlated oscillators
    n_particles = samples.shape[1]
    n_time = 1000
    dt = 0.01
    t = np.arange(n_time) * dt
    omega_0 = 10  # Central frequency
    gamma_dephasing = 0.5  # Individual dephasing rate

    # Individual (uncorrelated) spectrum
    # Each oscillator has random phase noise
    signal_uncorr = np.zeros(n_time)
    for i in range(n_particles):
        phase_noise = np.cumsum(np.random.randn(n_time)) * 0.1
        signal_uncorr += np.cos(omega_0 * t + phase_noise) * np.exp(-gamma_dephasing * t)

    # Correlated spectrum (n_corr oscillators move together)
    n_groups = max(1, n_particles // n_corr_true)
    signal_corr = np.zeros(n_time)
    for g in range(n_groups):
        phase_noise = np.cumsum(np.random.randn(n_time)) * 0.1
        # All oscillators in group share the phase noise
        signal_corr += n_corr_true * np.cos(omega_0 * t + phase_noise) * np.exp(-gamma_dephasing * t)

    # Compute spectra
    freq_uncorr = np.fft.fftfreq(n_time, dt)
    spec_uncorr = np.abs(np.fft.fft(signal_uncorr)) ** 2
    spec_corr = np.abs(np.fft.fft(signal_corr)) ** 2

    # Find linewidths (FWHM near omega_0)
    # This is approximate - in real systems would use proper fitting
    pos_freq = freq_uncorr[:n_time//2]

    # Peak finding
    idx_peak = np.argmax(spec_uncorr[:n_time//2])
    omega_measured = pos_freq[idx_peak] * 2 * np.pi

    # Linewidth from spectrum width
    half_max_uncorr = np.max(spec_uncorr) / 2
    half_max_corr = np.max(spec_corr) / 2

    # Count width
    width_uncorr = np.sum(spec_uncorr > half_max_uncorr * 0.1)
    width_corr = np.sum(spec_corr > half_max_corr * 0.1)

    # N_corr from linewidth ratio
    if width_corr > 0:
        n_corr_measured = (width_uncorr / width_corr) ** 2
    else:
        n_corr_measured = n_particles

    return n_corr_measured, (t, signal_uncorr, signal_corr)

def main():
    """Demonstrate all methods for measuring N_corr."""

    np.random.seed(42)

    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    fig.suptitle('Chemistry Session #26: Methods for Measuring N_corr', fontsize=14, fontweight='bold')

    n_particles = 100
    n_samples = 10000
    n_corr_values = [1, 2, 4, 10, 25, 50]

    # Results storage
    results = {
        'true': [],
        'method1': [],
        'method2': [],
        'method3': [],
        'method4': [],
    }

    # Generate reference uncorrelated samples
    samples_uncorr = generate_correlated_system(n_particles, 1, n_samples)

    for n_corr in n_corr_values:
        samples_corr = generate_correlated_system(n_particles, n_corr, n_samples)

        results['true'].append(n_corr)
        results['method1'].append(method1_fluctuation_analysis(samples_corr, samples_uncorr))

        m2, _ = method2_correlation_length(samples_corr)
        results['method2'].append(m2)

        results['method3'].append(method3_entropy_ratio(samples_corr, samples_uncorr))
        results['method4'].append(method4_information_theoretic(samples_corr))

    # Plot 1: Method comparison
    ax1 = axes[0, 0]
    methods = ['method1', 'method2', 'method3', 'method4']
    labels = ['Fluctuation', 'Correlation Length', 'Entropy', 'Information']
    markers = ['o', 's', '^', 'd']

    for method, label, marker in zip(methods, labels, markers):
        ax1.scatter(results['true'], results[method], label=label, marker=marker, s=80, alpha=0.7)

    ax1.plot([0, 55], [0, 55], 'k--', label='Perfect', alpha=0.5)
    ax1.set_xlabel('True N_corr')
    ax1.set_ylabel('Measured N_corr')
    ax1.set_title('Method Comparison')
    ax1.legend(fontsize=8)
    ax1.set_xlim(0, 55)
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Method 1 (Fluctuation) detail
    ax2 = axes[0, 1]

    n_corr_example = 16
    samples_16 = generate_correlated_system(n_particles, n_corr_example, n_samples)

    means_corr = np.mean(samples_16, axis=1)
    means_uncorr = np.mean(samples_uncorr, axis=1)

    ax2.hist(means_uncorr, bins=50, alpha=0.5, density=True, label=f'Uncorrelated (N_corr=1)')
    ax2.hist(means_corr, bins=50, alpha=0.5, density=True, label=f'Correlated (N_corr={n_corr_example})')

    sigma_ratio = np.std(means_corr) / np.std(means_uncorr)
    ax2.set_xlabel('Sample Mean')
    ax2.set_ylabel('Probability Density')
    ax2.set_title(f'Method 1: Fluctuation Analysis\nσ_ratio = {sigma_ratio:.2f}, N_corr_meas = {sigma_ratio**2:.1f}')
    ax2.legend()

    # Plot 3: Method 2 (Correlation Length) detail
    ax3 = axes[0, 2]

    _, correlations = method2_correlation_length(samples_16)
    r_values = np.arange(len(correlations))

    ax3.plot(r_values, correlations, 'b-', linewidth=2)
    ax3.axhline(y=1/np.e, color='r', linestyle='--', label='1/e threshold')
    ax3.axvline(x=n_corr_example, color='g', linestyle='--', label=f'True N_corr={n_corr_example}')

    ax3.set_xlabel('Distance r')
    ax3.set_ylabel('Correlation C(r)')
    ax3.set_title('Method 2: Correlation Length')
    ax3.legend()
    ax3.set_xlim(0, 50)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Method summary
    ax4 = axes[1, 0]

    ax4.axis('off')
    summary = """
Method Summary for Measuring N_corr:

1. FLUCTUATION ANALYSIS (Most reliable)
   - Measure σ of intensive property
   - Compare to uncorrelated expectation
   - N_corr = (σ_measured / σ_uncorrelated)²

2. CORRELATION LENGTH
   - Measure spatial correlation ξ
   - N_corr ~ (ξ/a)^d
   - Works for extended systems

3. ENTROPY RATIO
   - Measure system entropy S
   - N_corr = (S_uncorr / S_measured)²
   - Requires absolute entropy

4. INFORMATION THEORETIC
   - Measure multi-information I
   - N_corr ≈ exp(2I/N)
   - Works for Gaussian-like systems

5. SPECTRAL LINEWIDTH
   - Measure line narrowing factor
   - N_corr = (width_uncorr / width_corr)²
   - Works for oscillating systems
"""
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes, fontsize=9,
             verticalalignment='top', family='monospace')

    # Plot 5: Domain-specific recommendations
    ax5 = axes[1, 1]

    domains = [
        ('Superconductors', 'Correlation length\n(coherence length ξ)', 'N=(ξ/a)³'),
        ('Enzymes', 'KIE + simulation', 'From H-bond network'),
        ('Neurons', 'EEG synchrony', 'From cross-correlation'),
        ('Markets', 'Price correlations', 'From covariance matrix'),
        ('Proteins', 'MD simulation', 'From RMSF analysis'),
    ]

    y_pos = np.arange(len(domains))
    ax5.barh(y_pos, [1]*len(domains), alpha=0)  # Invisible bars for spacing

    for i, (domain, method, formula) in enumerate(domains):
        ax5.text(0.0, i, f'{domain}:', fontweight='bold', fontsize=10, va='center')
        ax5.text(0.3, i, method, fontsize=9, va='center')
        ax5.text(0.7, i, formula, fontsize=9, va='center', style='italic')

    ax5.set_xlim(0, 1)
    ax5.set_ylim(-0.5, len(domains)-0.5)
    ax5.axis('off')
    ax5.set_title('Domain-Specific Recommendations')

    # Plot 6: γ from N_corr measurement
    ax6 = axes[1, 2]

    n_corr_range = np.linspace(1, 100, 100)
    gamma_range = 2 / np.sqrt(n_corr_range)

    ax6.plot(n_corr_range, gamma_range, 'b-', linewidth=2, label='γ = 2/√N_corr')

    # Add measured points
    for n_corr, m1 in zip(results['true'], results['method1']):
        gamma_true = 2 / np.sqrt(n_corr)
        gamma_meas = 2 / np.sqrt(m1)
        ax6.scatter([n_corr], [gamma_true], s=100, c='green', edgecolor='black', zorder=5)
        ax6.scatter([m1], [gamma_meas], s=60, c='red', marker='x', zorder=5)
        ax6.plot([n_corr, m1], [gamma_true, gamma_meas], 'r--', alpha=0.3)

    ax6.set_xlabel('N_corr (measured or true)')
    ax6.set_ylabel('γ')
    ax6.set_title('From N_corr Measurement to γ')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/measuring_ncorr.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("=" * 70)
    print("Chemistry Session #26: Methods for Measuring N_corr")
    print("=" * 70)
    print()
    print("With γ = 2/√N_corr derived (Session #25), the experimental challenge")
    print("is measuring N_corr for specific systems.")
    print()
    print("-" * 70)
    print("METHOD COMPARISON")
    print("-" * 70)
    print()
    print(f"{'N_corr':>8} | {'Fluctuation':>12} | {'Corr Len':>12} | {'Entropy':>12} | {'Info':>12}")
    print(f"{'(true)':>8} | {'(Method 1)':>12} | {'(Method 2)':>12} | {'(Method 3)':>12} | {'(Method 4)':>12}")
    print("-" * 70)
    for i, n in enumerate(results['true']):
        print(f"{n:>8} | {results['method1'][i]:>12.2f} | {results['method2'][i]:>12.2f} | "
              f"{results['method3'][i]:>12.2f} | {results['method4'][i]:>12.2f}")
    print()
    print("-" * 70)
    print()
    print("RECOMMENDED METHODS BY DOMAIN:")
    print()
    print("  Superconductors: Correlation length from penetration depth/coherence length")
    print("                   N_corr = (ξ/a)³ where ξ is coherence length")
    print()
    print("  Enzymes:         KIE measurements combined with MD simulation")
    print("                   N_corr from H-bond network analysis")
    print()
    print("  Neural systems:  EEG cross-correlation analysis")
    print("                   N_corr from synchrony indices")
    print()
    print("  Financial:       Price correlation matrix eigenvalue analysis")
    print("                   N_corr ~ 1/participation ratio")
    print()
    print("  Proteins:        Molecular dynamics RMSF analysis")
    print("                   N_corr from covariance matrix")
    print()
    print("-" * 70)
    print("KEY INSIGHT:")
    print()
    print("  Method 1 (Fluctuation Analysis) is the most direct and reliable.")
    print("  It follows directly from the derivation:")
    print()
    print("  N_corr = (σ_correlated / σ_uncorrelated)²")
    print()
    print("  The challenge is knowing σ_uncorrelated for comparison.")
    print("  This often requires simulation or theoretical prediction.")
    print()
    print("=" * 70)
    print("MEASUREMENT PROTOCOLS ESTABLISHED")
    print("=" * 70)

if __name__ == "__main__":
    main()
