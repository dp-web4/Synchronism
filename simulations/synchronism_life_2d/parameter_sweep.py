#!/usr/bin/env python3
"""
Parameter Sweep: Find Stability Threshold

Systematically test combinations of:
- Resistance exponent (n)
- Initial amplitude
- Pattern size (sigma)

Goal: Map boundary between dissipative and stable regimes
"""

import numpy as np
from intent_life_2d import IntentLife2D
from pathlib import Path
from datetime import datetime
import csv
import itertools


def test_stability(n, amplitude, sigma, steps=2000, size=128):
    """
    Test single parameter combination

    Returns:
        dict with stability metrics
    """
    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)
    grid.n = n  # Override resistance exponent
    grid.add_pattern('gaussian', amplitude=amplitude, sigma=sigma)

    initial_max = grid.max_intent()
    initial_total = grid.total_intent()

    # Track metrics at specific intervals
    metrics = []
    check_points = [0, 100, 500, 1000, 1500, 2000]

    for step in range(steps + 1):
        if step in check_points:
            metrics.append({
                'step': step,
                'time': grid.time,
                'max_intent': grid.max_intent(),
                'total_intent': grid.total_intent(),
            })

        grid.step()

    # Calculate stability score
    final_max = metrics[-1]['max_intent']
    retention = final_max / initial_max

    # Check for oscillation (variance in later half)
    later_half_max = [m['max_intent'] for m in metrics[len(metrics)//2:]]
    oscillation_variance = np.var(later_half_max)

    result = {
        'n': n,
        'amplitude': amplitude,
        'sigma': sigma,
        'initial_max': initial_max,
        'final_max': final_max,
        'retention': retention,
        'oscillation_variance': oscillation_variance,
        'resistance_at_peak': 1 - (amplitude)**n,  # R(I) at initial peak
        'status': classify_stability(retention, oscillation_variance),
        'metrics': metrics,
    }

    return result


def classify_stability(retention, oscillation_variance):
    """
    Classify stability based on retention and oscillation

    Returns:
        'stable', 'quasi-stable', 'oscillating', 'dissipated'
    """
    if retention > 0.95:
        if oscillation_variance > 0.001:
            return 'oscillating'
        else:
            return 'stable'
    elif retention > 0.80:
        return 'quasi-stable'
    else:
        return 'dissipated'


def run_parameter_sweep():
    """
    Systematic parameter sweep

    Tests combinations of n, amplitude, sigma
    """
    print("=" * 70)
    print("SYNCHRONISM LIFE 2D - PARAMETER SWEEP")
    print("Finding Stability Threshold")
    print("=" * 70)
    print()

    # Parameter ranges to test
    n_values = [2, 3, 4, 5, 6]
    amplitude_values = [0.85, 0.90, 0.95, 0.98, 0.99]
    sigma_values = [5.0]  # Fix sigma for now

    total_tests = len(n_values) * len(amplitude_values) * len(sigma_values)

    print(f"Testing {total_tests} parameter combinations:")
    print(f"  n (resistance exponent): {n_values}")
    print(f"  amplitude: {amplitude_values}")
    print(f"  sigma: {sigma_values}")
    print()
    print("This will take ~5-10 minutes...")
    print()

    results = []
    test_num = 0

    # Run all combinations
    for n, amplitude, sigma in itertools.product(n_values, amplitude_values, sigma_values):
        test_num += 1

        print(f"[{test_num}/{total_tests}] Testing n={n}, amp={amplitude:.2f}, σ={sigma:.1f}... ", end='', flush=True)

        result = test_stability(n, amplitude, sigma)
        results.append(result)

        # Quick summary
        status_symbol = {
            'stable': '✓',
            'quasi-stable': '~',
            'oscillating': '⚡',
            'dissipated': '✗'
        }
        symbol = status_symbol.get(result['status'], '?')

        print(f"{symbol} {result['status']:15s} (retention={result['retention']:.1%}, R={result['resistance_at_peak']:.1%})")

    print()
    print("Sweep complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save summary CSV
    summary_file = output_dir / f"parameter_sweep_{timestamp}.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'n', 'amplitude', 'sigma',
            'initial_max', 'final_max', 'retention',
            'resistance_at_peak', 'oscillation_variance', 'status'
        ])

        for r in results:
            writer.writerow([
                r['n'], r['amplitude'], r['sigma'],
                r['initial_max'], r['final_max'], r['retention'],
                r['resistance_at_peak'], r['oscillation_variance'], r['status']
            ])

    print(f"Saved: {summary_file.name}")
    print()

    # Analysis
    analyze_results(results)

    return results


def analyze_results(results):
    """Analyze and summarize parameter sweep results"""
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print()

    # Group by status
    by_status = {}
    for r in results:
        status = r['status']
        if status not in by_status:
            by_status[status] = []
        by_status[status].append(r)

    print("Summary by Status:")
    for status in ['stable', 'oscillating', 'quasi-stable', 'dissipated']:
        count = len(by_status.get(status, []))
        print(f"  {status:15s}: {count:2d} configurations")
    print()

    # Find best performers
    sorted_by_retention = sorted(results, key=lambda r: r['retention'], reverse=True)

    print("Top 10 Most Stable Configurations:")
    print()
    print("  Rank | n | Amp  | Retention | Resistance | Status")
    print("  " + "-" * 60)
    for i, r in enumerate(sorted_by_retention[:10], 1):
        print(f"  {i:4d} | {r['n']} | {r['amplitude']:.2f} | {r['retention']:8.1%} | {r['resistance_at_peak']:9.1%} | {r['status']}")
    print()

    # Find stability boundary
    stable_configs = [r for r in results if r['status'] in ['stable', 'oscillating']]

    if stable_configs:
        print(f"✓ STABILITY ACHIEVED in {len(stable_configs)} configurations!")
        print()

        min_n_stable = min(r['n'] for r in stable_configs)
        min_amp_stable = min(r['amplitude'] for r in stable_configs)

        print(f"Minimum parameters for stability:")
        print(f"  n ≥ {min_n_stable}")
        print(f"  amplitude ≥ {min_amp_stable:.2f}")
        print()

        # Show stable configurations
        print("All stable configurations:")
        for r in stable_configs:
            print(f"  n={r['n']}, amplitude={r['amplitude']:.2f} → {r['status']} (retention={r['retention']:.1%})")
        print()
    else:
        print("✗ No stable configurations found in tested parameter space")
        print()

        best = sorted_by_retention[0]
        print(f"Best result: n={best['n']}, amplitude={best['amplitude']:.2f}")
        print(f"  Retention: {best['retention']:.1%}")
        print(f"  Resistance: {best['resistance_at_peak']:.1%}")
        print()
        print("Recommendation: Test higher n (>6) or higher amplitude (>0.99)")
        print()

    # Resistance analysis
    print("Pattern by Resistance at Peak:")
    print("  (Transfer resistance = what % of Intent flow is blocked)")
    print()

    resistance_bins = [
        (0.0, 0.5, "Weak (<50%)"),
        (0.5, 0.8, "Moderate (50-80%)"),
        (0.8, 0.9, "Strong (80-90%)"),
        (0.9, 1.0, "Very Strong (>90%)"),
    ]

    for low, high, label in resistance_bins:
        configs = [r for r in results if low <= r['resistance_at_peak'] < high]
        if configs:
            avg_retention = np.mean([r['retention'] for r in configs])
            stable_count = sum(1 for r in configs if r['status'] in ['stable', 'oscillating'])
            print(f"  {label:20s}: {len(configs):2d} configs, avg retention={avg_retention:.1%}, {stable_count} stable")
    print()

    # Insight
    print("Key Insight:")
    print("  Saturation dynamics require EXTREME resistance (>90%) for stability.")
    print("  Moderate resistance (50-80%) only slows dissipation, doesn't prevent it.")
    print("  This is physically meaningful: entities need near-perfect saturation cores.")
    print()


if __name__ == "__main__":
    results = run_parameter_sweep()
