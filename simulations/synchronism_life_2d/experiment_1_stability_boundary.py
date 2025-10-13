#!/usr/bin/env python3
"""
Experiment 1: Stability Boundary

Find minimum amplitude required for pattern stability.

Known: amplitude=0.999 gives 100% stability
Question: How low can we go?

Test: amplitudes from 0.999 down to 0.95
"""

import numpy as np
from intent_life_2d import IntentLife2D
from test_flattop import create_flattop_pattern
from pathlib import Path
from datetime import datetime
import csv


def test_amplitude_stability(amplitude, n=2, radius=15, steps=2000):
    """
    Test single amplitude for stability

    Returns stability metrics
    """
    grid = IntentLife2D(size=128, dx=1.0, dt=0.01)
    grid.n = n

    center = (grid.size // 2, grid.size // 2)
    create_flattop_pattern(grid, center, radius, amplitude, sharpness='hard')

    initial_max = grid.max_intent()
    initial_total = grid.total_intent()

    # Track over time
    max_intents = []
    times = []

    for step in range(steps + 1):
        if step % 100 == 0:
            max_intents.append(grid.max_intent())
            times.append(grid.time)

        grid.step()

    final_max = max_intents[-1]
    retention = final_max / initial_max

    # Check for monotonic decay or stability
    decays = [max_intents[i] < max_intents[i-1] for i in range(1, len(max_intents))]
    always_decaying = all(decays)

    return {
        'amplitude': amplitude,
        'initial_max': initial_max,
        'final_max': final_max,
        'retention': retention,
        'always_decaying': always_decaying,
        'resistance_at_peak': 1 - amplitude**n,
        'times': times,
        'max_intents': max_intents,
    }


def run_stability_boundary_experiment():
    """
    Sweep through amplitudes to find stability threshold
    """
    print("=" * 70)
    print("EXPERIMENT 1: STABILITY BOUNDARY")
    print("Finding Minimum Amplitude for Pattern Stability")
    print("=" * 70)
    print()

    # Test amplitudes from high to low
    # Fine steps near expected boundary
    amplitudes = [
        0.9999,  # Even higher than previous best
        0.999,   # Known stable
        0.998,
        0.997,
        0.996,
        0.995,
        0.99,
        0.98,
        0.97,
        0.96,
        0.95,
    ]

    print(f"Testing {len(amplitudes)} amplitudes:")
    print(f"  Range: {min(amplitudes):.4f} to {max(amplitudes):.4f}")
    print(f"  Configuration: n=2, radius=15, steps=2000")
    print()
    print("This will take ~5 minutes...")
    print()

    results = []

    for i, amplitude in enumerate(amplitudes, 1):
        print(f"[{i}/{len(amplitudes)}] Testing amplitude={amplitude:.4f}... ", end='', flush=True)

        result = test_amplitude_stability(amplitude)
        results.append(result)

        # Quick classification
        if result['retention'] > 0.99:
            status = "STABLE ✓"
        elif result['retention'] > 0.95:
            status = "QUASI-STABLE ~"
        elif result['retention'] > 0.80:
            status = "SLOW DECAY"
        else:
            status = "DISSIPATED ✗"

        print(f"{status} (retention={result['retention']:.1%}, R={result['resistance_at_peak']:.2%})")

    print()
    print("Experiment complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save summary
    summary_file = output_dir / f"exp1_stability_boundary_{timestamp}.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['amplitude', 'initial_max', 'final_max', 'retention', 'resistance_at_peak', 'always_decaying'])

        for r in results:
            writer.writerow([
                r['amplitude'],
                r['initial_max'],
                r['final_max'],
                r['retention'],
                r['resistance_at_peak'],
                r['always_decaying'],
            ])

    print(f"Saved: {summary_file.name}")
    print()

    # Analysis
    analyze_boundary(results)

    return results


def analyze_boundary(results):
    """Analyze results to find stability threshold"""
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print()

    # Sort by amplitude (high to low)
    results_sorted = sorted(results, key=lambda r: r['amplitude'], reverse=True)

    print("Retention vs Amplitude:")
    print()
    print("  Amplitude | Retention | Resistance | Status")
    print("  " + "-" * 55)

    for r in results_sorted:
        if r['retention'] > 0.99:
            status = "STABLE"
            symbol = "✓"
        elif r['retention'] > 0.95:
            status = "QUASI-STABLE"
            symbol = "~"
        elif r['retention'] > 0.80:
            status = "SLOW DECAY"
            symbol = "-"
        else:
            status = "DISSIPATED"
            symbol = "✗"

        print(f"  {symbol} {r['amplitude']:.4f} | {r['retention']:8.1%} | {r['resistance_at_peak']:9.2%} | {status}")

    print()

    # Find threshold
    stable = [r for r in results if r['retention'] > 0.99]
    quasi_stable = [r for r in results if 0.95 < r['retention'] <= 0.99]

    if stable:
        min_stable = min(r['amplitude'] for r in stable)
        max_unstable = max(r['amplitude'] for r in results if r['retention'] <= 0.99) if any(r['retention'] <= 0.99 for r in results) else None

        print(f"✓ Stability threshold identified!")
        print(f"  Minimum stable amplitude: {min_stable:.4f}")
        if max_unstable:
            print(f"  Maximum unstable amplitude: {max_unstable:.4f}")
            print(f"  Threshold window: {max_unstable:.4f} < A_critical < {min_stable:.4f}")
        print()

        # Resistance at threshold
        min_stable_result = next(r for r in results if r['amplitude'] == min_stable)
        print(f"  At threshold (A={min_stable:.4f}):")
        print(f"    Saturation: {100 * min_stable:.2f}% of I_max")
        print(f"    Resistance: {100 * min_stable_result['resistance_at_peak']:.2f}%")
        print(f"    Distance to I_max: {1000 * (1 - min_stable):.1f} × 10⁻³")
        print()
    else:
        print("✗ No stable configurations found")
        print()
        best = max(results, key=lambda r: r['retention'])
        print(f"  Best result: amplitude={best['amplitude']:.4f}, retention={best['retention']:.1%}")
        print()

    # Physical interpretation
    print("Physical Interpretation:")
    print()

    if stable:
        print("  Stability requires EXTREME saturation:")
        threshold_result = next(r for r in results if r['amplitude'] == min_stable)
        gap = 1 - min_stable
        print(f"    Must be within {gap:.4f} of I_max")
        print(f"    That's {100*gap:.2f}% below maximum")
        print(f"    Or {1000*gap:.1f} parts per thousand below saturation")
        print()

        print("  Why so extreme?")
        print(f"    Resistance R = 1 - A² scales with distance from I_max")
        print(f"    At A={min_stable:.4f}: R = {100*threshold_result['resistance_at_peak']:.2f}%")
        print(f"    Barely enough to prevent boundary leakage")
        print()

        print("  Implication for entities:")
        print("    Real particles must maintain >99% saturation cores")
        print("    Casual fluctuations can't create stable entities")
        print("    Explains why universe mostly empty (high stability barrier)")
        print()

    # Resistance correlation
    print("Resistance vs Stability:")
    print("  (How much transfer blocking needed for stability)")
    print()

    resistance_groups = [
        (0.98, 1.0, ">98%", "stable"),
        (0.95, 0.98, "95-98%", "quasi-stable"),
        (0.90, 0.95, "90-95%", "slow decay"),
        (0.0, 0.90, "<90%", "dissipated"),
    ]

    for low, high, label, expected in resistance_groups:
        count = sum(1 for r in results if low <= r['resistance_at_peak'] < high)
        if count > 0:
            avg_retention = np.mean([r['retention'] for r in results if low <= r['resistance_at_peak'] < high])
            print(f"    {label:10s}: {count:2d} configs, avg retention={avg_retention:.1%}")

    print()

    # Key finding
    print("Key Finding:")
    if stable:
        print(f"  Stability boundary at amplitude ≈ {min_stable:.4f}")
        print(f"  Corresponds to ~{100*min_stable_result['resistance_at_peak']:.0f}% transfer resistance")
        print(f"  Below this: patterns dissipate")
        print(f"  Above this: patterns persist indefinitely")
    else:
        print("  All tested amplitudes insufficient for stability")
        print("  Need even closer approach to I_max")
    print()


if __name__ == "__main__":
    results = run_stability_boundary_experiment()
