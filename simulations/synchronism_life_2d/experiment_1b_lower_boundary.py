#!/usr/bin/env python3
"""
Experiment 1b: Lower Stability Boundary

First test showed stability down to 0.95!
Now test even lower amplitudes to find actual threshold.
"""

import numpy as np
from experiment_1_stability_boundary import test_amplitude_stability
from pathlib import Path
from datetime import datetime
import csv


def run_lower_boundary():
    """Test lower amplitudes"""
    print("=" * 70)
    print("EXPERIMENT 1b: LOWER STABILITY BOUNDARY")
    print("=" * 70)
    print()

    # Test lower amplitudes
    amplitudes = [
        0.95,   # Known stable
        0.94,
        0.93,
        0.92,
        0.91,
        0.90,
        0.85,
        0.80,
        0.75,
        0.70,
    ]

    print(f"Testing {len(amplitudes)} lower amplitudes")
    print()

    results = []

    for i, amplitude in enumerate(amplitudes, 1):
        print(f"[{i}/{len(amplitudes)}] Testing amplitude={amplitude:.2f}... ", end='', flush=True)

        result = test_amplitude_stability(amplitude)
        results.append(result)

        if result['retention'] > 0.99:
            status = "STABLE ✓"
        elif result['retention'] > 0.95:
            status = "QUASI ~"
        elif result['retention'] > 0.80:
            status = "SLOW"
        else:
            status = "DECAY ✗"

        print(f"{status} (ret={result['retention']:.1%}, R={result['resistance_at_peak']:.1%})")

    print()

    # Save
    output_dir = Path(__file__).parent / "output"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    csv_file = output_dir / f"exp1b_lower_boundary_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['amplitude', 'retention', 'resistance_at_peak'])
        for r in results:
            writer.writerow([r['amplitude'], r['retention'], r['resistance_at_peak']])

    print(f"Saved: {csv_file.name}")
    print()

    # Find threshold
    print("THRESHOLD ANALYSIS")
    print("=" * 70)
    print()

    stable = [r for r in results if r['retention'] > 0.99]
    unstable = [r for r in results if r['retention'] <= 0.99]

    if stable and unstable:
        min_stable = min(r['amplitude'] for r in stable)
        max_unstable = max(r['amplitude'] for r in unstable)

        print(f"✓ Threshold found between:")
        print(f"    {max_unstable:.2f} (unstable, retention={next(r for r in unstable if r['amplitude']==max_unstable)['retention']:.1%})")
        print(f"    {min_stable:.2f} (stable, retention={next(r for r in stable if r['amplitude']==min_stable)['retention']:.1%})")
        print()
        print(f"  Critical amplitude: {max_unstable:.2f} < A_crit < {min_stable:.2f}")
    elif stable:
        print("✓ All tested amplitudes stable!")
        print(f"  Stable down to: {min(r['amplitude'] for r in results):.2f}")
    else:
        print("✗ No stable configurations found")

    print()

    return results


if __name__ == "__main__":
    run_lower_boundary()
