#!/usr/bin/env python3
"""
Experiment 3: Ring Patterns

Test hollow ring structure for oscillation modes.
Question: Do rings breathe, rotate, or stay static?

Rings might be more dynamic than solid disks.
"""

import numpy as np
from intent_life_2d import IntentLife2D
from test_flattop import create_flattop_pattern
from pathlib import Path
from datetime import datetime
import csv


def create_ring_pattern(grid, center, inner_radius, outer_radius, amplitude):
    """Create hollow ring pattern"""
    cx, cy = center

    x = np.arange(grid.size)
    y = np.arange(grid.size)
    X, Y = np.meshgrid(x, y, indexing='ij')

    dx = np.minimum(np.abs(X - cx), grid.size - np.abs(X - cx))
    dy = np.minimum(np.abs(Y - cy), grid.size - np.abs(Y - cy))
    r = np.sqrt(dx**2 + dy**2)

    # Ring: inner_radius < r < outer_radius
    ring_mask = (r >= inner_radius) & (r <= outer_radius)

    grid.I[ring_mask] += amplitude * grid.I_max
    grid.I = np.clip(grid.I, 0.0, grid.I_max)


def run_ring_experiment():
    """Test ring pattern behavior"""
    print("=" * 70)
    print("EXPERIMENT 3: RING PATTERNS")
    print("Testing for Oscillation Modes")
    print("=" * 70)
    print()

    size = 128
    amplitude = 0.85
    inner_radius = 12
    outer_radius = 18
    width = outer_radius - inner_radius

    print(f"Configuration:")
    print(f"  Ring inner radius: {inner_radius}")
    print(f"  Ring outer radius: {outer_radius}")
    print(f"  Width: {width}")
    print(f"  Amplitude: {amplitude:.2f}")
    print()

    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)
    grid.n = 2

    center = (size // 2, size // 2)
    create_ring_pattern(grid, center, inner_radius, outer_radius, amplitude)

    print(f"Initial total Intent: {grid.total_intent():.2f}")
    print(f"Initial max Intent: {grid.max_intent():.4f}")
    print()

    # Track metrics
    metrics = {
        'time': [],
        'max_intent': [],
        'total_intent': [],
        'mean_intent': [],
        'std_intent': [],  # Variance might indicate oscillation
    }

    steps = 3000
    save_every = 50

    print("Running simulation...")
    print("  Watching for:")
    print("    - Breathing (max Intent oscillates)")
    print("    - Collapse (ring fills inward)")
    print("    - Dispersion (ring spreads outward)")
    print()

    for step in range(steps + 1):
        if step % save_every == 0:
            metrics['time'].append(grid.time)
            metrics['max_intent'].append(grid.max_intent())
            metrics['total_intent'].append(grid.total_intent())
            metrics['mean_intent'].append(np.mean(grid.I[grid.I > 0]))
            metrics['std_intent'].append(np.std(grid.I))

            if step % 500 == 0:
                print(f"Step {step:4d} / {steps} (t={grid.time:.1f})")
                print(f"  Max Intent: {metrics['max_intent'][-1]:.4f}")
                print(f"  Total Intent: {metrics['total_intent'][-1]:.2f}")
                print(f"  Std Dev: {metrics['std_intent'][-1]:.4f}")

        grid.step()

    print()
    print("Simulation complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    csv_file = output_dir / f"exp3_ring_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'max_intent', 'total_intent', 'mean_intent', 'std_intent'])
        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['max_intent'][i],
                metrics['total_intent'][i],
                metrics['mean_intent'][i],
                metrics['std_intent'][i],
            ])

    print(f"Saved: {csv_file.name}")

    state_file = output_dir / f"exp3_ring_final_{timestamp}.npy"
    grid.save_snapshot(state_file)
    print(f"Saved: {state_file.name}")
    print()

    # Analysis
    analyze_ring_behavior(metrics)

    return grid, metrics


def analyze_ring_behavior(metrics):
    """Analyze ring dynamics"""
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print()

    initial_max = metrics['max_intent'][0]
    final_max = metrics['max_intent'][-1]
    retention = final_max / initial_max

    print(f"Stability:")
    print(f"  Initial max: {initial_max:.4f}")
    print(f"  Final max:   {final_max:.4f}")
    print(f"  Retention:   {retention:.1%}")
    print()

    if retention > 0.99:
        print("  ✓ Ring is stable")
    else:
        print("  ✗ Ring dissipating")

    print()

    # Check for oscillation
    max_series = np.array(metrics['max_intent'])
    mean_value = np.mean(max_series)
    deviations = max_series - mean_value
    oscillation_amplitude = np.std(deviations)

    print(f"Oscillation Check:")
    print(f"  Mean max Intent: {mean_value:.4f}")
    print(f"  Std deviation: {oscillation_amplitude:.6f}")
    print()

    if oscillation_amplitude > 0.001:
        print("  ⚡ OSCILLATION DETECTED")
        print(f"    Amplitude: {oscillation_amplitude:.6f}")
        print("    Ring is breathing!")

        # Try to detect period
        # Simple approach: find peaks
        from_idx = len(max_series) // 4  # Skip transient
        series_late = max_series[from_idx:]

        # Autocorrelation to find period
        # (simplified - just check if there's periodicity)
        print()
    else:
        print("  - No significant oscillation")
        print("    Ring is static")

    print()

    # Growth/collapse
    initial_total = metrics['total_intent'][0]
    final_total = metrics['total_intent'][-1]
    total_change = final_total - initial_total

    print(f"Mass Change:")
    print(f"  Initial total: {initial_total:.2f}")
    print(f"  Final total:   {final_total:.2f}")
    print(f"  Change:        {total_change:+.2f} ({100 * total_change / initial_total:+.1f}%)")
    print()

    if total_change > initial_total * 0.1:
        print("  → Ring GROWING (accumulating Intent)")
    elif total_change < -initial_total * 0.1:
        print("  → Ring SHRINKING (losing Intent)")
    else:
        print("  → Mass roughly conserved")

    print()

    # Conclusion
    print("Conclusion:")
    if retention > 0.99 and oscillation_amplitude < 0.001:
        print("  Ring behaves like static flat-top disk")
        print("  No special oscillation modes apparent")
    elif retention > 0.99 and oscillation_amplitude > 0.001:
        print("  ✓ Ring exhibits oscillation (breathing mode!)")
        print("    First non-static pattern found!")
    else:
        print("  Ring unstable, dissipating")

    print()


if __name__ == "__main__":
    grid, metrics = run_ring_experiment()
