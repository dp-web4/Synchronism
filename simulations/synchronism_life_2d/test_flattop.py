#!/usr/bin/env python3
"""
Test Flat-Top Pattern Stability

Hypothesis: Gaussian patterns dissipate because they have gradients.
Test: Uniformly saturated core with sharp boundary (step function)

This is more like a quantum particle (box wavefunction) than Gaussian.
"""

import numpy as np
from intent_life_2d import IntentLife2D
from pathlib import Path
from datetime import datetime
import csv


def create_flattop_pattern(grid, center, radius, amplitude, sharpness='hard'):
    """
    Create flat-top pattern: saturated core with sharp boundary

    Args:
        grid: IntentLife2D instance
        center: (x, y) center coordinates
        radius: Core radius
        amplitude: Intent value in core (fraction of I_max)
        sharpness: 'hard' (step function) or 'soft' (tanh transition)
    """
    cx, cy = center

    x = np.arange(grid.size)
    y = np.arange(grid.size)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Distance from center (with periodic wrapping)
    dx = np.minimum(np.abs(X - cx), grid.size - np.abs(X - cx))
    dy = np.minimum(np.abs(Y - cy), grid.size - np.abs(Y - cy))
    r = np.sqrt(dx**2 + dy**2)

    if sharpness == 'hard':
        # Step function: 1 inside radius, 0 outside
        pattern = (r <= radius).astype(float) * amplitude * grid.I_max

    elif sharpness == 'soft':
        # Smooth transition using tanh
        # tanh((radius - r) / width) goes from 1 (inside) to 0 (outside)
        width = 2.0  # Transition width
        pattern = amplitude * grid.I_max * 0.5 * (1 + np.tanh((radius - r) / width))

    grid.I += pattern
    grid.I = np.clip(grid.I, 0.0, grid.I_max)


def test_flattop(n=2, amplitude=0.999, radius=15, sharpness='hard', steps=3000):
    """
    Test flat-top pattern stability
    """
    print("=" * 70)
    print("FLAT-TOP PATTERN STABILITY TEST")
    print("=" * 70)
    print(f"Configuration:")
    print(f"  n (resistance exponent): {n}")
    print(f"  amplitude: {amplitude:.4f}")
    print(f"  radius: {radius}")
    print(f"  sharpness: {sharpness}")
    print(f"  steps: {steps}")
    print()

    grid = IntentLife2D(size=128, dx=1.0, dt=0.01)
    grid.n = n

    # Create flat-top pattern
    center = (grid.size // 2, grid.size // 2)
    create_flattop_pattern(grid, center, radius, amplitude, sharpness)

    initial_max = grid.max_intent()
    initial_total = grid.total_intent()
    initial_resistance = 1 - amplitude**n

    print(f"Initial state:")
    print(f"  Max Intent: {initial_max:.6f}")
    print(f"  Total Intent: {initial_total:.2f}")
    print(f"  Resistance at peak: {initial_resistance:.1%}")
    print(f"  Core area: π×{radius}² ≈ {int(np.pi * radius**2)} cells")
    print()

    # Track metrics
    metrics = {
        'time': [],
        'max_intent': [],
        'total_intent': [],
        'patterns': [],
    }

    save_every = 100

    print("Running simulation...")
    print()

    for step in range(steps + 1):
        if step % save_every == 0:
            patterns = grid.detect_patterns(threshold=0.5)

            metrics['time'].append(grid.time)
            metrics['max_intent'].append(grid.max_intent())
            metrics['total_intent'].append(grid.total_intent())
            metrics['patterns'].append(len(patterns))

            retention = metrics['max_intent'][-1] / initial_max

            print(f"Step {step:4d} / {steps} (t={grid.time:.2f})")
            print(f"  Max Intent:   {metrics['max_intent'][-1]:.6f} (retention={retention:.1%})")
            print(f"  Total Intent: {metrics['total_intent'][-1]:.2f}")
            print(f"  Patterns:     {metrics['patterns'][-1]}")
            print()

        grid.step()

    print("Simulation complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save metrics
    csv_file = output_dir / f"flattop_n{n}_amp{amplitude:.4f}_{sharpness}_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'max_intent', 'total_intent', 'patterns'])
        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['max_intent'][i],
                metrics['total_intent'][i],
                metrics['patterns'][i],
            ])

    print(f"Saved: {csv_file.name}")

    # Save final state
    state_file = output_dir / f"flattop_n{n}_amp{amplitude:.4f}_{sharpness}_final_{timestamp}.npy"
    grid.save_snapshot(state_file)
    print(f"Saved: {state_file.name}")
    print()

    # Analysis
    final_max = metrics['max_intent'][-1]
    retention = final_max / initial_max

    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Initial max Intent: {initial_max:.6f}")
    print(f"Final max Intent:   {final_max:.6f}")
    print(f"Retention:          {retention:.1%}")
    print()

    if retention > 0.95:
        print("✓✓✓ STABLE: Pattern maintained >95% of peak!")
        print("    This is the first stable Intent pattern discovered!")
    elif retention > 0.90:
        print("✓✓ QUASI-STABLE: Pattern maintained >90%")
        print("    Near stability threshold")
    elif retention > 0.80:
        print("✓ SLOW DISSIPATION: Better than Gaussian")
        print("    But not fully stable")
    else:
        print("✗ DISSIPATED: Similar to Gaussian patterns")
        print("    Flat-top hypothesis not confirmed")
    print()

    # Compare to Gaussian
    print("Comparison to Gaussian (from parameter sweep):")
    print(f"  Gaussian n={n}, amp={amplitude:.2f}: ~{55.9 if n==2 else 47.1}% retention at t=20")
    print(f"  Flat-top n={n}, amp={amplitude:.4f}: {retention:.1%} retention at t={metrics['time'][-1]:.1f}")
    print()

    return grid, metrics


if __name__ == "__main__":
    import sys

    # Parse arguments
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    amplitude = float(sys.argv[2]) if len(sys.argv) > 2 else 0.999
    sharpness = sys.argv[3] if len(sys.argv) > 3 else 'hard'

    # Run test
    grid, metrics = test_flattop(
        n=n,
        amplitude=amplitude,
        radius=15,
        sharpness=sharpness,
        steps=3000  # Longer than Gaussian tests
    )
