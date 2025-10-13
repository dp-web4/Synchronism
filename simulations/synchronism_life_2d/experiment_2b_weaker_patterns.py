#!/usr/bin/env python3
"""
Experiment 2b: Two-Body with Weaker Patterns

Previous test: amplitude=0.90, patterns frozen, no interaction
This test: amplitude=0.75, less rigid, hopefully mobile
"""

import numpy as np
from experiment_2_two_body import run_two_body_experiment, analyze_interaction
from intent_life_2d import IntentLife2D
from test_flattop import create_flattop_pattern
from pathlib import Path
from datetime import datetime
import csv


def run_weaker_two_body():
    """Test with weaker saturation"""
    print("=" * 70)
    print("EXPERIMENT 2b: TWO-BODY WITH WEAKER PATTERNS")
    print("=" * 70)
    print()

    # Weaker amplitude for more mobility
    size = 128
    amplitude = 0.75  # Lower than 0.90
    radius = 10
    separation = 40   # Closer than 50

    print(f"Configuration (adjusted for mobility):")
    print(f"  Amplitude: {amplitude:.2f} (was 0.90)")
    print(f"  Separation: {separation} cells (was 50)")
    print(f"  Rationale: Weaker saturation → less frozen → can move")
    print()

    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)
    grid.n = 2

    # Create patterns
    center1 = (size // 2 - separation // 2, size // 2)
    center2 = (size // 2 + separation // 2, size // 2)

    create_flattop_pattern(grid, center1, radius, amplitude, 'hard')
    create_flattop_pattern(grid, center2, radius, amplitude, 'hard')

    initial_distance = separation
    print(f"Initial distance: {initial_distance} cells")
    print()

    # Track metrics
    from experiment_2_two_body import detect_two_patterns

    metrics = {
        'time': [],
        'distance': [],
        'pattern1_x': [],
        'pattern1_y': [],
        'pattern2_x': [],
        'pattern2_y': [],
        'velocity': [],
    }

    steps = 5000
    save_every = 100

    print("Running simulation...")
    print()

    for step in range(steps + 1):
        if step % save_every == 0:
            c1, c2 = detect_two_patterns(grid.I, threshold=0.5 * amplitude)

            if c1 and c2:
                distance = np.sqrt((c2[0] - c1[0])**2 + (c2[1] - c1[1])**2)

                metrics['time'].append(grid.time)
                metrics['distance'].append(distance)
                metrics['pattern1_x'].append(c1[0])
                metrics['pattern1_y'].append(c1[1])
                metrics['pattern2_x'].append(c2[0])
                metrics['pattern2_y'].append(c2[1])

                if len(metrics['distance']) > 1:
                    dt_save = save_every * grid.dt
                    velocity = (metrics['distance'][-1] - metrics['distance'][-2]) / dt_save
                else:
                    velocity = 0.0

                metrics['velocity'].append(velocity)

                if step % 500 == 0:
                    print(f"Step {step:4d} / {steps} (t={grid.time:.1f})")
                    print(f"  Distance: {distance:.2f} (Δ={distance - initial_distance:+.2f})")
                    print(f"  Velocity: {velocity:.4f} cells/time")
            else:
                print(f"Step {step}: Pattern detection failed")
                break

        grid.step()

    print()
    print("Simulation complete!")
    print()

    # Save
    output_dir = Path(__file__).parent / "output"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    csv_file = output_dir / f"exp2b_weaker_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'distance', 'velocity', 'p1_x', 'p1_y', 'p2_x', 'p2_y'])
        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['distance'][i],
                metrics['velocity'][i],
                metrics['pattern1_x'][i],
                metrics['pattern1_y'][i],
                metrics['pattern2_x'][i],
                metrics['pattern2_y'][i],
            ])

    print(f"Saved: {csv_file.name}")
    print()

    # Analysis
    analyze_interaction(metrics, initial_distance)

    return grid, metrics


if __name__ == "__main__":
    run_weaker_two_body()
