#!/usr/bin/env python3
"""
Experiment 2: Two-Body Interaction

Create two stable flat-top patterns separated by distance.
Question: Do they attract (proto-gravity)?

Test predictions from Section 4.5 (Field Effects):
- Saturated cores create gradients
- Other patterns experience transfer bias along gradients
- Should drift toward each other
"""

import numpy as np
from intent_life_2d import IntentLife2D
from test_flattop import create_flattop_pattern
from pathlib import Path
from datetime import datetime
import csv


def compute_center_of_mass(I):
    """Compute center of mass of Intent field"""
    total = np.sum(I)
    if total == 0:
        return None

    x = np.arange(I.shape[0])
    y = np.arange(I.shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')

    cx = np.sum(X * I) / total
    cy = np.sum(Y * I) / total

    return (cx, cy)


def detect_two_patterns(I, threshold=0.5, centers_guess=None):
    """
    Detect and track two separate patterns
    Returns centers of mass for each

    Simple method: Find two highest local maxima
    """
    # Find local maxima in Intent field
    # Simple approach: grid search for two peaks

    if I.max() < threshold:
        return None, None

    # Find first peak (global maximum)
    peak1_idx = np.unravel_index(np.argmax(I), I.shape)

    # Mask out region around first peak
    mask_radius = 20
    x, y = np.ogrid[:I.shape[0], :I.shape[1]]
    mask = ((x - peak1_idx[0])**2 + (y - peak1_idx[1])**2) > mask_radius**2

    # Find second peak in remaining regions
    I_masked = I.copy()
    I_masked[~mask] = 0

    if I_masked.max() < threshold:
        return None, None  # Only one pattern

    peak2_idx = np.unravel_index(np.argmax(I_masked), I.shape)

    # Compute center of mass around each peak
    def center_of_mass_local(I, peak, radius=15):
        x, y = np.ogrid[:I.shape[0], :I.shape[1]]
        mask = ((x - peak[0])**2 + (y - peak[1])**2) <= radius**2

        if not np.any(mask):
            return peak

        I_local = I * mask
        total = np.sum(I_local)

        if total == 0:
            return peak

        cx = np.sum(x * I_local) / total
        cy = np.sum(y * I_local) / total

        return (float(cx), float(cy))

    c1 = center_of_mass_local(I, peak1_idx)
    c2 = center_of_mass_local(I, peak2_idx)

    return c1, c2


def run_two_body_experiment():
    """
    Test gravitational-like attraction between patterns
    """
    print("=" * 70)
    print("EXPERIMENT 2: TWO-BODY INTERACTION")
    print("Testing Proto-Gravitational Attraction")
    print("=" * 70)
    print()

    # Configuration
    size = 128
    amplitude = 0.90  # Moderate saturation (known stable)
    radius = 10
    separation = 50   # Initial separation

    print(f"Configuration:")
    print(f"  Grid: {size} × {size}")
    print(f"  Pattern amplitude: {amplitude:.2f}")
    print(f"  Pattern radius: {radius}")
    print(f"  Initial separation: {separation} cells")
    print()

    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)
    grid.n = 2

    # Create two patterns
    center1 = (size // 2 - separation // 2, size // 2)
    center2 = (size // 2 + separation // 2, size // 2)

    print(f"Creating patterns:")
    print(f"  Pattern 1 at: {center1}")
    print(f"  Pattern 2 at: {center2}")
    print()

    create_flattop_pattern(grid, center1, radius, amplitude, 'hard')
    create_flattop_pattern(grid, center2, radius, amplitude, 'hard')

    initial_distance = np.sqrt((center2[0] - center1[0])**2 + (center2[1] - center1[1])**2)

    print(f"Initial distance: {initial_distance:.2f} cells")
    print()

    # Track metrics over time
    metrics = {
        'time': [],
        'distance': [],
        'pattern1_x': [],
        'pattern1_y': [],
        'pattern2_x': [],
        'pattern2_y': [],
        'total_intent': [],
        'velocity': [],  # Rate of distance change
    }

    steps = 5000  # Longer run to see interaction
    save_every = 100

    print("Running simulation...")
    print(f"  Steps: {steps} (t={steps * grid.dt:.1f})")
    print(f"  Watching for mutual attraction...")
    print()

    for step in range(steps + 1):
        if step % save_every == 0:
            # Track pattern positions
            c1, c2 = detect_two_patterns(grid.I, threshold=0.5 * amplitude)

            if c1 and c2:
                distance = np.sqrt((c2[0] - c1[0])**2 + (c2[1] - c1[1])**2)

                metrics['time'].append(grid.time)
                metrics['distance'].append(distance)
                metrics['pattern1_x'].append(c1[0])
                metrics['pattern1_y'].append(c1[1])
                metrics['pattern2_x'].append(c2[0])
                metrics['pattern2_y'].append(c2[1])
                metrics['total_intent'].append(grid.total_intent())

                # Compute velocity (rate of approach)
                if len(metrics['distance']) > 1:
                    dt_save = save_every * grid.dt
                    velocity = (metrics['distance'][-1] - metrics['distance'][-2]) / dt_save
                else:
                    velocity = 0.0

                metrics['velocity'].append(velocity)

                print(f"Step {step:4d} / {steps} (t={grid.time:.1f})")
                print(f"  Distance: {distance:.2f} cells (Δ={distance - initial_distance:+.2f})")
                print(f"  Velocity: {velocity:.4f} cells/time (negative = approaching)")
                print(f"  Pattern 1: ({c1[0]:.1f}, {c1[1]:.1f})")
                print(f"  Pattern 2: ({c2[0]:.1f}, {c2[1]:.1f})")
                print()
            else:
                print(f"Step {step:4d} / {steps} - WARNING: Patterns merged or split!")
                break

        grid.step()

    print("Simulation complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    csv_file = output_dir / f"exp2_two_body_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'distance', 'velocity', 'p1_x', 'p1_y', 'p2_x', 'p2_y', 'total_intent'])

        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['distance'][i],
                metrics['velocity'][i],
                metrics['pattern1_x'][i],
                metrics['pattern1_y'][i],
                metrics['pattern2_x'][i],
                metrics['pattern2_y'][i],
                metrics['total_intent'][i],
            ])

    print(f"Saved: {csv_file.name}")

    # Save final state
    state_file = output_dir / f"exp2_two_body_final_{timestamp}.npy"
    grid.save_snapshot(state_file)
    print(f"Saved: {state_file.name}")
    print()

    # Analysis
    analyze_interaction(metrics, initial_distance)

    return grid, metrics


def analyze_interaction(metrics, initial_distance):
    """Analyze whether patterns attracted"""
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print()

    if len(metrics['distance']) < 2:
        print("Insufficient data")
        return

    final_distance = metrics['distance'][-1]
    distance_change = final_distance - initial_distance

    print(f"Distance Evolution:")
    print(f"  Initial: {initial_distance:.2f} cells")
    print(f"  Final:   {final_distance:.2f} cells")
    print(f"  Change:  {distance_change:+.2f} cells ({100 * distance_change / initial_distance:+.1f}%)")
    print()

    # Average velocity
    avg_velocity = np.mean(metrics['velocity'][1:])  # Skip first (zero)

    print(f"Average velocity: {avg_velocity:.4f} cells/time")
    if avg_velocity < 0:
        print("  → Patterns APPROACHING (attractive force!)")
    elif avg_velocity > 0:
        print("  → Patterns SEPARATING (repulsive force?)")
    else:
        print("  → No net motion")
    print()

    # Check for consistent trend
    late_velocities = metrics['velocity'][len(metrics['velocity'])//2:]
    avg_late = np.mean(late_velocities) if late_velocities else 0

    print(f"Late-stage velocity: {avg_late:.4f} cells/time")
    print()

    # Verdict
    if distance_change < -1.0:  # Significant approach
        print("✓✓✓ ATTRACTION CONFIRMED!")
        print("    Patterns drifted toward each other")
        print("    Validates saturation gradient hypothesis")
        print()

        # Estimate force magnitude
        approach_rate = -avg_velocity
        print(f"Approach rate: {approach_rate:.4f} cells/time")
        print(f"At separation {initial_distance:.0f} cells")
        print()

        # Check for 1/r² (need multiple separations for this)
        print("To test F ∝ 1/r², need experiments at different separations")
        print()

    elif distance_change < -0.1:  # Small approach
        print("~ WEAK ATTRACTION")
        print("  Patterns drifted slightly toward each other")
        print("  Effect present but weak")
        print()

    elif abs(distance_change) < 0.1:  # Essentially no change
        print("✗ NO SIGNIFICANT INTERACTION")
        print("  Distance unchanged")
        print("  Either:")
        print("    - Interaction too weak at this separation")
        print("    - Patterns too stable (frozen)")
        print("    - Time scale too short")
        print()

    else:  # Separation
        print("? UNEXPECTED: SEPARATION")
        print("  Patterns moved apart")
        print("  Possible explanations:")
        print("    - Boundary growth pushes patterns apart")
        print("    - Repulsive interaction")
        print("    - Numerical artifacts")
        print()

    # Pattern trajectory
    if len(metrics['pattern1_x']) > 1:
        p1_displacement_x = metrics['pattern1_x'][-1] - metrics['pattern1_x'][0]
        p1_displacement_y = metrics['pattern1_y'][-1] - metrics['pattern1_y'][0]
        p2_displacement_x = metrics['pattern2_x'][-1] - metrics['pattern2_x'][0]
        p2_displacement_y = metrics['pattern2_y'][-1] - metrics['pattern2_y'][0]

        print(f"Pattern Trajectories:")
        print(f"  Pattern 1 displacement: ({p1_displacement_x:+.2f}, {p1_displacement_y:+.2f})")
        print(f"  Pattern 2 displacement: ({p2_displacement_x:+.2f}, {p2_displacement_y:+.2f})")
        print()

        # Check if moving toward each other
        # Expect opposite-sign x-displacement (approaching)
        if np.sign(p1_displacement_x) != np.sign(p2_displacement_x) and abs(p1_displacement_x) > 0.1:
            print("  ✓ Patterns moving toward each other along x-axis")
        else:
            print("  ? No clear approach pattern")

        print()

    # Physical interpretation
    print("Physical Interpretation:")
    print()

    if distance_change < -0.5:
        print("  Saturation gradients create attractive bias:")
        print("    - Each pattern has high Intent (saturated core)")
        print("    - Creates gradient extending outward")
        print("    - Other pattern in gradient region")
        print("    - Transfer bias toward gradient source")
        print("    - Appears as 'gravitational' attraction")
        print()
        print("  This validates Section 4.5 of whitepaper!")
    else:
        print("  If no attraction seen:")
        print("    - Patterns may be too stable (Intent locked)")
        print("    - Need longer time scale")
        print("    - Or initial separation too large")
        print("    - Or amplitude too high (boundary resistance blocks flow)")

    print()


if __name__ == "__main__":
    grid, metrics = run_two_body_experiment()
