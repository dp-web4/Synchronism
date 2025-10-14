#!/usr/bin/env python3
"""
Glider & Oscillator Search

Systematic exploration for dynamic patterns:
- Moving patterns (gliders)
- Oscillating patterns (blinkers, pulsars)
- Rotating patterns

Not looking for static stability - looking for interesting dynamics!
"""

import numpy as np
from intent_life_2d import IntentLife2D
from pathlib import Path
from datetime import datetime
import csv


def create_pattern(grid, pattern_type, center, saturation=0.80):
    """
    Create various test patterns

    Inspired by Game of Life shapes but adapted for continuous Intent field
    """
    cx, cy = center
    size = grid.size
    I_val = saturation * grid.I_max

    patterns = {
        # Asymmetric (glider candidates)
        'L_shape': [(0, 0), (1, 0), (0, 1)],
        'glider_3': [(0, 1), (1, 2), (2, 0), (2, 1), (2, 2)],  # Like GoL glider
        'arrow': [(0, 1), (1, 0), (1, 1), (1, 2), (2, 1)],  # Pointed
        'wedge': [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)],

        # Multi-lobed (oscillator candidates)
        'dumbbell': [(0, 2), (1, 2), (2, 2), (5, 2), (6, 2), (7, 2)],  # Two blobs
        'triple': [(0, 5), (1, 5), (5, 0), (5, 1), (5, 5)],  # Three peaks
        'diagonal_pair': [(0, 0), (1, 0), (5, 5), (6, 5)],

        # Rotating (angular momentum candidates)
        'off_center_disk': None,  # Special handling
        'spiral_seed': [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2)],

        # Wave packets
        'pulse_line': [(i, 2) for i in range(8)],  # Horizontal line
        'pulse_vertical': [(2, i) for i in range(8)],  # Vertical line
    }

    if pattern_type == 'off_center_disk':
        # Disk but offset from center
        for i in range(-3, 4):
            for j in range(-3, 4):
                if i**2 + j**2 <= 9:  # radius 3
                    x = (cx + i + 2) % size  # Offset by 2
                    y = (cy + j) % size
                    grid.I[x, y] = I_val

    elif pattern_type in patterns:
        offsets = patterns[pattern_type]
        for dx, dy in offsets:
            x = (cx + dx) % size
            y = (cy + dy) % size
            grid.I[x, y] = I_val


def detect_motion(position_history, threshold=0.5):
    """
    Detect if pattern is moving

    Returns: (is_moving, velocity, direction)
    """
    if len(position_history) < 10:
        return False, 0.0, None

    # Check if center of mass has moved significantly
    initial = position_history[0]
    final = position_history[-1]

    if initial is None or final is None:
        return False, 0.0, None

    displacement = np.sqrt((final[0] - initial[0])**2 + (final[1] - initial[1])**2)

    if displacement > threshold:
        time_elapsed = len(position_history) * 0.1  # Sample every 10 steps
        velocity = displacement / time_elapsed

        direction = np.array([final[0] - initial[0], final[1] - initial[1]])
        direction = direction / np.linalg.norm(direction)

        return True, velocity, direction

    return False, 0.0, None


def detect_oscillation(max_history, min_period=5, max_period=100):
    """
    Detect periodic behavior in max Intent

    Simple approach: look for repeated pattern
    """
    if len(max_history) < max_period:
        return False, None

    series = np.array(max_history)

    # Check for periodicity using autocorrelation
    for period in range(min_period, min(max_period, len(series) // 2)):
        # Compare series with itself shifted by period
        correlation = np.corrcoef(series[:-period], series[period:])[0, 1]

        if correlation > 0.95:  # Strong correlation
            return True, period

    return False, None


def compute_pattern_center(I, threshold=0.1):
    """Compute center of mass of Intent field"""
    mask = I > threshold
    if not np.any(mask):
        return None

    total = np.sum(I[mask])
    if total == 0:
        return None

    x = np.arange(I.shape[0])
    y = np.arange(I.shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')

    cx = np.sum(X[mask] * I[mask]) / total
    cy = np.sum(Y[mask] * I[mask]) / total

    return (float(cx), float(cy))


def test_pattern(pattern_type, saturation=0.80, size=64, steps=1000):
    """
    Test single pattern for dynamic behavior

    Returns dict with motion/oscillation detection results
    """
    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)
    grid.n = 2

    center = (size // 2, size // 2)
    create_pattern(grid, pattern_type, center, saturation)

    # Track over time
    position_history = []
    max_history = []

    sample_every = 10

    for step in range(steps + 1):
        if step % sample_every == 0:
            pos = compute_pattern_center(grid.I, threshold=0.1 * saturation)
            position_history.append(pos)
            max_history.append(grid.max_intent())

        grid.step()

    # Analyze results
    is_moving, velocity, direction = detect_motion(position_history)
    is_oscillating, period = detect_oscillation(max_history)

    result = {
        'pattern': pattern_type,
        'saturation': saturation,
        'is_moving': is_moving,
        'velocity': velocity if is_moving else 0.0,
        'direction': direction,
        'is_oscillating': is_oscillating,
        'period': period,
        'final_max': max_history[-1],
        'initial_max': max_history[0],
        'retention': max_history[-1] / max_history[0] if max_history[0] > 0 else 0,
    }

    return result


def run_glider_search():
    """
    Systematic search for gliders and oscillators

    Test many pattern types at various saturations
    """
    print("=" * 70)
    print("GLIDER & OSCILLATOR SEARCH")
    print("Looking for Dynamic Patterns")
    print("=" * 70)
    print()

    # Pattern types to test
    patterns = [
        'L_shape', 'glider_3', 'arrow', 'wedge',
        'dumbbell', 'triple', 'diagonal_pair',
        'off_center_disk', 'spiral_seed',
        'pulse_line', 'pulse_vertical',
    ]

    # Saturation levels to try
    saturations = [0.70, 0.80, 0.85]

    print(f"Testing {len(patterns)} patterns × {len(saturations)} saturations")
    print(f"  = {len(patterns) * len(saturations)} total tests")
    print()
    print("This will take ~5 minutes...")
    print()

    results = []

    for pattern in patterns:
        for saturation in saturations:
            print(f"Testing {pattern:20s} at {saturation:.0%}... ", end='', flush=True)

            result = test_pattern(pattern, saturation, size=64, steps=1000)
            results.append(result)

            # Quick summary
            status = []
            if result['is_moving']:
                status.append(f"⚡MOVES v={result['velocity']:.3f}")
            if result['is_oscillating']:
                status.append(f"🔄OSC period={result['period']}")
            if result['retention'] > 0.95:
                status.append("✓stable")
            elif result['retention'] < 0.70:
                status.append("✗dissipated")

            if status:
                print(" | ".join(status))
            else:
                print("static")

    print()
    print("Search complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_file = output_dir / f"glider_search_{timestamp}.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pattern', 'saturation', 'is_moving', 'velocity',
            'is_oscillating', 'period', 'retention'
        ])
        writer.writeheader()

        for r in results:
            writer.writerow({
                'pattern': r['pattern'],
                'saturation': r['saturation'],
                'is_moving': r['is_moving'],
                'velocity': r['velocity'],
                'is_oscillating': r['is_oscillating'],
                'period': r['period'],
                'retention': r['retention'],
            })

    print(f"Saved: {csv_file.name}")
    print()

    # Analysis
    analyze_search_results(results)

    return results


def analyze_search_results(results):
    """Summarize what was found"""
    print("=" * 70)
    print("SEARCH RESULTS")
    print("=" * 70)
    print()

    # Count by behavior
    moving = [r for r in results if r['is_moving']]
    oscillating = [r for r in results if r['is_oscillating']]
    stable_static = [r for r in results if r['retention'] > 0.95 and not r['is_moving'] and not r['is_oscillating']]
    dissipated = [r for r in results if r['retention'] < 0.70]

    print(f"Summary:")
    print(f"  {len(moving):2d} patterns showed MOTION (glider candidates!)")
    print(f"  {len(oscillating):2d} patterns showed OSCILLATION")
    print(f"  {len(stable_static):2d} patterns were stable but static")
    print(f"  {len(dissipated):2d} patterns dissipated")
    print()

    if moving:
        print("🎉 GLIDERS FOUND!")
        print()
        print("  Pattern                | Saturation | Velocity | Direction")
        print("  " + "-" * 65)
        for r in sorted(moving, key=lambda x: -x['velocity']):
            dir_str = f"({r['direction'][0]:.2f}, {r['direction'][1]:.2f})" if r['direction'] is not None else "—"
            print(f"  {r['pattern']:20s} | {r['saturation']:9.0%} | {r['velocity']:8.4f} | {dir_str}")
        print()
    else:
        print("✗ No gliders found")
        print("  Tried asymmetric patterns but none moved significantly")
        print()

    if oscillating:
        print("🔄 OSCILLATORS FOUND!")
        print()
        print("  Pattern                | Saturation | Period | Retention")
        print("  " + "-" * 60)
        for r in sorted(oscillating, key=lambda x: x['period'] if x['period'] else 999):
            print(f"  {r['pattern']:20s} | {r['saturation']:9.0%} | {r['period']:6d} | {r['retention']:8.1%}")
        print()
    else:
        print("✗ No oscillators found")
        print("  Tried multi-lobed patterns but none showed periodic behavior")
        print()

    # Recommendations
    print("Next Steps:")
    if moving or oscillating:
        print("  ✓ Found dynamic patterns! Investigate further:")
        if moving:
            print(f"    - Study glider mechanics")
            print(f"    - Test persistence over longer runs")
            print(f"    - Try collisions between gliders")
        if oscillating:
            print(f"    - Characterize oscillation modes")
            print(f"    - Test amplitude/phase relationships")
    else:
        print("  - Try different pattern types (more asymmetric?)")
        print("  - Try different saturation levels")
        print("  - Try giving patterns initial 'momentum'")
        print("  - Try larger patterns (current all small)")
        print("  - Consider: maybe dynamics require lower saturation?")

    print()


if __name__ == "__main__":
    results = run_glider_search()
