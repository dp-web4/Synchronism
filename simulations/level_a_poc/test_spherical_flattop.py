#!/usr/bin/env python3
"""
3D Spherical Flat-Top Pattern Test

Extending 2D findings to 3D:
- 2D disk stable at amplitude ≥0.70
- Prediction: 3D sphere should be more stable (6 vs 4 neighbors)

Question: Is this prediction supported?
"""

import numpy as np
from intent_simulation import IntentGrid
from pathlib import Path
from datetime import datetime
import csv


def create_spherical_flattop(grid, center, radius, amplitude):
    """
    Create spherical flat-top pattern (3D analog of 2D disk)

    Uniform saturation inside sphere, sharp boundary
    """
    cx, cy, cz = center

    x = np.arange(grid.size)
    y = np.arange(grid.size)
    z = np.arange(grid.size)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Distance from center (periodic wrapping)
    dx = np.minimum(np.abs(X - cx), grid.size - np.abs(X - cx))
    dy = np.minimum(np.abs(Y - cy), grid.size - np.abs(Y - cy))
    dz = np.minimum(np.abs(Z - cz), grid.size - np.abs(Z - cz))

    r = np.sqrt(dx**2 + dy**2 + dz**2)

    # Spherical region
    sphere_mask = (r <= radius)

    grid.I[sphere_mask] = amplitude * grid.I_max


def test_spherical_stability(amplitude, radius=10, size=64, steps=2000):
    """Test single spherical pattern for stability"""

    print(f"Testing spherical flat-top: amplitude={amplitude:.2f}, radius={radius}")

    grid = IntentGrid(size=size, dx=1.0, dt=0.01)
    grid.n = 2  # Same as 2D tests

    # Create sphere at center
    center = (size // 2, size // 2, size // 2)
    create_spherical_flattop(grid, center, radius, amplitude)

    initial_max = grid.max_intent()
    initial_total = grid.total_intent()

    print(f"  Initial: max={initial_max:.4f}, total={initial_total:.2f}")

    # Evolve
    for step in range(steps):
        grid.step_saturating()

    final_max = grid.max_intent()
    final_total = grid.total_intent()
    retention = final_max / initial_max

    print(f"  Final:   max={final_max:.4f}, total={final_total:.2f}")
    print(f"  Retention: {retention:.1%}")

    if retention > 0.99:
        print(f"  ✓ STABLE")
    elif retention > 0.95:
        print(f"  ~ QUASI-STABLE")
    else:
        print(f"  ✗ DISSIPATED")

    print()

    return {
        'amplitude': amplitude,
        'radius': radius,
        'initial_max': initial_max,
        'final_max': final_max,
        'retention': retention,
    }


def run_3d_comparison():
    """
    Compare 3D stability to 2D results

    2D findings:
    - Stable at amplitude ≥0.70
    - Minimum radius ~10 cells

    3D prediction:
    - Should be more stable (more neighbors)
    - Possibly stable at lower amplitudes (~0.60?)
    """

    print("=" * 70)
    print("3D SPHERICAL FLAT-TOP STABILITY TEST")
    print("Comparing to 2D Results")
    print("=" * 70)
    print()

    print("2D Results (for reference):")
    print("  Amplitude 0.70: 99.5% retention (stable)")
    print("  Amplitude 0.80: 99.9% retention (stable)")
    print("  Amplitude 0.90: 100.0% retention (stable)")
    print("  Minimum radius: ~10 cells")
    print()

    print("Testing 3D spheres:")
    print()

    # Test range of amplitudes
    amplitudes = [0.90, 0.80, 0.70, 0.65, 0.60]

    results = []

    for amp in amplitudes:
        result = test_spherical_stability(
            amplitude=amp,
            radius=10,
            size=64,
            steps=2000
        )
        results.append(result)

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_file = output_dir / f"3d_spherical_stability_{timestamp}.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['amplitude', 'radius', 'initial_max', 'final_max', 'retention'])
        for r in results:
            writer.writerow([
                r['amplitude'],
                r['radius'],
                r['initial_max'],
                r['final_max'],
                r['retention'],
            ])

    print(f"Saved: {csv_file.name}")
    print()

    # Analysis
    analyze_3d_vs_2d(results)

    return results


def analyze_3d_vs_2d(results):
    """Compare 3D results to 2D expectations"""

    print("=" * 70)
    print("3D vs 2D COMPARISON")
    print("=" * 70)
    print()

    print("Results Summary:")
    print()
    print("  Amplitude | 3D Retention | 2D Retention | Comparison")
    print("  " + "-" * 60)

    # Known 2D results (approximate)
    ref_2d = {
        0.90: 1.000,
        0.80: 0.999,
        0.70: 0.995,
        0.65: None,  # Not tested in 2D
        0.60: None,
    }

    stable_3d = []
    stable_2d_count = 3  # 0.70, 0.80, 0.90 all stable in 2D

    for r in sorted(results, key=lambda x: -x['amplitude']):
        amp = r['amplitude']
        ret_3d = r['retention']
        ret_2d = ref_2d.get(amp, None)

        if ret_3d > 0.99:
            status_3d = "✓ stable"
            stable_3d.append(amp)
        else:
            status_3d = "✗ dissipated"

        if ret_2d is not None:
            status_2d = f"{ret_2d:.1%}"
            if ret_3d > ret_2d:
                comparison = "3D > 2D"
            elif ret_3d < ret_2d - 0.01:
                comparison = "3D < 2D"
            else:
                comparison = "3D ≈ 2D"
        else:
            status_2d = "not tested"
            comparison = "—"

        print(f"  {amp:.2f}     | {ret_3d:11.1%}  | {status_2d:12s} | {comparison}")

    print()

    # Verdict
    print("Key Findings:")
    print()

    min_stable_3d = min(stable_3d) if stable_3d else None

    if min_stable_3d and min_stable_3d < 0.70:
        print(f"  ✓ PREDICTION SUPPORTED")
        print(f"    3D stable at amplitude={min_stable_3d:.2f}")
        print(f"    2D stable at amplitude≥0.70")
        print(f"    3D is MORE stable than 2D (as expected)")
        print()
        print(f"  Rationale confirmed:")
        print(f"    - 6 neighbors vs 4")
        print(f"    - More resistance paths")
        print(f"    - Spherical symmetry advantage")
    elif min_stable_3d == 0.70:
        print(f"  ~ PREDICTION PARTIALLY SUPPORTED")
        print(f"    3D stable at amplitude≥0.70")
        print(f"    2D stable at amplitude≥0.70")
        print(f"    No clear difference (or difference below 0.70)")
    else:
        print(f"  ? UNEXPECTED RESULT")
        if min_stable_3d:
            print(f"    3D threshold at {min_stable_3d:.2f}")
            print(f"    Different from 2D (0.70)")
        else:
            print(f"    3D LESS STABLE than 2D!")
            print(f"    All tested amplitudes dissipated")
            print(f"    2D stable at 0.70, 3D not stable at 0.90")

    print()

    # Why 3D less stable?
    if not stable_3d or (min_stable_3d and min_stable_3d > 0.70):
        print("Why is 3D LESS stable than expected?")
        print()
        print("  Initial hypothesis was WRONG:")
        print("    'More neighbors → more resistance → more stable'")
        print()
        print("  Reality: More neighbors → more outflow paths!")
        print()
        print("  Diffusion rate comparison:")
        print("    2D Laplacian: ∇²I = (I[i+1] + I[i-1] + I[j+1] + I[j-1] - 4*I) / dx²")
        print("    3D Laplacian: ∇²I = (6 neighbor terms - 6*I) / dx²")
        print()
        print("  The '6' in 3D vs '4' in 2D means:")
        print("    - More dimensions to diffuse into")
        print("    - Faster equilibration")
        print("    - Harder to maintain concentration")
        print()
        print("  Surface-to-volume effect:")
        print("    2D: S/V ~ 2πr / πr² = 2/r")
        print("    3D: S/V ~ 4πr² / (4πr³/3) = 3/r")
        print()
        print("    3D has HIGHER surface-to-volume ratio!")
        print("    More boundary leakage relative to core")
        print()
        print("  Lesson: Dimensionality increase makes confinement HARDER")
        print("          (Like trying to hold water - easier in tray than in cloud)")

    print()


if __name__ == "__main__":
    results = run_3d_comparison()
