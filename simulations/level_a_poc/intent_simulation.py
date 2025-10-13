#!/usr/bin/env python3
"""
Synchronism Intent Flow Simulation - Level A: Proof of Concept

Tests whether saturation resistance enables stable Intent patterns
that would otherwise dissipate under linear diffusion.

Grid: 64³ cells (smaller than planned 128³ for faster iteration)
Physics: Nonlinear diffusion with saturation resistance
Goal: Demonstrate pattern persistence vs dissipation
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import csv

# Optional matplotlib for visualization
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available - will save data to CSV only")


class IntentGrid:
    """3D Intent field with saturation dynamics"""

    def __init__(self, size=64, dx=1.0, dt=0.01):
        """
        Args:
            size: Grid dimension (size^3 cells)
            dx: Spatial step (in Planck lengths, normalized to 1.0)
            dt: Time step (in Planck times, must satisfy stability)
        """
        self.size = size
        self.dx = dx
        self.dt = dt

        # Intent field: I(x,y,z,t)
        self.I = np.zeros((size, size, size), dtype=np.float64)

        # Physical parameters
        self.D0 = 1.0        # Base diffusion coefficient
        self.I_max = 1.0     # Saturation maximum
        self.n = 2           # Resistance exponent (n=2 gives smooth transition)

        # Stability check (for explicit Euler)
        max_dt = (dx**2) / (6 * self.D0)  # 3D diffusion stability
        if dt > max_dt:
            print(f"WARNING: dt={dt} exceeds stability limit {max_dt:.6f}")

        self.time = 0.0
        self.step_count = 0

    def saturation_resistance(self, I):
        """
        Resistance function R(I) = [1 - (I/I_max)^n]

        Returns transfer resistance factor (0 = blocked, 1 = free flow)
        """
        # Clip to avoid numerical issues at boundary
        I_normalized = np.clip(I / self.I_max, 0.0, 0.999)
        return 1.0 - I_normalized**self.n

    def saturation_diffusion(self, I):
        """
        Saturation-dependent diffusion coefficient
        D(I) = D0 * R(I)
        """
        return self.D0 * self.saturation_resistance(I)

    def compute_laplacian(self, field):
        """
        Compute discrete Laplacian using 6-neighbor stencil
        ∇²f = (f[i+1] + f[i-1] + f[j+1] + f[j-1] + f[k+1] + f[k-1] - 6*f[i,j,k]) / dx²

        Uses periodic boundary conditions (wrap-around universe)
        """
        laplacian = np.zeros_like(field)

        # X-direction
        laplacian += np.roll(field, 1, axis=0)   # i-1
        laplacian += np.roll(field, -1, axis=0)  # i+1

        # Y-direction
        laplacian += np.roll(field, 1, axis=1)   # j-1
        laplacian += np.roll(field, -1, axis=1)  # j+1

        # Z-direction
        laplacian += np.roll(field, 1, axis=2)   # k-1
        laplacian += np.roll(field, -1, axis=2)  # k+1

        laplacian -= 6.0 * field
        laplacian /= self.dx**2

        return laplacian

    def step_linear(self):
        """
        Single time step with LINEAR diffusion (no saturation)

        For comparison: shows how patterns dissipate without saturation

        ∂I/∂t = D0 * ∇²I
        """
        laplacian = self.compute_laplacian(self.I)
        self.I += self.dt * self.D0 * laplacian

        # Enforce non-negativity
        self.I = np.clip(self.I, 0.0, None)

        self.time += self.dt
        self.step_count += 1

    def step_saturating(self):
        """
        Single time step with SATURATING diffusion

        This is the core Synchronism mechanism.

        ∂I/∂t = ∇ · [D(I) × ∇I]

        Simplified to: ∂I/∂t ≈ D(I) × ∇²I  (when D varies slowly)
        """
        D = self.saturation_diffusion(self.I)
        laplacian = self.compute_laplacian(self.I)

        self.I += self.dt * D * laplacian

        # Enforce saturation limit
        self.I = np.clip(self.I, 0.0, self.I_max)

        self.time += self.dt
        self.step_count += 1

    def add_gaussian_pattern(self, center=None, amplitude=0.8, sigma=5.0):
        """
        Add Gaussian Intent concentration at specified location

        Args:
            center: (x,y,z) coordinates (defaults to grid center)
            amplitude: Peak Intent (fraction of I_max)
            sigma: Width parameter (in grid cells)
        """
        if center is None:
            center = (self.size // 2, self.size // 2, self.size // 2)

        cx, cy, cz = center

        # Create coordinate grids
        x = np.arange(self.size)
        y = np.arange(self.size)
        z = np.arange(self.size)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        # Compute distance from center (with periodic wrapping)
        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))
        dz = np.minimum(np.abs(Z - cz), self.size - np.abs(Z - cz))

        r_squared = dx**2 + dy**2 + dz**2

        # Gaussian profile
        pattern = amplitude * self.I_max * np.exp(-r_squared / (2 * sigma**2))

        self.I += pattern
        self.I = np.clip(self.I, 0.0, self.I_max)

    def total_intent(self):
        """Total Intent in grid (conservation check)"""
        return np.sum(self.I)

    def max_intent(self):
        """Maximum Intent at any cell"""
        return np.max(self.I)

    def pattern_coherence(self, center, radius):
        """
        Measure coherence of pattern around center

        Returns fraction of Intent within radius of center
        """
        cx, cy, cz = center

        x = np.arange(self.size)
        y = np.arange(self.size)
        z = np.arange(self.size)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        # Distance from center
        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))
        dz = np.minimum(np.abs(Z - cz), self.size - np.abs(Z - cz))

        r = np.sqrt(dx**2 + dy**2 + dz**2)

        # Intent within radius vs total
        intent_in_radius = np.sum(self.I[r <= radius])
        total = np.sum(self.I)

        return intent_in_radius / total if total > 0 else 0.0

    def get_slice(self, axis='z', index=None):
        """Get 2D slice through grid for visualization"""
        if index is None:
            index = self.size // 2

        if axis == 'z':
            return self.I[:, :, index]
        elif axis == 'y':
            return self.I[:, index, :]
        elif axis == 'x':
            return self.I[index, :, :]


def run_comparison_experiment(size=64, steps=1000, save_every=100):
    """
    Run side-by-side comparison: linear vs saturating diffusion

    Tests core hypothesis: Saturation enables pattern persistence
    """
    print("=" * 60)
    print("Synchronism Intent Simulation - Level A: Proof of Concept")
    print("=" * 60)
    print(f"Grid size: {size}³ = {size**3:,} cells")
    print(f"Time steps: {steps}")
    print()

    # Create two grids with identical initial conditions
    grid_linear = IntentGrid(size=size, dx=1.0, dt=0.01)
    grid_saturating = IntentGrid(size=size, dx=1.0, dt=0.01)

    # Add identical Gaussian pattern at center
    center = (size // 2, size // 2, size // 2)
    amplitude = 0.8  # 80% of I_max
    sigma = 5.0      # Width

    print("Initial conditions:")
    print(f"  Pattern: Gaussian at center {center}")
    print(f"  Amplitude: {amplitude * grid_linear.I_max:.2f} (80% of I_max)")
    print(f"  Width: {sigma:.1f} cells")
    print()

    grid_linear.add_gaussian_pattern(center, amplitude, sigma)
    grid_saturating.add_gaussian_pattern(center, amplitude, sigma)

    I0_linear = grid_linear.total_intent()
    I0_saturating = grid_saturating.total_intent()

    print(f"Initial Intent: {I0_linear:.6f}")
    print()

    # Track metrics over time
    metrics = {
        'time': [],
        'linear_total': [],
        'linear_max': [],
        'linear_coherence': [],
        'saturating_total': [],
        'saturating_max': [],
        'saturating_coherence': [],
    }

    # Output directory
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("Running simulation...")
    print()

    # Main simulation loop
    for step in range(steps + 1):
        # Record metrics
        if step % save_every == 0:
            metrics['time'].append(grid_linear.time)

            metrics['linear_total'].append(grid_linear.total_intent())
            metrics['linear_max'].append(grid_linear.max_intent())
            metrics['linear_coherence'].append(grid_linear.pattern_coherence(center, radius=10))

            metrics['saturating_total'].append(grid_saturating.total_intent())
            metrics['saturating_max'].append(grid_saturating.max_intent())
            metrics['saturating_coherence'].append(grid_saturating.pattern_coherence(center, radius=10))

            print(f"Step {step:4d} / {steps} (t={grid_linear.time:.2f})")
            print(f"  Linear:     max={metrics['linear_max'][-1]:.4f}, coherence={metrics['linear_coherence'][-1]:.3f}")
            print(f"  Saturating: max={metrics['saturating_max'][-1]:.4f}, coherence={metrics['saturating_coherence'][-1]:.3f}")
            print()

        # Evolve both grids
        grid_linear.step_linear()
        grid_saturating.step_saturating()

    print("Simulation complete!")
    print()

    # Save results
    print("Saving results...")

    # Always save CSV data
    save_metrics_csv(metrics, output_dir)

    # Generate plots if matplotlib available
    if HAS_MATPLOTLIB:
        plot_comparison_results(metrics, output_dir)
        plot_field_slices(grid_linear, grid_saturating, output_dir)
    else:
        print("  (Skipping plots - matplotlib not available)")

    print(f"Results saved to: {output_dir}")
    print()

    # Summary
    print("=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print()
    print("Linear Diffusion (no saturation):")
    print(f"  Final max Intent:   {metrics['linear_max'][-1]:.6f} (started at {metrics['linear_max'][0]:.6f})")
    print(f"  Final coherence:    {metrics['linear_coherence'][-1]:.4f} (started at {metrics['linear_coherence'][0]:.4f})")
    print(f"  Pattern dissipated: {100 * (1 - metrics['linear_coherence'][-1] / metrics['linear_coherence'][0]):.1f}%")
    print()
    print("Saturating Diffusion (with saturation resistance):")
    print(f"  Final max Intent:   {metrics['saturating_max'][-1]:.6f} (started at {metrics['saturating_max'][0]:.6f})")
    print(f"  Final coherence:    {metrics['saturating_coherence'][-1]:.4f} (started at {metrics['saturating_coherence'][0]:.4f})")
    print(f"  Pattern retained:   {100 * (metrics['saturating_coherence'][-1] / metrics['saturating_coherence'][0]):.1f}%")
    print()

    if metrics['saturating_coherence'][-1] > 0.9 * metrics['saturating_coherence'][0]:
        print("✓ HYPOTHESIS CONFIRMED: Saturation enables pattern persistence")
    else:
        print("✗ HYPOTHESIS UNCLEAR: Pattern still decaying with saturation")
    print()

    return metrics, grid_linear, grid_saturating


def save_metrics_csv(metrics, output_dir):
    """Save metrics to CSV file"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = output_dir / f"metrics_{timestamp}.csv"

    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        writer.writerow([
            'time',
            'linear_max', 'linear_coherence', 'linear_total',
            'saturating_max', 'saturating_coherence', 'saturating_total'
        ])

        # Data rows
        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['linear_max'][i],
                metrics['linear_coherence'][i],
                metrics['linear_total'][i],
                metrics['saturating_max'][i],
                metrics['saturating_coherence'][i],
                metrics['saturating_total'][i],
            ])

    print(f"  Saved: {filename.name}")


def plot_comparison_results(metrics, output_dir):
    """Generate comparison plots"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Synchronism Intent Simulation: Linear vs Saturating Diffusion', fontsize=14, fontweight='bold')

    time = metrics['time']

    # Max Intent over time
    ax = axes[0, 0]
    ax.plot(time, metrics['linear_max'], 'b-', label='Linear (no saturation)', linewidth=2)
    ax.plot(time, metrics['saturating_max'], 'r-', label='Saturating (with resistance)', linewidth=2)
    ax.set_xlabel('Time (Planck units)')
    ax.set_ylabel('Maximum Intent')
    ax.set_title('Peak Intent Concentration')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pattern coherence over time
    ax = axes[0, 1]
    ax.plot(time, metrics['linear_coherence'], 'b-', label='Linear', linewidth=2)
    ax.plot(time, metrics['saturating_coherence'], 'r-', label='Saturating', linewidth=2)
    ax.set_xlabel('Time (Planck units)')
    ax.set_ylabel('Coherence (Intent within radius)')
    ax.set_title('Pattern Coherence')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])

    # Total Intent conservation
    ax = axes[1, 0]
    ax.plot(time, metrics['linear_total'], 'b-', label='Linear', linewidth=2)
    ax.plot(time, metrics['saturating_total'], 'r-', label='Saturating', linewidth=2)
    ax.set_xlabel('Time (Planck units)')
    ax.set_ylabel('Total Intent')
    ax.set_title('Intent Conservation Check')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Dissipation comparison
    ax = axes[1, 1]
    linear_retention = np.array(metrics['linear_coherence']) / metrics['linear_coherence'][0]
    saturating_retention = np.array(metrics['saturating_coherence']) / metrics['saturating_coherence'][0]

    ax.plot(time, linear_retention, 'b-', label='Linear', linewidth=2)
    ax.plot(time, saturating_retention, 'r-', label='Saturating', linewidth=2)
    ax.axhline(y=1.0, color='g', linestyle='--', alpha=0.5, label='Perfect retention')
    ax.set_xlabel('Time (Planck units)')
    ax.set_ylabel('Coherence Retention (fraction of initial)')
    ax.set_title('Pattern Stability Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])

    plt.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = output_dir / f"comparison_metrics_{timestamp}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {filename.name}")


def plot_field_slices(grid_linear, grid_saturating, output_dir):
    """Plot 2D slices through final fields"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Final Intent Field (Z-slice through center)', fontsize=14, fontweight='bold')

    # Linear diffusion result
    ax = axes[0]
    slice_linear = grid_linear.get_slice('z')
    im = ax.imshow(slice_linear.T, origin='lower', cmap='hot', vmin=0, vmax=grid_linear.I_max)
    ax.set_title('Linear Diffusion\n(Pattern dissipated)')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='Intent')

    # Saturating diffusion result
    ax = axes[1]
    slice_saturating = grid_saturating.get_slice('z')
    im = ax.imshow(slice_saturating.T, origin='lower', cmap='hot', vmin=0, vmax=grid_saturating.I_max)
    ax.set_title('Saturating Diffusion\n(Pattern persists)')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='Intent')

    plt.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = output_dir / f"field_slices_{timestamp}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {filename.name}")


if __name__ == "__main__":
    # Run Level A proof of concept
    metrics, grid_linear, grid_saturating = run_comparison_experiment(
        size=64,
        steps=1000,
        save_every=100
    )
