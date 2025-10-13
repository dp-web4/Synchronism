#!/usr/bin/env python3
"""
Synchronism Life 2D - Simplified Intent Dynamics

Inspired by Conway's Game of Life and artificial life research.
2D grid for faster iteration and easier visualization.
Continuous Intent values (not binary) but simplified dynamics.

Goal: Discover stable Intent patterns (analogous to GoL's gliders, oscillators, etc.)
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import csv
import json

# Optional matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import LinearSegmentedColormap
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available - limited visualization")


class IntentLife2D:
    """2D Intent field with saturation dynamics"""

    def __init__(self, size=128, dx=1.0, dt=0.01):
        """
        Args:
            size: Grid dimension (size x size cells)
            dx: Spatial step
            dt: Time step
        """
        self.size = size
        self.dx = dx
        self.dt = dt

        # Intent field: I(x,y,t)
        self.I = np.zeros((size, size), dtype=np.float64)

        # Physical parameters
        self.D0 = 1.0        # Base diffusion coefficient
        self.I_max = 1.0     # Saturation maximum
        self.n = 3           # Resistance exponent (higher = sharper cutoff)

        # Stability check
        max_dt = (dx**2) / (4 * self.D0)  # 2D diffusion stability
        if dt > max_dt:
            print(f"WARNING: dt={dt} exceeds stability limit {max_dt:.6f}")

        self.time = 0.0
        self.step_count = 0

        # Pattern tracking
        self.pattern_history = []

    def saturation_resistance(self, I):
        """R(I) = [1 - (I/I_max)^n]"""
        I_normalized = np.clip(I / self.I_max, 0.0, 0.999)
        return 1.0 - I_normalized**self.n

    def compute_laplacian(self, field):
        """
        2D Laplacian using 4-neighbor stencil
        ∇²f = (f[i+1,j] + f[i-1,j] + f[i,j+1] + f[i,j-1] - 4*f[i,j]) / dx²
        Periodic boundaries
        """
        laplacian = np.zeros_like(field)

        laplacian += np.roll(field, 1, axis=0)   # i-1
        laplacian += np.roll(field, -1, axis=0)  # i+1
        laplacian += np.roll(field, 1, axis=1)   # j-1
        laplacian += np.roll(field, -1, axis=1)  # j+1

        laplacian -= 4.0 * field
        laplacian /= self.dx**2

        return laplacian

    def step(self):
        """
        Evolve one timestep with saturation dynamics
        ∂I/∂t ≈ D(I) × ∇²I
        """
        D = self.D0 * self.saturation_resistance(self.I)
        laplacian = self.compute_laplacian(self.I)

        self.I += self.dt * D * laplacian
        self.I = np.clip(self.I, 0.0, self.I_max)

        self.time += self.dt
        self.step_count += 1

    def add_pattern(self, pattern_type, center=None, amplitude=0.95, **kwargs):
        """
        Add various pattern types for testing

        Args:
            pattern_type: 'gaussian', 'ring', 'square', 'random', 'glider_test'
            center: (x, y) coordinates
            amplitude: Peak Intent (fraction of I_max)
            **kwargs: Pattern-specific parameters
        """
        if center is None:
            center = (self.size // 2, self.size // 2)

        cx, cy = center

        if pattern_type == 'gaussian':
            sigma = kwargs.get('sigma', 5.0)
            self._add_gaussian(cx, cy, amplitude, sigma)

        elif pattern_type == 'ring':
            radius = kwargs.get('radius', 10.0)
            width = kwargs.get('width', 3.0)
            self._add_ring(cx, cy, amplitude, radius, width)

        elif pattern_type == 'square':
            size = kwargs.get('size', 10)
            self._add_square(cx, cy, amplitude, size)

        elif pattern_type == 'glider_test':
            # Asymmetric pattern that might propagate
            self._add_glider_test(cx, cy, amplitude)

        elif pattern_type == 'random':
            density = kwargs.get('density', 0.1)
            self._add_random(amplitude, density)

    def _add_gaussian(self, cx, cy, amplitude, sigma):
        """Gaussian blob"""
        x = np.arange(self.size)
        y = np.arange(self.size)
        X, Y = np.meshgrid(x, y, indexing='ij')

        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))

        r_squared = dx**2 + dy**2
        pattern = amplitude * self.I_max * np.exp(-r_squared / (2 * sigma**2))

        self.I += pattern
        self.I = np.clip(self.I, 0.0, self.I_max)

    def _add_ring(self, cx, cy, amplitude, radius, width):
        """Ring pattern (might oscillate?)"""
        x = np.arange(self.size)
        y = np.arange(self.size)
        X, Y = np.meshgrid(x, y, indexing='ij')

        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))

        r = np.sqrt(dx**2 + dy**2)

        # Gaussian ring
        pattern = amplitude * self.I_max * np.exp(-((r - radius)**2) / (2 * width**2))

        self.I += pattern
        self.I = np.clip(self.I, 0.0, self.I_max)

    def _add_square(self, cx, cy, amplitude, size):
        """Square pattern (test symmetry breaking)"""
        half = size // 2
        x_min = max(0, cx - half)
        x_max = min(self.size, cx + half)
        y_min = max(0, cy - half)
        y_max = min(self.size, cy + half)

        self.I[x_min:x_max, y_min:y_max] += amplitude * self.I_max
        self.I = np.clip(self.I, 0.0, self.I_max)

    def _add_glider_test(self, cx, cy, amplitude):
        """
        Asymmetric pattern inspired by GoL glider
        Tests if Intent patterns can propagate
        """
        # Small asymmetric blob
        offsets = [
            (0, 1),   # .#.
            (1, 2),   # ..#
            (2, 0),   # ###
            (2, 1),
            (2, 2),
        ]

        for dx, dy in offsets:
            x = (cx + dx) % self.size
            y = (cy + dy) % self.size
            self.I[x, y] += amplitude * self.I_max

        self.I = np.clip(self.I, 0.0, self.I_max)

    def _add_random(self, amplitude, density):
        """Random Intent field"""
        random_field = np.random.random((self.size, self.size))
        mask = random_field < density
        self.I[mask] += amplitude * self.I_max * random_field[mask]
        self.I = np.clip(self.I, 0.0, self.I_max)

    def total_intent(self):
        """Total Intent (conservation check)"""
        return np.sum(self.I)

    def max_intent(self):
        """Maximum Intent"""
        return np.max(self.I)

    def detect_patterns(self, threshold=0.5):
        """
        Detect coherent Intent concentrations
        Returns list of patterns with properties
        """
        try:
            from scipy import ndimage
        except ImportError:
            return []  # Skip if scipy not available

        # Binary mask of high-Intent regions
        mask = self.I > (threshold * self.I_max)

        # Connected component labeling
        labeled, num_patterns = ndimage.label(mask)

        patterns = []
        for i in range(1, num_patterns + 1):
            pattern_mask = (labeled == i)

            if np.sum(pattern_mask) < 4:  # Skip tiny patterns
                continue

            center = ndimage.center_of_mass(self.I, labeled, i)

            pattern = {
                'id': i,
                'total_intent': np.sum(self.I[pattern_mask]),
                'area': np.sum(pattern_mask),
                'center': center,
                'max_intent': np.max(self.I[pattern_mask]),
                'mean_intent': np.mean(self.I[pattern_mask]),
            }
            patterns.append(pattern)

        return patterns

    def save_snapshot(self, filename):
        """Save current state"""
        np.save(filename, self.I)


def run_pattern_test(pattern_type='gaussian', steps=2000, size=128, save_animation=False):
    """
    Test stability of specific pattern type

    Returns:
        grid: Final IntentLife2D state
        metrics: Time series data
    """
    print("=" * 60)
    print(f"Synchronism Life 2D - Pattern Test: {pattern_type}")
    print("=" * 60)
    print(f"Grid: {size} × {size} = {size**2:,} cells")
    print(f"Steps: {steps}")
    print()

    grid = IntentLife2D(size=size, dx=1.0, dt=0.01)

    # Add pattern
    if pattern_type == 'gaussian':
        grid.add_pattern('gaussian', amplitude=0.95, sigma=5.0)
        print("Pattern: Gaussian blob (amplitude=0.95, sigma=5)")

    elif pattern_type == 'ring':
        grid.add_pattern('ring', amplitude=0.95, radius=15.0, width=3.0)
        print("Pattern: Ring (radius=15, width=3)")

    elif pattern_type == 'glider_test':
        grid.add_pattern('glider_test', amplitude=0.95)
        print("Pattern: Glider-inspired asymmetric shape")

    elif pattern_type == 'random':
        grid.add_pattern('random', amplitude=0.7, density=0.15)
        print("Pattern: Random Intent field (density=0.15)")

    print()

    # Track metrics
    metrics = {
        'time': [],
        'total_intent': [],
        'max_intent': [],
        'num_patterns': [],
        'snapshots': [],  # Store field snapshots for animation
    }

    save_every = 100
    snapshot_every = 50 if save_animation else steps + 1

    print("Running simulation...")
    print()

    for step in range(steps + 1):
        if step % save_every == 0:
            patterns = grid.detect_patterns(threshold=0.5)

            metrics['time'].append(grid.time)
            metrics['total_intent'].append(grid.total_intent())
            metrics['max_intent'].append(grid.max_intent())
            metrics['num_patterns'].append(len(patterns))

            print(f"Step {step:4d} / {steps} (t={grid.time:.2f})")
            print(f"  Total Intent: {metrics['total_intent'][-1]:.2f}")
            print(f"  Max Intent:   {metrics['max_intent'][-1]:.4f}")
            print(f"  Patterns:     {metrics['num_patterns'][-1]}")
            print()

        if step % snapshot_every == 0 and save_animation:
            metrics['snapshots'].append(grid.I.copy())

        grid.step()

    print("Simulation complete!")
    print()

    # Save results
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save metrics CSV
    csv_file = output_dir / f"{pattern_type}_metrics_{timestamp}.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'total_intent', 'max_intent', 'num_patterns'])
        for i in range(len(metrics['time'])):
            writer.writerow([
                metrics['time'][i],
                metrics['total_intent'][i],
                metrics['max_intent'][i],
                metrics['num_patterns'][i],
            ])

    print(f"Saved: {csv_file.name}")

    # Save final state
    state_file = output_dir / f"{pattern_type}_final_{timestamp}.npy"
    grid.save_snapshot(state_file)
    print(f"Saved: {state_file.name}")

    # Generate plots if matplotlib available
    if HAS_MATPLOTLIB:
        plot_results(grid, metrics, pattern_type, output_dir, timestamp)

        if save_animation and len(metrics['snapshots']) > 1:
            create_animation(metrics['snapshots'], pattern_type, output_dir, timestamp)

    print()
    return grid, metrics


def plot_results(grid, metrics, pattern_type, output_dir, timestamp):
    """Generate result plots"""
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # Final field
    ax = fig.add_subplot(gs[0:2, 0])
    im = ax.imshow(grid.I.T, origin='lower', cmap='hot', vmin=0, vmax=grid.I_max)
    ax.set_title(f'Final Intent Field (t={grid.time:.1f})')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='Intent')

    # Metrics over time
    time = metrics['time']

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(time, metrics['max_intent'], 'r-', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Max Intent')
    ax.set_title('Peak Concentration')
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(time, metrics['num_patterns'], 'b-', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Pattern Count')
    ax.set_title('Number of Coherent Patterns')
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[2, :])
    ax.plot(time, metrics['total_intent'], 'g-', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Total Intent')
    ax.set_title('Conservation Check')
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Synchronism Life 2D: {pattern_type}', fontsize=14, fontweight='bold')

    filename = output_dir / f"{pattern_type}_results_{timestamp}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saved: {filename.name}")


def create_animation(snapshots, pattern_type, output_dir, timestamp):
    """Create animation of field evolution"""
    fig, ax = plt.subplots(figsize=(8, 8))

    im = ax.imshow(snapshots[0].T, origin='lower', cmap='hot', vmin=0, vmax=1.0, animated=True)
    ax.set_title('Intent Field Evolution')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='Intent')

    def update(frame):
        im.set_array(snapshots[frame].T)
        ax.set_title(f'Intent Field (frame {frame}/{len(snapshots)-1})')
        return [im]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots), interval=50, blit=True)

    filename = output_dir / f"{pattern_type}_animation_{timestamp}.gif"
    anim.save(filename, writer='pillow', fps=20)
    plt.close()

    print(f"Saved: {filename.name}")


if __name__ == "__main__":
    import sys

    pattern_type = sys.argv[1] if len(sys.argv) > 1 else 'gaussian'

    # Run pattern test
    grid, metrics = run_pattern_test(
        pattern_type=pattern_type,
        steps=2000,
        size=128,
        save_animation=False  # Set True to generate animations (slow)
    )

    # Summary
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Pattern type: {pattern_type}")
    print(f"Initial max Intent: {metrics['max_intent'][0]:.4f}")
    print(f"Final max Intent:   {metrics['max_intent'][-1]:.4f}")
    print(f"Retention:          {100 * metrics['max_intent'][-1] / metrics['max_intent'][0]:.1f}%")
    print(f"Final patterns:     {metrics['num_patterns'][-1]}")
    print()

    if metrics['max_intent'][-1] > 0.9 * metrics['max_intent'][0]:
        print("✓ STABLE: Pattern persisted with >90% retention")
    elif metrics['max_intent'][-1] > 0.7 * metrics['max_intent'][0]:
        print("~ QUASI-STABLE: Slow dissipation")
    else:
        print("✗ DISSIPATED: Pattern did not stabilize")
    print()
