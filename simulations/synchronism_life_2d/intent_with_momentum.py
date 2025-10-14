#!/usr/bin/env python3
"""
Intent Dynamics with Momentum/Tension

Key change: Intent has VELOCITY not just concentration
- Gradients create FORCE (tension)
- Force causes ACCELERATION
- Enables oscillation and waves

This is wave equation with damping, not diffusion equation.
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import csv


class IntentWithMomentum:
    """2D Intent field with momentum dynamics"""

    def __init__(self, size=128, dx=1.0, dt=0.005):
        """
        Args:
            size: Grid dimension
            dx: Spatial step
            dt: Time step (needs to be smaller for stability!)
        """
        self.size = size
        self.dx = dx
        self.dt = dt

        # Intent field
        self.I = np.zeros((size, size), dtype=np.float64)

        # NEW: Velocity field (Intent flow velocity)
        self.V = np.zeros((size, size), dtype=np.float64)

        # Parameters
        self.I_max = 1.0
        self.damping = 0.1    # Damping coefficient (prevents infinite oscillation)
        self.tension = 1.0    # How strongly gradients create force

        # For stability: wave equation CFL condition
        # dt < dx / sqrt(tension)
        max_dt = dx / np.sqrt(self.tension)
        if dt > max_dt:
            print(f"WARNING: dt={dt} may be unstable, max recommended: {max_dt:.6f}")

        self.time = 0.0
        self.step_count = 0

    def compute_gradient_force(self):
        """
        Compute force on each cell from Intent gradient

        F = -∇I (force points down gradient, toward lower Intent)

        Returns: Force field (same shape as I)
        """
        # Central differences for gradient
        dI_dx = (np.roll(self.I, -1, axis=0) - np.roll(self.I, 1, axis=0)) / (2 * self.dx)
        dI_dy = (np.roll(self.I, -1, axis=1) - np.roll(self.I, 1, axis=1)) / (2 * self.dx)

        # Force magnitude (gradient magnitude)
        # Negative because force points DOWN gradient
        force = -np.sqrt(dI_dx**2 + dI_dy**2)

        return force

    def compute_laplacian(self):
        """Laplacian for wave equation"""
        laplacian = np.zeros_like(self.I)

        # 4-neighbor stencil
        laplacian += np.roll(self.I, 1, axis=0)
        laplacian += np.roll(self.I, -1, axis=0)
        laplacian += np.roll(self.I, 1, axis=1)
        laplacian += np.roll(self.I, -1, axis=1)
        laplacian -= 4.0 * self.I
        laplacian /= self.dx**2

        return laplacian

    def step_wave_equation(self):
        """
        Evolve using wave equation with damping

        ∂²I/∂t² = c² ∇²I - γ ∂I/∂t

        Where:
        - c² = tension (wave speed squared)
        - γ = damping (energy dissipation)

        Rewritten as two first-order equations:
        ∂I/∂t = V
        ∂V/∂t = c² ∇²I - γ V
        """
        # Compute acceleration from Laplacian (wave propagation)
        laplacian = self.compute_laplacian()
        acceleration = self.tension * laplacian

        # Subtract damping
        acceleration -= self.damping * self.V

        # Update velocity
        self.V += self.dt * acceleration

        # Update Intent
        self.I += self.dt * self.V

        # Enforce non-negativity (Intent can't be negative)
        # When hitting zero, reverse velocity (reflection)
        mask_negative = self.I < 0
        self.I[mask_negative] = 0
        self.V[mask_negative] = -self.V[mask_negative] * 0.5  # Partial reflection

        # Enforce saturation maximum
        mask_saturated = self.I > self.I_max
        self.I[mask_saturated] = self.I_max
        self.V[mask_saturated] = 0  # Stop at saturation

        self.time += self.dt
        self.step_count += 1

    def add_pulse(self, center, amplitude=0.5, sigma=5.0, velocity_x=0.0, velocity_y=0.0):
        """
        Add Gaussian pulse with optional initial velocity

        This creates a wave packet that should travel
        """
        cx, cy = center

        x = np.arange(self.size)
        y = np.arange(self.size)
        X, Y = np.meshgrid(x, y, indexing='ij')

        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))
        r_squared = dx**2 + dy**2

        # Gaussian profile
        pulse = amplitude * self.I_max * np.exp(-r_squared / (2 * sigma**2))

        self.I += pulse

        # Add velocity (momentum) - give entire pulse uniform velocity
        if velocity_x != 0 or velocity_y != 0:
            # Simple approach: velocity field proportional to pulse amplitude
            # Entire wave packet moves together with specified velocity
            # This is analogous to p = mv where pulse is the "mass"
            velocity_magnitude = np.sqrt(velocity_x**2 + velocity_y**2)
            self.V += pulse * velocity_magnitude

    def add_wave_packet_modulated(self, center, amplitude=0.3, sigma=5.0,
                                   wavelength=10.0, direction=(1, 0)):
        """
        Add wave packet with phase modulation (oscillatory structure)

        This has built-in wavelength and direction, like a real photon:
        I(x,y) = A exp(-r²/2σ²) cos(k·r)

        Args:
            center: (cx, cy) center position
            amplitude: Peak amplitude
            sigma: Gaussian width
            wavelength: Wavelength of oscillation
            direction: (dx, dy) direction of propagation
        """
        cx, cy = center
        dx_dir, dy_dir = direction

        # Normalize direction
        dir_mag = np.sqrt(dx_dir**2 + dy_dir**2)
        if dir_mag > 0:
            dx_dir /= dir_mag
            dy_dir /= dir_mag

        x = np.arange(self.size)
        y = np.arange(self.size)
        X, Y = np.meshgrid(x, y, indexing='ij')

        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))

        # Distance from center
        r_squared = dx**2 + dy**2

        # Position along propagation direction
        # (this determines phase)
        r_parallel = (X - cx) * dx_dir + (Y - cy) * dy_dir

        # Wave number
        k = 2 * np.pi / wavelength

        # Gaussian envelope × oscillatory phase
        envelope = amplitude * self.I_max * np.exp(-r_squared / (2 * sigma**2))
        phase = np.cos(k * r_parallel)

        self.I += envelope * phase

        # Initial velocity: 90° out of phase (sine instead of cosine)
        # This creates the traveling wave initial condition
        velocity = (2 * np.pi / wavelength)  # ω/k for wave
        phase_velocity = np.sin(k * r_parallel)
        self.V += envelope * phase_velocity * velocity

    def add_oscillator_seed(self, center, amplitude=0.6, size=5):
        """
        Create pattern that might oscillate

        Central high region with velocity perturbation
        """
        cx, cy = center

        # Create square region
        for i in range(-size, size+1):
            for j in range(-size, size+1):
                x = (cx + i) % self.size
                y = (cy + j) % self.size
                self.I[x, y] = amplitude * self.I_max

        # Add small velocity perturbation (kick it into oscillation)
        self.V[cx, cy] = 0.1

    def total_energy(self):
        """
        Total energy: kinetic + potential

        KE = (1/2) Σ V²
        PE = (1/2) Σ I²  (stored in Intent field)
        """
        kinetic = 0.5 * np.sum(self.V**2)
        potential = 0.5 * np.sum(self.I**2)
        return kinetic + potential


def test_wave_packet():
    """
    Test 1: Can we create a traveling wave packet (photon analog)?
    """
    print("=" * 70)
    print("TEST 1: TRAVELING WAVE PACKET (Photon Analog)")
    print("=" * 70)
    print()

    grid = IntentWithMomentum(size=128, dx=1.0, dt=0.005)
    grid.damping = 0.02  # Lower damping for traveling waves

    # Create pulse with velocity - smaller amplitude to avoid saturation
    center = (30, 64)
    grid.add_pulse(center, amplitude=0.3, sigma=5.0, velocity_x=0.5)

    print(f"Initial pulse at {center} with velocity_x=2.0")
    print(f"Watching for: Travel to right while maintaining shape")
    print()

    # Track center of mass
    centers = []
    energies = []

    steps = 2000
    sample_every = 100

    for step in range(steps + 1):
        if step % sample_every == 0:
            # Find center of mass
            total = np.sum(grid.I)
            if total > 0.1:
                x = np.arange(grid.size)
                y = np.arange(grid.size)
                X, Y = np.meshgrid(x, y, indexing='ij')
                cx = np.sum(X * grid.I) / total
                cy = np.sum(Y * grid.I) / total
                centers.append((cx, cy))
            else:
                centers.append(None)

            energies.append(grid.total_energy())

            max_I = np.max(grid.I)
            max_V = np.max(np.abs(grid.V))
            energy = energies[-1]

            print(f"Step {step:4d} (t={grid.time:.2f}): ", end='')
            if centers[-1]:
                print(f"center=({centers[-1][0]:.1f}, {centers[-1][1]:.1f}), ", end='')
            print(f"max_I={max_I:.3f}, max_V={max_V:.3f}, energy={energy:.2f}")

        grid.step_wave_equation()

    print()

    # Analysis
    if len([c for c in centers if c is not None]) > 10:
        valid_centers = [c for c in centers if c is not None]
        displacement_x = valid_centers[-1][0] - valid_centers[0][0]
        displacement_y = valid_centers[-1][1] - valid_centers[0][1]

        print(f"Results:")
        print(f"  X displacement: {displacement_x:.1f} cells")
        print(f"  Y displacement: {displacement_y:.1f} cells")
        print(f"  Time elapsed: {grid.time:.2f}")
        print(f"  Average velocity: {displacement_x / grid.time:.3f} cells/time")
        print()

        if displacement_x > 20:
            print("  ✓ WAVE PACKET TRAVELED!")
            print("    Photon-like behavior observed")
        else:
            print("  ✗ Minimal travel, packet dispersed or damped")
    else:
        print("  ✗ Pattern dissipated completely")

    print()
    return grid


def test_oscillator():
    """
    Test 2: Can we create stable oscillator (electron analog)?
    """
    print("=" * 70)
    print("TEST 2: STABLE OSCILLATOR (Electron Analog)")
    print("=" * 70)
    print()

    grid = IntentWithMomentum(size=128, dx=1.0, dt=0.005)
    grid.damping = 0.05  # Lower damping for sustained oscillation

    # Create oscillator seed
    center = (64, 64)
    grid.add_oscillator_seed(center, amplitude=0.6, size=5)

    print(f"Oscillator seed at {center}")
    print(f"Damping = {grid.damping:.3f} (low for sustained oscillation)")
    print(f"Watching for: Periodic oscillation of Intent at center")
    print()

    # Track Intent at center over time
    center_values = []
    energies = []

    steps = 3000
    sample_every = 20

    for step in range(steps + 1):
        if step % sample_every == 0:
            center_values.append(grid.I[center[0], center[1]])
            energies.append(grid.total_energy())

            if step % 500 == 0:
                print(f"Step {step:4d} (t={grid.time:.2f}): center_I={center_values[-1]:.4f}, energy={energies[-1]:.2f}")

        grid.step_wave_equation()

    print()

    # Detect oscillation
    series = np.array(center_values[10:])  # Skip initial transient

    # Simple period detection: autocorrelation
    if len(series) > 50:
        mean = np.mean(series)
        variance = np.var(series)

        print(f"Results:")
        print(f"  Mean Intent at center: {mean:.4f}")
        print(f"  Variance: {variance:.6f}")
        print()

        if variance > 0.001:
            print("  ✓ OSCILLATION DETECTED!")
            print("    Electron-like behavior: localized oscillating pattern")

            # Try to find period
            for period in range(5, 100):
                if period * 2 > len(series):
                    break
                corr = np.corrcoef(series[:-period], series[period:])[0, 1]
                if corr > 0.9:
                    print(f"    Period ≈ {period * sample_every * grid.dt:.2f} time units")
                    break
        else:
            print("  ✗ No oscillation, pattern settled to equilibrium")
    else:
        print("  ? Insufficient data")

    print()
    return grid


def test_binding():
    """
    Test 3: Can two oscillators bind at stable separation (atom analog)?
    """
    print("=" * 70)
    print("TEST 3: BINDING (Atom Analog)")
    print("=" * 70)
    print()
    print("Two oscillators at different saturations (nucleus + electron analog)")
    print()

    grid = IntentWithMomentum(size=128, dx=1.0, dt=0.005)
    grid.damping = 0.05

    # Heavy "nucleus" (high saturation, large)
    center1 = (54, 64)
    grid.add_oscillator_seed(center1, amplitude=0.8, size=8)

    # Light "electron" (lower saturation, small)
    center2 = (74, 64)
    grid.add_oscillator_seed(center2, amplitude=0.5, size=3)

    initial_separation = 20

    print(f"Pattern 1 (nucleus): center {center1}, amplitude=0.8, size=8")
    print(f"Pattern 2 (electron): center {center2}, amplitude=0.5, size=3")
    print(f"Initial separation: {initial_separation} cells")
    print()

    # Track separation over time
    separations = []

    steps = 3000
    sample_every = 100

    for step in range(steps + 1):
        if step % sample_every == 0:
            # Find centers of mass for each pattern
            # (Simplified: just use max locations)
            peak1_idx = np.unravel_index(np.argmax(grid.I[:70, :]), grid.I[:70, :].shape)
            peak2_idx = np.unravel_index(np.argmax(grid.I[70:, :]), grid.I[70:, :].shape)
            peak2_idx = (peak2_idx[0] + 70, peak2_idx[1])

            separation = np.sqrt((peak2_idx[0] - peak1_idx[0])**2 + (peak2_idx[1] - peak1_idx[1])**2)
            separations.append(separation)

            if step % 500 == 0:
                print(f"Step {step:4d} (t={grid.time:.2f}): separation={separation:.1f} cells")

        grid.step_wave_equation()

    print()

    # Analysis
    final_separation = separations[-1]
    delta = final_separation - initial_separation

    print(f"Results:")
    print(f"  Initial separation: {initial_separation:.1f}")
    print(f"  Final separation: {final_separation:.1f}")
    print(f"  Change: {delta:+.1f} cells")
    print()

    if abs(delta) < 2:
        print("  ✓ STABLE BINDING!")
        print("    Patterns maintained separation (atom-like)")
    elif delta < -5:
        print("  → Patterns attracted, may have merged")
    else:
        print("  ✗ No clear binding, patterns dispersed")

    print()
    return grid


if __name__ == "__main__":
    print("Intent Dynamics with Momentum/Tension")
    print("Testing for particle-like behaviors")
    print()

    # Run all three tests
    grid1 = test_wave_packet()
    grid2 = test_oscillator()
    grid3 = test_binding()

    print("=" * 70)
    print("ALL TESTS COMPLETE")
    print("=" * 70)
    print()
    print("These tests probe whether wave equation dynamics")
    print("can produce photons (traveling waves), electrons (oscillators),")
    print("and atoms (bound systems).")
    print()
    print("This is the right question for Synchronism at particle scale!")
