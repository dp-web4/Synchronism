#!/usr/bin/env python3
"""
Intent Dynamics with Two-Component Wave (Complex Field)

Key insight: Scalar wave equation doesn't support traveling waves naturally
Need TWO coupled fields (like Re/Im or E/B) for true traveling waves

This implements:
ψ = I + iθ (complex wavefunction approach)
∂I/∂t = -∇²θ
∂θ/∂t = ∇²I

Or equivalently with velocity/momentum:
∂I/∂t = V_I
∂θ/∂t = V_θ
∂V_I/∂t = tension × ∇²θ - damping × V_I
∂V_θ/∂t = tension × ∇²I - damping × V_θ
"""

import numpy as np


class IntentTwoComponent:
    """2D Intent field with two coupled components for traveling waves"""

    def __init__(self, size=128, dx=1.0, dt=0.005):
        """
        Args:
            size: Grid dimension
            dx: Spatial step
            dt: Time step
        """
        self.size = size
        self.dx = dx
        self.dt = dt

        # Two components (like real/imaginary or E/B)
        self.I = np.zeros((size, size), dtype=np.float64)  # "Real" part
        self.theta = np.zeros((size, size), dtype=np.float64)  # "Imaginary" part

        # Velocities for each component
        self.V_I = np.zeros((size, size), dtype=np.float64)
        self.V_theta = np.zeros((size, size), dtype=np.float64)

        # Parameters
        self.I_max = 1.0
        self.damping = 0.05
        self.tension = 1.0

        # CFL check
        max_dt = dx / np.sqrt(self.tension)
        if dt > max_dt:
            print(f"WARNING: dt={dt} may be unstable, max recommended: {max_dt:.6f}")

        self.time = 0.0
        self.step_count = 0

    def compute_laplacian(self, field):
        """Laplacian with periodic boundaries"""
        laplacian = np.zeros_like(field)
        laplacian += np.roll(field, 1, axis=0)
        laplacian += np.roll(field, -1, axis=0)
        laplacian += np.roll(field, 1, axis=1)
        laplacian += np.roll(field, -1, axis=1)
        laplacian -= 4.0 * field
        laplacian /= self.dx**2
        return laplacian

    def step_coupled_wave(self):
        """
        Evolve two-component system

        This creates traveling wave solutions:
        ∂I/∂t = V_I
        ∂θ/∂t = V_θ
        ∂V_I/∂t = tension × ∇²θ - damping × V_I
        ∂V_θ/∂t = -tension × ∇²I - damping × V_θ

        Note the negative sign - this creates circulation/rotation
        """
        # Compute Laplacians
        lap_I = self.compute_laplacian(self.I)
        lap_theta = self.compute_laplacian(self.theta)

        # Accelerations (note: cross-coupled!)
        acc_I = self.tension * lap_theta - self.damping * self.V_I
        acc_theta = -self.tension * lap_I - self.damping * self.V_theta

        # Update velocities
        self.V_I += self.dt * acc_I
        self.V_theta += self.dt * acc_theta

        # Update fields
        self.I += self.dt * self.V_I
        self.theta += self.dt * self.V_theta

        # Boundary conditions (allow both positive and negative)
        # Just clip magnitude
        magnitude = np.sqrt(self.I**2 + self.theta**2)
        mask_too_big = magnitude > self.I_max
        if np.any(mask_too_big):
            scale = self.I_max / magnitude[mask_too_big]
            self.I[mask_too_big] *= scale
            self.theta[mask_too_big] *= scale

        self.time += self.dt
        self.step_count += 1

    def add_wave_packet(self, center, amplitude=0.3, sigma=5.0,
                       velocity_x=0.0, velocity_y=0.0):
        """
        Create wave packet with proper initial conditions for traveling wave

        For traveling in +x direction:
        I(x,y,0) = A exp(-r²/2σ²)
        θ(x,y,0) = 0
        V_I(x,y,0) = 0
        V_θ(x,y,0) = velocity × exp(-r²/2σ²)

        The theta velocity creates the traveling wave
        """
        cx, cy = center

        x = np.arange(self.size)
        y = np.arange(self.size)
        X, Y = np.meshgrid(x, y, indexing='ij')

        dx = np.minimum(np.abs(X - cx), self.size - np.abs(X - cx))
        dy = np.minimum(np.abs(Y - cy), self.size - np.abs(Y - cy))
        r_squared = dx**2 + dy**2

        # Gaussian envelope
        envelope = amplitude * self.I_max * np.exp(-r_squared / (2 * sigma**2))

        # Initial state: all in I component
        self.I += envelope

        # Initial momentum: in theta velocity
        if velocity_x != 0 or velocity_y != 0:
            velocity_mag = np.sqrt(velocity_x**2 + velocity_y**2)
            self.V_theta += envelope * velocity_mag

    def total_energy(self):
        """Total energy in both components"""
        kinetic = 0.5 * (np.sum(self.V_I**2) + np.sum(self.V_theta**2))
        potential = 0.5 * (np.sum(self.I**2) + np.sum(self.theta**2))
        return kinetic + potential

    def amplitude(self):
        """Complex amplitude |ψ| = sqrt(I² + θ²)"""
        return np.sqrt(self.I**2 + self.theta**2)


def test_two_component_wave():
    """Test if two-component system supports traveling waves"""
    print("=" * 70)
    print("TEST: TWO-COMPONENT TRAVELING WAVE")
    print("=" * 70)
    print()
    print("Using coupled fields I and θ (like real/imaginary)")
    print("This should support traveling wave solutions")
    print()

    grid = IntentTwoComponent(size=128, dx=1.0, dt=0.005)
    grid.damping = 0.02  # Low damping

    # Create wave packet
    center = (30, 64)
    velocity = 0.5
    grid.add_wave_packet(center, amplitude=0.3, sigma=5.0, velocity_x=velocity)

    print(f"Initial wave packet at {center}")
    print(f"Target velocity: {velocity} cells/time")
    print()

    # Track center of mass
    centers = []
    energies = []

    steps = 4000
    sample_every = 200

    for step in range(steps + 1):
        if step % sample_every == 0:
            # Compute amplitude
            amp = grid.amplitude()

            # Find center of mass
            total = np.sum(amp)
            if total > 0.1:
                x = np.arange(grid.size)
                y = np.arange(grid.size)
                X, Y = np.meshgrid(x, y, indexing='ij')
                cx = np.sum(X * amp) / total
                cy = np.sum(Y * amp) / total
                centers.append((cx, cy))
            else:
                centers.append(None)

            energies.append(grid.total_energy())

            max_amp = np.max(amp)
            max_I = np.max(np.abs(grid.I))
            max_theta = np.max(np.abs(grid.theta))
            energy = energies[-1]

            print(f"Step {step:4d} (t={grid.time:.2f}): ", end='')
            if centers[-1]:
                print(f"center=({centers[-1][0]:.1f}, {centers[-1][1]:.1f}), ", end='')
            print(f"amp={max_amp:.3f}, I={max_I:.3f}, θ={max_theta:.3f}, E={energy:.2f}")

        grid.step_coupled_wave()

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

        if grid.time > 0:
            avg_velocity = displacement_x / grid.time
            print(f"  Average velocity: {avg_velocity:.3f} cells/time")
            print(f"  Target velocity: {velocity:.3f} cells/time")
            print()

            if displacement_x > 20:
                print("  ✓ WAVE PACKET TRAVELED!")
                print("    Photon-like behavior achieved!")
                if abs(avg_velocity - velocity) / velocity < 0.2:
                    print(f"    Velocity matches target within 20%")
            elif displacement_x > 5:
                print("  ≈ Some travel observed, but not fully coherent")
            else:
                print("  ✗ Minimal travel")
    else:
        print("  ✗ Pattern dissipated")

    print()
    return grid


if __name__ == "__main__":
    print("Two-Component Intent Dynamics")
    print("Testing if coupled fields enable traveling waves")
    print()

    grid = test_two_component_wave()

    print("=" * 70)
    print("This tests whether TWO coupled fields (like E/B or Re/Im)")
    print("can produce traveling wave packets where scalar field could not.")
    print("=" * 70)
