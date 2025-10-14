#!/usr/bin/env python3
"""
Test Phase-Modulated Wave Packet

This tests whether a wave packet with built-in oscillatory structure
can propagate as a traveling wave.

Key difference from simple Gaussian:
- Simple Gaussian: I(x,y) = A exp(-r²/2σ²) [no direction]
- Phase-modulated: I(x,y) = A exp(-r²/2σ²) cos(kx) [has wavelength & direction]
"""

import numpy as np
from intent_with_momentum import IntentWithMomentum


def test_phase_modulated_wave():
    """
    Test if phase-modulated wave packet propagates
    """
    print("=" * 70)
    print("TEST: PHASE-MODULATED WAVE PACKET")
    print("=" * 70)
    print()
    print("Testing wave packet with built-in wavelength and direction")
    print("I(x,y) = A exp(-r²/2σ²) cos(k·r)")
    print()

    grid = IntentWithMomentum(size=128, dx=1.0, dt=0.005)
    grid.damping = 0.001  # Minimal damping for traveling waves

    # Create phase-modulated wave packet
    center = (30, 64)
    amplitude = 0.3
    sigma = 8.0
    wavelength = 8.0  # Wavelength ~ pulse width for optimal propagation

    grid.add_wave_packet_modulated(
        center=center,
        amplitude=amplitude,
        sigma=sigma,
        wavelength=wavelength,
        direction=(1, 0)  # +x direction
    )

    print(f"Wave packet at {center}")
    print(f"Amplitude: {amplitude}, Width: {sigma:.1f}, Wavelength: {wavelength:.1f}")
    print(f"Direction: +x")
    print()

    # Track center of mass
    centers = []
    energies = []

    steps = 4000
    sample_every = 200

    for step in range(steps + 1):
        if step % sample_every == 0:
            # Find center of mass (using |I| to handle oscillating sign)
            I_abs = np.abs(grid.I)
            total = np.sum(I_abs)

            if total > 0.1:
                x = np.arange(grid.size)
                y = np.arange(grid.size)
                X, Y = np.meshgrid(x, y, indexing='ij')
                cx = np.sum(X * I_abs) / total
                cy = np.sum(Y * I_abs) / total
                centers.append((cx, cy))
            else:
                centers.append(None)

            energies.append(grid.total_energy())

            max_I = np.max(np.abs(grid.I))
            max_V = np.max(np.abs(grid.V))
            energy = energies[-1]

            print(f"Step {step:4d} (t={grid.time:.2f}): ", end='')
            if centers[-1]:
                print(f"center=({centers[-1][0]:.1f}, {centers[-1][1]:.1f}), ", end='')
            print(f"max_I={max_I:.3f}, max_V={max_V:.3f}, E={energy:.2f}")

        grid.step_wave_equation()

    print()

    # Analysis
    if len([c for c in centers if c is not None]) > 10:
        valid_centers = [c for c in centers if c is not None]
        displacement_x = valid_centers[-1][0] - valid_centers[0][0]
        displacement_y = valid_centers[-1][1] - valid_centers[0][1]

        print(f"Results:")
        print(f"  Initial center: ({valid_centers[0][0]:.1f}, {valid_centers[0][1]:.1f})")
        print(f"  Final center: ({valid_centers[-1][0]:.1f}, {valid_centers[-1][1]:.1f})")
        print(f"  X displacement: {displacement_x:.1f} cells")
        print(f"  Y displacement: {displacement_y:.1f} cells")
        print(f"  Time elapsed: {grid.time:.2f}")
        print()

        if grid.time > 0:
            avg_velocity = displacement_x / grid.time
            print(f"  Average velocity: {avg_velocity:.3f} cells/time")
            print()

            if displacement_x > 30:
                print("  ✓ SIGNIFICANT TRAVEL!")
                print("    Phase-modulated wave packet propagated successfully!")
                print("    Photon-like behavior achieved!")
            elif displacement_x > 10:
                print("  ≈ Some travel observed")
                print("    Partial success - packet moved but may have dispersed")
            elif displacement_x > 2:
                print("  ≈ Minimal travel")
                print("    Wave packet showed slight drift")
            else:
                print("  ✗ No significant travel")
                print("    Center remained essentially stationary")
    else:
        print("  ✗ Pattern dissipated completely")

    print()
    return grid


if __name__ == "__main__":
    print("Phase-Modulated Wave Packet Test")
    print()
    print("Approach: Create wave packet with oscillatory phase structure")
    print("Theory: Built-in wavelength provides direction and momentum")
    print()

    grid = test_phase_modulated_wave()

    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print()
    print("Real photons are EM waves with:")
    print("- Oscillatory E and B fields")
    print("- Defined wavelength λ")
    print("- Propagation velocity c")
    print()
    print("This test checks if adding oscillatory structure (wavelength)")
    print("enables propagation in our Intent dynamics model.")
    print()
    print("If successful: Demonstrates photon-like traveling wave capability")
    print("If unsuccessful: May need different formulation or larger scale")
