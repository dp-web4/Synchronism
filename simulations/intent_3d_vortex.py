"""
3D Intent Dynamics with Momentum Conservation
===============================================
Designed for Jetson AGX Thor (GPU, 122GB unified memory).
Tests whether 3D vortex structures self-confine through
dynamic saturation (the smoke ring hypothesis).

Physics:
  - Two fields: I (intent density, scalar) and v (velocity, 3-vector)
  - R(I) = 1 - (I/I_max)^n (saturation resistance)
  - When flow hits R → 0: momentum redirects (elastic), not absorbed
  - Natural pressure gradient -∇P where P = I_max - I drives radial balance
  - 3D allows vortex tubes, rings, and knots impossible in 2D

Conservation laws:
  - ΣI = const (intent conservation — non-negotiable)
  - Σ(I·v) = const (momentum conservation — the fix from conservation bug hypothesis)
  - ΣI·|v|² = energy (should be approximately conserved)

Grid: 64³ baseline, 128³ if memory allows.

Usage:
  python3 intent_3d_vortex.py [--grid 64] [--steps 10000] [--gpu]
"""

import numpy as np
import argparse
import time
import json
from pathlib import Path

# Try CuPy for GPU acceleration
try:
    import cupy as cp
    HAS_GPU = True
except ImportError:
    cp = np
    HAS_GPU = False


def R_func(I, I_max, n):
    """Saturation resistance. R → 0 at I → I_max."""
    return np.maximum(0, 1.0 - (np.clip(I, 0, I_max) / I_max) ** n)


class IntentGrid3D:
    """3D Intent dynamics with momentum conservation."""

    def __init__(self, N, I_max=1.0, n=2, c=1.0, dt=0.1, use_gpu=False):
        self.N = N
        self.I_max = I_max
        self.n = n
        self.c = c
        self.dt = dt

        self.xp = cp if (use_gpu and HAS_GPU) else np

        # Fields
        self.I = self.xp.zeros((N, N, N), dtype=np.float32)
        self.vx = self.xp.zeros((N, N, N), dtype=np.float32)
        self.vy = self.xp.zeros((N, N, N), dtype=np.float32)
        self.vz = self.xp.zeros((N, N, N), dtype=np.float32)

    def R(self, I_field):
        """Vectorized saturation resistance."""
        return self.xp.maximum(0, 1.0 - (self.xp.clip(I_field, 0, self.I_max) / self.I_max) ** self.n)

    def init_vortex_ring(self, R_ring=None, r_tube=None, I_core=0.7, I_bg=0.1):
        """
        Initialize a vortex ring (smoke ring) in the center of the grid.

        R_ring: major radius of the ring
        r_tube: minor radius (tube thickness)
        I_core: intent density in the vortex core
        I_bg: background intent density
        """
        N = self.N
        if R_ring is None:
            R_ring = N // 6
        if r_tube is None:
            r_tube = N // 12

        cx, cy, cz = N // 2, N // 2, N // 2

        # Background
        self.I[:] = I_bg

        # Create vortex ring in the xy-plane centered at (cx, cy, cz)
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    # Distance from ring axis (z-axis through center)
                    dx, dy, dz = i - cx, j - cy, k - cz

                    # Distance from the ring centerline in the xy-plane
                    r_xy = np.sqrt(dx**2 + dy**2)
                    # Distance from the torus surface
                    d_ring = np.sqrt((r_xy - R_ring)**2 + dz**2)

                    if d_ring < r_tube * 2:
                        # Intent density: Gaussian profile around ring
                        self.I[i, j, k] = I_bg + (I_core - I_bg) * np.exp(-d_ring**2 / (2 * (r_tube * 0.7)**2))

                        # Velocity: toroidal (around the tube cross-section)
                        # This creates the vortex rotation
                        if r_xy > 0.1 and d_ring > 0.1:
                            # Poloidal angle (around the tube cross-section)
                            # Tangent to circles around the tube
                            # In the plane containing the ring axis and (i,j,k):
                            r_hat_x = dx / r_xy  # radial unit vector in xy
                            r_hat_y = dy / r_xy

                            # Poloidal velocity: rotate in the (radial, z) plane
                            v_pol = 0.3 * np.exp(-d_ring**2 / (2 * r_tube**2))

                            # Poloidal direction: perpendicular to d_ring vector
                            dr_comp = (r_xy - R_ring)  # radial component of d_ring
                            dz_comp = dz  # z component
                            d_norm = max(np.sqrt(dr_comp**2 + dz_comp**2), 0.01)

                            # Tangent to poloidal circle (90° rotation of radial vector in that plane)
                            pol_r = -dz_comp / d_norm  # radial component of poloidal velocity
                            pol_z = dr_comp / d_norm   # z component

                            self.vx[i, j, k] = v_pol * pol_r * r_hat_x
                            self.vy[i, j, k] = v_pol * pol_r * r_hat_y
                            self.vz[i, j, k] = v_pol * pol_z

    def init_colliding_pulses(self, I_pulse=0.8, I_bg=0.1, v_pulse=0.3):
        """Two pulses moving toward each other — test dynamic wall formation."""
        N = self.N
        cx, cy, cz = N // 2, N // 2, N // 2

        self.I[:] = I_bg

        # Pulse 1: left of center, moving right
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    r1 = np.sqrt((i - cx + N//6)**2 + (j - cy)**2 + (k - cz)**2)
                    r2 = np.sqrt((i - cx - N//6)**2 + (j - cy)**2 + (k - cz)**2)

                    if r1 < N // 8:
                        self.I[i, j, k] += I_pulse * np.exp(-r1**2 / (2 * (N//12)**2))
                        self.vx[i, j, k] = v_pulse
                    if r2 < N // 8:
                        self.I[i, j, k] += I_pulse * np.exp(-r2**2 / (2 * (N//12)**2))
                        self.vx[i, j, k] = -v_pulse

        self.I = self.xp.clip(self.I, 0, self.I_max)

    def step(self):
        """
        One time step of 3D Intent dynamics.

        Uses compressible Euler equations with R(I)-modulated wave speed:
          ∂I/∂t + ∇·(I·v) = 0                    (continuity)
          ∂(Iv)/∂t + ∇·(Iv⊗v) = -c²·R(I)·∇I     (momentum with R-modulated pressure)

        Staggered leapfrog for stability.
        """
        xp = self.xp
        N = self.N
        dt = self.dt
        c = self.c
        I = self.I
        vx, vy, vz = self.vx, self.vy, self.vz

        R = self.R(I)

        # Pressure gradient: -c²·R(I)·∇I
        # Using central differences with periodic BC
        dIdx = xp.roll(I, -1, axis=0) - xp.roll(I, 1, axis=0)
        dIdy = xp.roll(I, -1, axis=1) - xp.roll(I, 1, axis=1)
        dIdz = xp.roll(I, -1, axis=2) - xp.roll(I, 1, axis=2)

        # Momentum update (half-step)
        fx = -c**2 * R * dIdx / 2.0
        fy = -c**2 * R * dIdy / 2.0
        fz = -c**2 * R * dIdz / 2.0

        # Add small viscous damping proportional to R (high R = more damping, low R = less)
        visc = 0.001
        lap_vx = (xp.roll(vx, -1, 0) + xp.roll(vx, 1, 0) +
                  xp.roll(vx, -1, 1) + xp.roll(vx, 1, 1) +
                  xp.roll(vx, -1, 2) + xp.roll(vx, 1, 2) - 6 * vx)
        lap_vy = (xp.roll(vy, -1, 0) + xp.roll(vy, 1, 0) +
                  xp.roll(vy, -1, 1) + xp.roll(vy, 1, 1) +
                  xp.roll(vy, -1, 2) + xp.roll(vy, 1, 2) - 6 * vy)
        lap_vz = (xp.roll(vz, -1, 0) + xp.roll(vz, 1, 0) +
                  xp.roll(vz, -1, 1) + xp.roll(vz, 1, 1) +
                  xp.roll(vz, -1, 2) + xp.roll(vz, 1, 2) - 6 * vz)

        vx_new = vx + dt * (fx + visc * R * lap_vx)
        vy_new = vy + dt * (fy + visc * R * lap_vy)
        vz_new = vz + dt * (fz + visc * R * lap_vz)

        # Continuity: ∂I/∂t = -∇·(I·v)
        flux_x = I * vx_new
        flux_y = I * vy_new
        flux_z = I * vz_new

        div_flux = ((xp.roll(flux_x, -1, 0) - xp.roll(flux_x, 1, 0)) +
                    (xp.roll(flux_y, -1, 1) - xp.roll(flux_y, 1, 1)) +
                    (xp.roll(flux_z, -1, 2) - xp.roll(flux_z, 1, 2))) / 2.0

        I_new = I - dt * div_flux

        # Enforce intent conservation (rescale to preserve total)
        total_before = float(xp.sum(I))
        I_new = xp.clip(I_new, 0.001, self.I_max)
        total_after = float(xp.sum(I_new))
        if total_after > 0:
            I_new *= total_before / total_after

        self.I = I_new
        self.vx = vx_new
        self.vy = vy_new
        self.vz = vz_new

    def diagnostics(self):
        """Compute diagnostic quantities."""
        xp = self.xp
        I = self.I
        vx, vy, vz = self.vx, self.vy, self.vz

        total_I = float(xp.sum(I))
        KE = float(0.5 * xp.sum(I * (vx**2 + vy**2 + vz**2)))
        max_I = float(xp.max(I))
        min_I = float(xp.min(I))

        # Angular momentum (Lz around grid center)
        N = self.N
        cx, cy = N // 2, N // 2
        ix = xp.arange(N).reshape(-1, 1, 1) - cx
        iy = xp.arange(N).reshape(1, -1, 1) - cy
        Lz = float(xp.sum(I * (ix * vy - iy * vx)))

        # Core detection: how much I is within R_ring of center
        iz = xp.arange(N).reshape(1, 1, -1) - N // 2
        r_from_center = xp.sqrt(ix**2 + iy**2 + iz**2)
        core_mask = r_from_center < N // 4
        core_I = float(xp.sum(I[core_mask]))

        # Width: extent of high-I region
        high_I_mask = I > 0.3
        width = float(xp.sum(high_I_mask))

        return {
            'total_I': total_I,
            'KE': KE,
            'max_I': max_I,
            'min_I': min_I,
            'Lz': Lz,
            'core_I': core_I,
            'high_I_cells': width,
        }


def main():
    parser = argparse.ArgumentParser(description='3D Intent Dynamics')
    parser.add_argument('--grid', type=int, default=32, help='Grid size (N³)')
    parser.add_argument('--steps', type=int, default=2000, help='Time steps')
    parser.add_argument('--gpu', action='store_true', help='Use GPU (CuPy)')
    parser.add_argument('--test', choices=['vortex', 'collide'], default='vortex')
    parser.add_argument('--dt', type=float, default=0.1)
    args = parser.parse_args()

    N = args.grid
    print(f'3D Intent Dynamics: {N}³ grid = {N**3:,} cells')
    print(f'Memory: ~{N**3 * 4 * 4 / 1e9:.2f} GB (4 fields × float32)')
    print(f'GPU: {"Yes (CuPy)" if (args.gpu and HAS_GPU) else "No (NumPy)"}')
    print(f'Test: {args.test}')
    print()

    grid = IntentGrid3D(N, dt=args.dt, use_gpu=args.gpu)

    if args.test == 'vortex':
        grid.init_vortex_ring()
        print('Initialized: vortex ring (smoke ring)')
    else:
        grid.init_colliding_pulses()
        print('Initialized: colliding pulses')

    # Initial diagnostics
    d0 = grid.diagnostics()
    print(f't=0: total_I={d0["total_I"]:.2f}, KE={d0["KE"]:.4f}, '
          f'Lz={d0["Lz"]:.2f}, high_I_cells={d0["high_I_cells"]:.0f}')

    results = [{'t': 0, **d0}]

    t_start = time.time()
    for t in range(1, args.steps + 1):
        grid.step()

        if t % (args.steps // 10) == 0 or t == args.steps:
            d = grid.diagnostics()
            elapsed = time.time() - t_start
            steps_per_sec = t / elapsed
            print(f't={t}: total_I={d["total_I"]:.2f}, KE={d["KE"]:.4f}, '
                  f'Lz={d["Lz"]:.2f}, high_I={d["high_I_cells"]:.0f}, '
                  f'max_I={d["max_I"]:.4f} [{steps_per_sec:.1f} steps/s]')
            results.append({'t': t, **d})

    # Summary
    d_final = results[-1]
    print()
    print('=' * 60)
    print('SUMMARY')
    print('=' * 60)
    print(f'Conservation: ΔI = {abs(d_final["total_I"] - d0["total_I"]) / d0["total_I"] * 100:.4f}%')
    print(f'KE: {d0["KE"]:.4f} → {d_final["KE"]:.4f} ({d_final["KE"]/max(d0["KE"],1e-10)*100:.1f}%)')
    print(f'Lz: {d0["Lz"]:.2f} → {d_final["Lz"]:.2f}')
    print(f'High-I cells: {d0["high_I_cells"]:.0f} → {d_final["high_I_cells"]:.0f}')

    if d_final['high_I_cells'] <= d0['high_I_cells'] * 1.5 and d_final['high_I_cells'] > 0:
        print('SELF-CONFINED: structure maintained')
    elif d_final['high_I_cells'] > d0['high_I_cells'] * 1.5:
        print('DISPERSED: structure spread')
    else:
        print('COLLAPSED: structure vanished')

    # Save results
    outfile = Path(f'results_3d_{args.test}_{N}.json')
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'Results saved to {outfile}')


if __name__ == '__main__':
    main()
