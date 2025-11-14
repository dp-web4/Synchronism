#!/usr/bin/env python3
"""
Generate Synthetic SPARC-like Galaxy Data
==========================================

Creates realistic galaxy rotation curve data mimicking SPARC format.
Used for testing Synchronism validation pipeline before real data access.

Key features:
- Exponential disk profiles (realistic for spirals)
- NFW-like dark matter halos (standard cosmology)
- Observational uncertainties
- Multiple galaxy types (various masses, scales)

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-14
Session: #15
"""

import numpy as np
import os
from typing import Dict, List


class SyntheticSPARCGenerator:
    """Generate synthetic galaxy data in SPARC format."""

    # Physical constants
    G = 4.3e-6  # (km/s)²·kpc/M_☉

    def __init__(self, output_dir: str = "sparc_data_cache"):
        """Initialize generator with output directory."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def generate_galaxy(
        self,
        name: str,
        M_disk: float,  # 10^10 M_☉
        R_disk: float,  # kpc
        M_halo: float,  # 10^10 M_☉
        R_halo: float,  # kpc
        n_points: int = 50
    ) -> str:
        """
        Generate one synthetic galaxy with exponential disk + NFW halo.

        Args:
            name: Galaxy identifier
            M_disk: Disk mass in 10^10 M_☉
            R_disk: Disk scale length in kpc
            M_halo: Halo mass in 10^10 M_☉
            R_halo: Halo scale radius in kpc
            n_points: Number of radial points

        Returns:
            Path to generated data file
        """
        # Radial grid (0.1 to 5 × R_disk)
        r = np.linspace(0.5, 5.0 * R_disk, n_points)

        # Exponential disk profile
        # Units: M_disk in 10^10 M_☉, R_disk in kpc
        # Result in M_☉/kpc², convert to M_☉/pc² by dividing by 10^6
        Sigma_disk = (M_disk * 1e10) / (2 * np.pi * R_disk**2 * 1e6) * np.exp(-r / R_disk)

        # NFW halo (simplified - use rotation curve directly)
        # v_halo² = G M_halo r / (r + R_halo)²
        v_halo_sq = self.G * M_halo * 1e10 * r / (r + R_halo)**2

        # Disk rotation curve (exponential disk)
        # Approximation: v_disk² ≈ G M_disk r² / R_disk³ for r << R_disk
        x = r / (2 * R_disk)
        I0 = lambda x: np.where(x < 10,
                                1 + x**2 * (1 + x**2 * (1 + x**2 * 0.25)),
                                np.exp(x) / np.sqrt(2 * np.pi * x))
        I1 = lambda x: np.where(x < 10,
                                x * (1 + x**2 * (1 + x**2 * 0.5)),
                                np.exp(x) / np.sqrt(2 * np.pi * x))

        K0 = lambda x: np.where(x < 10,
                                -np.log(x) - 0.5772 + x**2 * (0.25 + x**2 * 0.0625),
                                np.sqrt(np.pi / (2*x)) * np.exp(-x))
        K1 = lambda x: np.where(x < 10,
                                1/x + x * (0.5 + x**2 * 0.125),
                                np.sqrt(np.pi / (2*x)) * np.exp(-x))

        # Freeman (1970) exact solution
        y = r / (2 * R_disk)
        v_disk_sq = (4 * self.G * M_disk * 1e10 / R_disk) * y**2 * (
            I0(y) * K0(y) - I1(y) * K1(y)
        )**2

        # Total velocity
        v_total = np.sqrt(v_disk_sq + v_halo_sq)

        # Add observational noise (5% typical for SPARC)
        v_err = 0.05 * v_total + 3.0  # km/s (systematic + random)
        v_obs = v_total + np.random.normal(0, v_err)

        # Ensure positive velocities
        v_obs = np.abs(v_obs)

        # Gas component (simplified - 10% of disk in outer regions)
        # Already in M_☉/pc² from Sigma_disk calculation above
        Sigma_gas = 0.1 * Sigma_disk * (r / R_disk)**0.5

        # Convert to surface brightness (assume M/L = 0.5)
        SB_disk = Sigma_disk / 0.5  # L_☉/pc²

        # Bulge (20% of galaxies have significant bulge)
        has_bulge = np.random.random() < 0.2
        if has_bulge:
            M_bulge = 0.2 * M_disk
            R_bulge = 0.1 * R_disk
            # Same unit conversion: M_☉/kpc² → M_☉/pc²
            Sigma_bulge = (M_bulge * 1e10) / (2 * np.pi * R_bulge**2 * 1e6) * np.exp(-r / R_bulge)
            SB_bulge = Sigma_bulge / 0.5
        else:
            SB_bulge = np.zeros_like(r)

        # Compute component velocities (for SPARC format)
        v_gas = np.sqrt(self.G * 2 * np.pi * r * Sigma_gas * r / 1e6)  # Thin disk
        v_disk = np.sqrt(v_disk_sq)
        if has_bulge:
            v_bulge = np.sqrt(self.G * M_bulge * 1e10 * r / (r + R_bulge)**2)
        else:
            v_bulge = np.zeros_like(r)

        # Write SPARC-format file
        filename = os.path.join(self.output_dir, f"{name}.dat")
        with open(filename, 'w') as f:
            # Header (SPARC-style)
            f.write(f"# Galaxy: {name}\n")
            f.write(f"# M_disk = {M_disk:.2f} × 10^10 M_☉\n")
            f.write(f"# R_disk = {R_disk:.2f} kpc\n")
            f.write(f"# M_halo = {M_halo:.2f} × 10^10 M_☉\n")
            f.write(f"# Synthetic data for Synchronism validation\n")
            f.write("# Columns: Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n")
            f.write("# Units: kpc km/s km/s km/s km/s km/s L_☉/pc² L_☉/pc²\n")

            # Data rows
            for i in range(len(r)):
                f.write(f"{r[i]:8.3f} {v_obs[i]:8.2f} {v_err[i]:7.2f} "
                       f"{v_gas[i]:8.2f} {v_disk[i]:8.2f} {v_bulge[i]:8.2f} "
                       f"{SB_disk[i]:10.3e} {SB_bulge[i]:10.3e}\n")

        print(f"Generated: {filename}")
        return filename

    def generate_sample(self) -> List[str]:
        """Generate diverse galaxy sample mimicking SPARC variety."""

        galaxies = [
            # (name, M_disk, R_disk, M_halo, R_halo)
            # Large spirals
            ("NGC2403", 1.5, 3.5, 15.0, 10.0),
            ("NGC3198", 2.0, 4.0, 20.0, 12.0),
            ("NGC6946", 2.5, 4.5, 25.0, 15.0),
            ("NGC7331", 3.0, 5.0, 30.0, 18.0),
            ("UGC2885", 4.0, 6.0, 40.0, 20.0),

            # Medium spirals
            ("NGC925", 1.0, 3.0, 10.0, 9.0),
            ("NGC7793", 0.8, 2.5, 8.0, 8.0),
            ("NGC3621", 1.2, 3.2, 12.0, 10.0),
            ("NGC300", 0.6, 2.0, 6.0, 7.0),

            # Dwarf irregulars
            ("DDO154", 0.1, 1.0, 1.0, 3.0),
            ("IC2574", 0.15, 1.2, 1.5, 4.0),
            ("NGC2976", 0.2, 1.5, 2.0, 5.0),

            # Low surface brightness
            ("UGC128", 0.5, 4.0, 5.0, 15.0),
            ("NGC5585", 0.4, 3.5, 4.0, 12.0),

            # Compact
            ("NGC7392", 0.8, 2.0, 8.0, 6.0),
            ("NGC1560", 0.6, 2.2, 6.0, 7.0),

            # With bulge
            ("NGC5055", 2.5, 4.0, 20.0, 12.0),
            ("NGC3521", 2.0, 3.8, 18.0, 11.0),
            ("NGC4736", 1.8, 3.5, 15.0, 10.0),
            ("NGC2841", 2.8, 4.5, 25.0, 14.0),
        ]

        files = []
        for name, M_disk, R_disk, M_halo, R_halo in galaxies:
            filename = self.generate_galaxy(name, M_disk, R_disk, M_halo, R_halo)
            files.append(filename)

        print(f"\nGenerated {len(files)} synthetic galaxies")
        return files


def main():
    """Generate synthetic SPARC sample."""
    print("Generating synthetic SPARC-like galaxy data...")
    print()

    generator = SyntheticSPARCGenerator()
    files = generator.generate_sample()

    print()
    print("Synthetic SPARC dataset ready for validation!")
    print(f"Data directory: {generator.output_dir}")
    print()


if __name__ == "__main__":
    main()
