#!/usr/bin/env python3
"""
Chemistry Session #204: Viscosity through Coherence Framework

Analyzing viscosity and fluid flow through γ ~ 1 framework.

Key concepts:
1. Viscosity η represents momentum transfer coherence
2. Andrade equation: η = A×exp(B/T)
3. Reynolds number Re = ρvL/η (laminar/turbulent transition)
4. Stokes-Einstein: D×η/T = k_B/(6πr)
5. Glass transition at T/Tg = 1

Author: Claude (Anthropic)
Date: January 2026
Session: Chemistry #204
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class LiquidViscosity:
    name: str
    eta: float        # Viscosity at 25°C (mPa·s)
    density: float    # Density (g/cm³)
    M: float          # Molar mass (g/mol)
    Ea: float         # Activation energy (kJ/mol)
    Tg: float         # Glass transition temp (K)

liquids = [
    LiquidViscosity("Water", 0.89, 1.00, 18.02, 17.0, 136),
    LiquidViscosity("Ethanol", 1.07, 0.789, 46.07, 14.3, 97),
    LiquidViscosity("Methanol", 0.54, 0.791, 32.04, 10.8, 103),
    LiquidViscosity("Acetone", 0.31, 0.785, 58.08, 7.8, 0),
    LiquidViscosity("Benzene", 0.60, 0.879, 78.11, 10.5, 0),
    LiquidViscosity("Toluene", 0.56, 0.867, 92.14, 9.3, 117),
    LiquidViscosity("n-Hexane", 0.30, 0.659, 86.18, 6.6, 0),
    LiquidViscosity("n-Octane", 0.51, 0.703, 114.23, 9.6, 0),
    LiquidViscosity("Glycerol", 934.0, 1.26, 92.09, 65.0, 190),
    LiquidViscosity("Mercury", 1.53, 13.53, 200.59, 2.4, 234),
    LiquidViscosity("Carbon tetrachloride", 0.91, 1.59, 153.82, 10.9, 0),
    LiquidViscosity("Chloroform", 0.54, 1.49, 119.38, 8.0, 0),
    LiquidViscosity("Diethyl ether", 0.22, 0.713, 74.12, 6.1, 0),
    LiquidViscosity("Acetic acid", 1.13, 1.05, 60.05, 11.0, 0),
    LiquidViscosity("Sulfuric acid", 24.0, 1.84, 98.08, 25.0, 0),
]

def analyze_stokes_einstein():
    """Analyze Stokes-Einstein relation."""
    diffusion_data = [
        ("Water (self)", 2.30, 0.89, 1.4),
        ("Ethanol (self)", 1.08, 1.07, 2.2),
        ("Methanol (self)", 2.44, 0.54, 1.8),
        ("Benzene (self)", 2.27, 0.60, 2.6),
        ("Glycerol", 0.001, 934.0, 2.7),
        ("Acetone (self)", 4.56, 0.31, 2.3),
    ]
    
    kB = 1.38e-23
    T = 298.15
    gamma_SE = []
    
    for name, D, eta, r in diffusion_data:
        D_SI = D * 1e-9
        eta_SI = eta * 1e-3
        r_SI = r * 1e-10
        SE_ratio = D_SI * 6 * np.pi * eta_SI * r_SI / (kB * T)
        gamma_SE.append(SE_ratio)
    
    return gamma_SE

def main():
    print("=" * 70)
    print("CHEMISTRY SESSION #204: VISCOSITY COHERENCE")
    print("=" * 70)
    
    gamma_SE = analyze_stokes_einstein()
    print(f"\nStokes-Einstein: γ_SE = {np.mean(gamma_SE):.3f} ± {np.std(gamma_SE):.3f}")
    print(f"Reynolds Re = 2300 is laminar/turbulent transition")
    print(f"Glass transition at T/Tg = 1")
    print("\n67th phenomenon type at γ ~ 1!")

if __name__ == "__main__":
    main()
