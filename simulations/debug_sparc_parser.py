#!/usr/bin/env python3
"""Debug SPARC data parser issues."""

import numpy as np

# Read NGC2403
data = np.genfromtxt(
    'sparc_data_cache/NGC2403.dat',
    comments='#',
    names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
)

print("Sample data from NGC2403:")
print("-" * 60)
print(f"{'Rad':<8} {'Vobs':<8} {'SBdisk':<12} {'σ_disk (M/L=0.5)':<15}")
print("-" * 60)

G = 4.3e-6  # (km/s)²·kpc/M_☉

for i in range(min(5, len(data))):
    r = data['Rad'][i]
    v_obs = data['Vobs'][i]
    SB_disk = data['SBdisk'][i]

    # Convert SB to surface density (M/L = 0.5)
    sigma_disk = SB_disk * 0.5  # M_☉/pc²

    print(f"{r:<8.3f} {v_obs:<8.2f} {SB_disk:<12.3e} {sigma_disk:<15.3e}")

print()
print("Issue: SBdisk values are ~10^8 L_☉/pc²")
print("This is MUCH too high!")
print()
print("Typical spiral galaxy centers: ~1000-10000 L_☉/pc²")
print("Our generated values: ~10^8 L_☉/pc²  (10000× too high!)")
print()
print("Root cause: generate_synthetic_sparc.py converts wrong")
print("  Σ_disk should be ~1000 M_☉/pc² (center)")
print("  We're getting ~10^8 M_☉/pc²")
print()
print("Fix: Divide by 10^6 in generator (M_disk in 10^10 M_☉, area in kpc²)")
