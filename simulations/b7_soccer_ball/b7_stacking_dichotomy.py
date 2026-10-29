"""B7 (vacuum Umklapp) truncated-axis check (CBP-Claude, 2026-10-29): where does momentum sit relative to the zone edge G/2 = pi*hbar/l_Pl
under the two readings of 'stacking' -- total momentum of a composite vs per-mode momentum of a single quantum?"""
hbar, lp, c, eV = 1.0546e-34, 1.616e-35, 2.998e8, 1.602e-19
edge = 3.14159265 * hbar / lp
print(f"zone edge G/2 = pi*hbar/l_Pl = {edge:.1f} kg m/s")
print("TOTAL-MOMENTUM reading (stacking across constituents):")
for name, p in [("baseball 0.145 kg @ 40 m/s", 0.145 * 40), ("walking human 70 kg @ 1.4 m/s", 70 * 1.4), ("car 1500 kg @ 30 m/s", 1500 * 30), ("Earth orbital", 5.97e24 * 2.98e4)]:
    print(f"  {name:32s} p/(G/2) = {p/edge:.2g}")
print("PER-MODE reading (single quantum, substrate frame):")
for name, E in [("LHC proton 6.8 TeV", 6.8e12), ("Amaterasu-class UHECR 2.4e20 eV", 2.4e20), ("Oh-My-God particle 3.2e20 eV", 3.2e20)]:
    print(f"  {name:32s} p/(G/2) = {E*eV/c/edge:.2g}")
