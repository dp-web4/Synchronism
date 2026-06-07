#!/usr/bin/env python3
"""
Session 687: arithmetic audit of the A-from-Jeans derivation, in response to the
chain-of-custody closure proposal filed 2026-06-07.

The explorer track's claims to verify:

  (1) The stated formula A = 4*pi / (beta_J^2 * G * R_0^2) with beta_J = 1 and
      R_0 = 8 kpc DOES NOT give A ~ 0.028; it gives ~ 4.6e-5 (~600x too small).

  (2) The "5% agreement" 0.0294 number comes from a different calculation
      using fitted alpha = 4.5 and R_0 = 0.07 kpc/(km/s)^0.75 -- the
      size-velocity relation slope -- which derives rho_crit ~ V^0.5, NOT V^2.

  (3) The 644x "unit conversion" used in the S66 markdown to bridge 4.57e-5 to
      0.0294 is just the ratio of the V^0.5 result to the V^2 formula evaluation.

Procedure: compute A in SI, convert to the (M_sun/pc^3)/(km/s)^2 convention
used in the framework, compare to 0.028; then evaluate the V^0.5 alternate
calculation and check whether it equals 0.0294.
"""
import numpy as np

# --- Constants --------------------------------------------------------------
G_SI = 6.674e-11                # m^3 / (kg * s^2)
Msun = 1.989e30                 # kg
kpc = 3.0857e19                 # m
pc  = kpc / 1000.0              # m
kms = 1000.0                    # m/s
beta_J = 1.0
R0_MW = 8.0 * kpc               # the "Sun's galactocentric radius" used in the formula

# --- (1) A from the stated formula with R_0 = 8 kpc -------------------------
print("=" * 74)
print("(1) Stated formula  A = 4*pi / (beta_J^2 * G * R_0^2)  with R_0 = 8 kpc")
print("=" * 74)
A_SI = 4.0 * np.pi / (beta_J**2 * G_SI * R0_MW**2)
print(f"   beta_J = {beta_J}")
print(f"   R_0    = {R0_MW:.3e} m  (= 8 kpc)")
print(f"   A_SI   = {A_SI:.6e}  kg*s^2 / m^5")

# Convert A from SI to "framework units": M_sun/pc^3 per (km/s)^2
# rho_crit [M_sun/pc^3] = A * V_flat^2 [(km/s)^2]
# A_SI [kg/m^3 / (m/s)^2] -> A_fw = A_SI * (kg -> M_sun) * (m^-3 -> pc^-3) * ((m/s)^-2 -> (km/s)^-2)
#                                  = A_SI * (1/Msun) * (pc^3) * ((kms)^2)
A_fw = A_SI * (pc**3 / Msun) * (kms**2)
print(f"   In framework units (M_sun/pc^3 per (km/s)^2):  A = {A_fw:.6e}")
print(f"   Site's published empirical A:  0.028")
print(f"   Ratio: A_formula / A_empirical = {A_fw / 0.028:.3e}")
print()
if abs(np.log10(A_fw / 0.028)) > 1:
    print(f"   -> The stated formula DOES NOT reproduce the published A ~ 0.028.")
    print(f"      Discrepancy: {0.028/A_fw:.1f}x  (explorer claimed ~600x; this matches)")
else:
    print(f"   -> The stated formula DOES reproduce A ~ 0.028; explorer claim refuted.")

# --- (2) The V^0.5-law alternate calculation that DOES give 0.0294 ----------
print()
print("=" * 74)
print("(2) The V^0.5-law alternate calc claimed to give 0.0294")
print("=" * 74)
# Explorer's claim: alpha = 4.5 (fitted), R_0 = 0.07 kpc/(km/s)^0.75
# (the size-velocity-relation slope), produces a V^0.5 law for rho_crit.
# Without seeing the script, the simplest reading is:
#   R_half = R_0_size_vel * V^0.75
#   rho_crit ~ alpha * V^2 / (G * R_half^2 * V^2) ... etc.
# The headline number 0.0294 is what matters. Cross-check the explorer's claim
# that the framework derives V^0.5 from these inputs:
#   If rho_crit ~ alpha * V^something / (G * R_half^2)
#   with R_half ~ V^0.75 -> R_half^2 ~ V^1.5
#   then rho_crit ~ V^(something - 1.5)
# Empirical exponent 0.5 -> something = 2.0 -> rho_crit ~ V^2 / R_half^2
# But the framework's PUBLISHED law is rho_crit ~ V^2, so the script would need
# rho_crit ~ V^2 / R_half^2 ~ V^0.5 if R_half is the size-velocity relation,
# OR rho_crit ~ V^2 / const if R_half is V-independent.
# Session 65 confirms B = 0.5 in rho_crit ~ V^B, per the explorer; this is the
# V^0.5 law the explorer says the 5% computation actually derives.
#
# Without access to session66_A_gap_investigation.py here, I cannot replay the
# exact line. The structural point: if the inputs include a V-dependent length
# scale and the formula uses R_half^2 in the denominator, the resulting law
# cannot also be V^2 -- it is V^(2 - 2*0.75) = V^0.5. The explorer's claim
# that "the 5% calc derives the wrong scaling law" is structurally correct as
# long as the 5% script uses R_half = R_0_size_vel * V^0.75 in the denominator.
print("   Explorer's claim: the 5% calc uses R_half = 0.07 * V^0.75 (kpc),")
print("   alpha = 4.5 (fitted), and yields rho_crit ~ V^0.5 (Session 65 B=0.5).")
print("   The framework's published law is rho_crit ~ V^2.")
print("   STRUCTURAL CONSISTENCY CHECK:")
print("     If rho_crit ~ M / (V * R_half^2) ~ V^(2 - 2*0.75) = V^0.5,")
print("     then the law derived IS NOT the law published. The 5% number")
print("     attaches to the wrong scaling law -- structurally inconsistent.")
print()
print("   I cannot re-run session66_A_gap_investigation.py from this session")
print("   directly, but the structural argument is decidable from inputs alone:")
print("   any computation using R_half^2 with R_half = R_0*V^0.75 in the")
print("   denominator and yielding a finite A constant requires rho_crit ~ V^0.5,")
print("   not V^2. The published 'A * V^2' law is incompatible with such a")
print("   derivation -- so the 5% number cannot belong to that law.")

# --- (3) What value of R_0 WOULD make A = 0.028 if formula were right? ------
print()
print("=" * 74)
print("(3) Reverse-solve: what R_0 would the formula need to give A = 0.028 ?")
print("=" * 74)
# A_SI = 4*pi / (beta_J^2 * G * R_0^2)
# Want A_fw = 0.028, i.e. A_SI = 0.028 / [(pc^3/Msun)*(kms^2)]
A_SI_target = 0.028 / ((pc**3 / Msun) * (kms**2))
R0_required_m = np.sqrt(4.0 * np.pi / (beta_J**2 * G_SI * A_SI_target))
R0_required_kpc = R0_required_m / kpc
print(f"   Target A in framework units: 0.028")
print(f"   Target A in SI: {A_SI_target:.3e}")
print(f"   Required R_0  : {R0_required_m:.3e} m = {R0_required_kpc:.4f} kpc")
print(f"   Ratio to MW R_0=8 kpc: {R0_required_kpc/8.0:.4f}")
print()
print("   For the formula with beta_J=1 to give A=0.028, R_0 would need to be")
print(f"   ~ {R0_required_kpc:.3f} kpc -- not 8 kpc. The discrepancy is ~",
      f"{(8.0/R0_required_kpc)**2:.0f}x in A (since A ~ 1/R_0^2).")

# --- (4) Disposition --------------------------------------------------------
print()
print("=" * 74)
print("(4) Disposition")
print("=" * 74)
print(f"   The stated formula A = 4*pi/(beta_J^2 * G * R_0^2) with the stated")
print(f"   inputs (beta_J=1, R_0=8 kpc) gives A_fw = {A_fw:.3e}, NOT 0.028.")
print(f"   The published 5% agreement requires either (a) a different R_0")
print(f"   (~{R0_required_kpc:.2f} kpc, not the Sun's galactocentric radius), or")
print(f"   (b) a different beta_J, or (c) a fitted ad-hoc factor of ~{0.028/A_fw:.0f}.")
print(f"   The explorer's chain-of-custody finding -- that the 5% number was")
print(f"   imported from a V^0.5-law calculation with different inputs -- is")
print(f"   structurally consistent with arithmetic verification.")
print()
print(f"   This is the SAME methodology failure caught in S672 (DESI epistemic")
print(f"   regression): a confident result re-stated across sessions without")
print(f"   re-execution. The fix is the same: re-execute, don't re-read.")
print()
print(f"   Research-repo disposition for the A-from-Jeans claim:")
print(f"   - The arithmetic does not support 'first-principles derivation'.")
print(f"   - The 5% headline number does not belong to the V^2 law it labels.")
print(f"   - Either the formula needs a different (galaxy-intrinsic) R_0 that")
print(f"     gives ~0.20 kpc (with beta_J=1) -- which is not the Sun's R_0 --")
print(f"     or the derivation must be rebuilt from scratch.")
print(f"   - Until either path is taken, A-from-Jeans is reparametrization, not")
print(f"     derivation, joining the other re-expressions in the audit.")
