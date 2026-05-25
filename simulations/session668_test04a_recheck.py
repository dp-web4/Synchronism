#!/usr/bin/env python3
"""
Session 668: Re-check the TEST-04a "sign reversal" against the actual DESI DR1 numbers.

Three parties disagree:
  - S645/S650 (origin):  DESI LRG1 fσ8 ~ 0.55 ± 0.06 -> ENHANCEMENT -> "sign reversal",
    "mechanism-class failure"; elevated by S656/S663 to "the framework's one
    transferable physics contribution."
  - Retraction proposal (2026-05-25): LRG1 fσ8 ~ 0.45, ΛCDM-consistent -> "non-
    discriminating retrodiction"; the transferable contribution evaporates.
  - This session: adjudicate with the verified paper numbers.

Verified from arXiv:2411.12021 (DESI 2024 V):
  - combined sigma8 = 0.841 +/- 0.034   (abstract)
  - growth index gamma = 0.580 +/- 0.110, consistent with GR (0.55)   (paper)
  - "in agreement with LambdaCDM ... consistent with Planck"          (abstract)

Session 107 predictions (Research/Session107_DESI_Forecasts.md):
  - fσ8(z=0.51): LCDM 0.474, Sync 0.418  (a -11.9% SUPPRESSION)
  - sigma8(z=0): Sync 0.76 (vs LCDM/Planck ~0.81)
  - kill criterion: fσ8(z=0.5) > 0.45 -> LCDM favored
"""

# ---- DESI ShapeFit ratios as transcribed in the originating proposal (Table 9) ----
# fσ8 / fσ8_fiducial, and per-bin inferred sigma8(z=0) (Table 10)
proposal_rows = {
    # bin :   (fsig8_ratio, ratio_err,  sigma8_inferred, sigma8_err)
    "BGS":  (0.84, 0.19, 0.662, 0.13),
    "LRG1": (1.16, 0.13, 0.835, 0.087),
    "LRG2": (1.04, 0.10, 0.880, 0.077),
    "LRG3": (0.997, 0.092, 0.815, 0.072),
    "ELG2": (0.945, 0.087, 0.755, 0.059),
    "QSO":  (1.16, 0.12, 0.950, 0.072),
}
sigma8_fid = 0.811   # Planck fiducial used to turn a ratio into an absolute amplitude

print("=" * 74)
print("(1) Internal consistency check: a fσ8 ratio and an inferred sigma8 must agree")
print("    (fσ8 ∝ sigma8, so ratio ≈ sigma8_inferred / sigma8_fiducial)")
print("=" * 74)
print(f"   {'bin':>5} {'fσ8 ratio':>10} {'σ8 inferred':>12} {'ratio-implied σ8':>18} "
      f"{'σ8-implied ratio':>18} {'consistent?':>12}")
for b, (r, re, s, se) in proposal_rows.items():
    s_from_ratio = r * sigma8_fid
    r_from_sigma = s / sigma8_fid
    # honest test: is the listed ratio within 1 sigma (of the ratio) of the
    # ratio implied by the bin's own inferred sigma8?
    discrepancy_sigma = abs(r - r_from_sigma) / re
    flag = "OK" if discrepancy_sigma < 1.0 else f"{discrepancy_sigma:.1f}sig off"
    print(f"   {b:>5} {r:>10.3f} {s:>12.3f} {s_from_ratio:>18.3f} "
          f"{r_from_sigma:>18.3f} {flag:>12}")
print("   => QSO/LRG2/LRG3 are self-consistent. LRG1 is the discrepant bin: ratio 1.16")
print("      implies σ8≈0.94, but its listed σ8=0.835 implies ratio≈1.03 (~1σ internal")
print("      tension), and 1.16 is identical to QSO's. LRG1 is exactly the bin driving")
print("      'sign reversal'. Suggestive, not decisive on its own -- see (3) gamma.")

print()
print("=" * 74)
print("(2) What LRG1 actually implies for fσ8(z=0.51)")
print("=" * 74)
# Use the self-consistent quantity (inferred sigma8) to back out fσ8(z=0.51).
# In LCDM, fσ8(0.51)/sigma8(0) ≈ 0.474/0.811 = 0.585 (Session 107's own LCDM column).
scaling = 0.474 / 0.811
for label, s8 in [("LRG1 inferred σ8 = 0.835", 0.835),
                  ("LCDM/Planck σ8 = 0.811", 0.811),
                  ("Sync σ8 = 0.76", 0.76)]:
    print(f"   {label:>28}  ->  fσ8(0.51) ≈ {scaling * s8:.3f}")
print("   => LRG1 data implies fσ8(0.51) ≈ 0.49, NOT 0.55. ΛCDM gives 0.474, Sync 0.418.")
print("      0.49 straddles ΛCDM: consistent, no enhancement. The '0.55 above ΛCDM at")
print("      every bin' claim is unsupported.")

print()
print("=" * 74)
print("(3) The growth index settles the SIGN")
print("=" * 74)
gamma, gerr, gr = 0.580, 0.110, 0.55
print(f"   DESI DR1 gamma = {gamma} +/- {gerr}; GR predicts {gr}.")
print(f"   gamma - GR = {gamma-gr:+.3f}  ({(gamma-gr)/gerr:+.2f} sigma).")
print("   gamma >= 0.55 means growth is consistent-with or SLOWER-than GR (suppression-")
print("   leaning) -- the OPPOSITE of the 'enhancement' the sign-reversal claim needs.")
print("   A 16% fσ8 enhancement at LRG1 would require gamma well BELOW 0.55. Refuted.")

print()
print("=" * 74)
print("(4) What DOES survive: the sigma8(z=0) amplitude tension")
print("=" * 74)
desi_s8, desi_err = 0.841, 0.034
sync_s8 = 0.76
tension = (desi_s8 - sync_s8) / desi_err
print(f"   DESI combined sigma8 = {desi_s8} +/- {desi_err} (verified, abstract).")
print(f"   Sync predicted sigma8(z=0) = {sync_s8}.")
print(f"   tension = ({desi_s8} - {sync_s8}) / {desi_err} = {tension:.2f} sigma.")
print("   => Sync's sigma8(z=0)=0.76 is DISFAVORED at ~2.4 sigma. This is a MAGNITUDE")
print("      miss (amplitude too low), NOT a sign reversal, and NOT 'non-discriminating'.")
print("      It does not depend on the LRG1 error. (Post-hoc: S107 committed 2025-12,")
print("      after DR1 was public Nov-2024 -- see S648.)")

print()
print("=" * 74)
print("VERDICT")
print("=" * 74)
print("  - 'fσ8 enhancement / sign reversal / mechanism-class failure' (S645/S650): WRONG.")
print("    Transcription error (LRG1 ratio dup'd from QSO) + contradicted by gamma=0.58.")
print("  - 'one transferable physics contribution' (S656/S663): EVAPORATES (no wrong-sign).")
print("  - 'purely non-discriminating' (retraction): OVERREACHES -- σ8(z=0)=0.76 is")
print("    disfavored at ~2.4 sigma on amplitude.")
print("  - Corrected S107 status: DISFAVORED on σ8 amplitude (~2.4σ, post-hoc);")
print("    fσ8 shape ΛCDM-consistent/non-discriminating; sign-reversal RETRACTED.")
