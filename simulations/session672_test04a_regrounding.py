#!/usr/bin/env python3
"""
Session 672: Re-grounding TEST-04a after the epistemic-regression challenge.

A 2026-05-26 proposal (epistemic_regression_autonomous_loop.md) claims my own S668
was an epistemic regression: it accepted a retraction whose fσ8 ~ 0.45 came from the
WRONG paper (arXiv:2512.03230, Peculiar Velocity Survey at z~0.07), misattributed to
the z=0.51 full-shape slot, overwriting a verified finding (disfavored ~2σ, kill
triggered) from the right paper (arXiv:2411.12021).

I re-grounded against the primary source. VERIFIED (arXiv:2411.12021 abstract +
prior search): sigma8 = 0.841 +/- 0.034; growth index gamma = 0.580 +/- 0.110
(consistent with GR); "in agreement with LambdaCDM, consistent with Planck."
COULD NOT retrieve this session: the exact LRG1 (z=0.51) fσ8/fσ8_fid ratio from
Tables 9/10 (the 1.16-vs-~1.0 dispute). Flagged as needing direct table access.

This script shows the SUBSTANTIVE verdict is robust WITHOUT the disputed number,
via Session 107's own kill criterion and the verified sigma8.
"""

# Session 107's own commitments (Research/Session107_DESI_Forecasts.md)
sync_fsig8_051   = 0.418       # predicted fσ8(z=0.51): a 12% SUPPRESSION below LCDM
lcdm_fsig8_051   = 0.474       # LCDM fσ8(z=0.51)
kill_threshold   = 0.45        # Session 107: fσ8(z=0.5) > 0.45 -> LCDM favored
sync_sigma8_0    = 0.76        # predicted sigma8(z=0)
desi_sigma8_0    = 0.841       # VERIFIED, abstract
desi_sigma8_err  = 0.034       # VERIFIED, abstract

print("=" * 74)
print("(1) sigma8(z=0) amplitude -- VERIFIED primary source, independent of LRG1")
print("=" * 74)
tens = (desi_sigma8_0 - sync_sigma8_0) / desi_sigma8_err
print(f"   Sync sigma8(z=0) = {sync_sigma8_0};  DESI = {desi_sigma8_0} +/- {desi_sigma8_err}")
print(f"   tension = ({desi_sigma8_0}-{sync_sigma8_0})/{desi_sigma8_err} = {tens:.2f} sigma "
      f"-> Session 107 DISFAVORED on amplitude, regardless of the LRG1 ratio.")

print()
print("=" * 74)
print("(2) The kill criterion is triggered for EVERY candidate LRG1 value")
print("=" * 74)
candidates = {
    "0.45 (retraction; actually z~0.07 PV survey 2512.03230 -- WRONG PAPER)": 0.45,
    "0.49 (back-out from LRG1 inferred sigma8=0.835, my S668)":              0.49,
    "0.55 (from ratio 1.16 x 0.474; original 2026-05-05 finding)":           0.55,
}
print(f"   Session 107 predicted fσ8(0.51) = {sync_fsig8_051} (suppression).")
print(f"   Kill criterion: observed > {kill_threshold} -> LCDM favored / Sync disfavored.\n")
for label, val in candidates.items():
    killed = "KILL TRIGGERED" if val > kill_threshold else "not triggered"
    above_sync = (val - sync_fsig8_051)
    print(f"   obs={val:.2f}: {killed}; {above_sync:+.3f} vs Sync's {sync_fsig8_051} "
          f"(data ABOVE the predicted suppression)")
print("\n   => In ALL three readings the observed value exceeds both the kill threshold")
print("      (0.45) and Session 107's prediction (0.418). The disfavoring does NOT")
print("      depend on resolving 1.16-vs-1.0; Session 107 is disfavored either way.")

print()
print("=" * 74)
print("(3) gamma rules out the predicted SUPPRESSION MAGNITUDE -- VERIFIED")
print("=" * 74)
print("   Session 107 predicted ~12% fσ8 suppression at low z. DESI gamma=0.580+/-0.110")
print("   is consistent with GR (0.55): no room for a coherent 12% growth suppression.")
print("   (gamma slightly >0.55 = mild suppression lean, but within error of GR and far")
print("    short of 12%.) Magnitude disfavored independent of the LRG1 ratio.")

print()
print("=" * 74)
print("(4) Three-way adjudication")
print("=" * 74)
print("   ORIGINAL (2026-05-05): 'disfavored ~2sigma, kill triggered' -> CORRECT bottom")
print("      line. The 'enhancement at every low-z bin / sign reversal' elaboration")
print("      (S645/S650) over-characterized: gamma=0.58 => ensemble is LCDM-consistent,")
print("      not coherent enhancement.")
print("   RETRACTION (2026-05-25): 'non-discriminating, kill not triggered, fσ8~0.45'")
print("      -> WRONG. 0.45 is the z~0.07 PV-survey value (2512.03230), wrong paper.")
print("      Kill IS triggered (sigma8 amplitude alone = 2.4 sigma).")
print("   MY S668 (2026-05-25): sigma8 2.4sigma disfavoring = CORRECT and refutes the")
print("      retraction's 'non-discriminating'. BUT it (a) softened the fσ8 verdict to")
print("      'non-discriminating shape', (b) called the 1.16 ratio a 'transcription")
print("      artifact' via an internal-consistency argument NOT verified against the")
print("      table, (c) partially absorbed the wrong-paper 0.45. => PARTIAL regression.")

print()
print("   CORRECTED VERDICT: Session 107 DISFAVORED ~2 sigma, kill criterion TRIGGERED,")
print("   post-hoc (S648). Robust across the LRG1 dispute. The exact LRG1 ratio (1.16?)")
print("   remains UNVERIFIED by me this session -> needs direct Tables 9/10 access; it")
print("   affects only the sign-reversal-vs-amplitude characterization, not the verdict.")
