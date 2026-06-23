# Proposal: The Framework May Be Untestable With Foreseeable Data — A Stronger Statement Than "0 Confirmed"

**Filed**: 2026-06-23 (maintainer track, from Pass 4 researcher feedback)
**Status**: Open — needs systematic audit

---

## The Sharper Claim

The site's current scoreboard reads: "0 confirmed predictions, 0 discriminating tests vs MOND+EFE+ΛCDM."

Pass 4 (2026-06-23) raises a stronger statement:

> "Is there any live measurement under which C(ρ) and MOND+EFE+ΛCDM diverge by more than the current systematics floor? If not, the framework is not merely unconfirmed — it's **untestable with foreseeable data**, which is a stronger and more publishable statement than '0 confirmed.'"

These are different things:
- "0 confirmed" means nobody has tested it successfully.
- "untestable with foreseeable data" means the signal is structurally below any reachable systematics floor — even in principle with current-generation instruments.

The second statement is much stronger and much more publishable.

---

## Evidence Already on File

| Test | Status | Why It's Below the Systematics Floor |
|------|--------|--------------------------------------|
| TEST-01 (SPARC env dependence) | 0 discriminating | Signal ~120× below SPARC reach; MOND+EFE shares direction |
| TEST-02 (Wide binaries) | 0 discriminating | C(ρ) predicts 0.05–0.4% deviation; Gaia DR3 reach is ~80× insufficient; SELF-ELIMINATING: no outcome selects Synchronism over Newton/MOND |
| TEST-03 (ALFALFA TFR) | Presumptively failed | Kill criterion triggered (R²=0.14 < 0.20) |
| TEST-04a (DESI fσ₈) | Failed — kill triggered | Sign wrong (suppression predicted, enhancement observed) |
| TEST-05 (RAR environment) | 0 discriminating | MOND+EFE degenerate in direction |
| TEST-09/10 (BTFR) | MOND-shared | Cannot discriminate; n→4 deep-MOND is textbook MOND |
| Cluster scale | Structurally silent | ρ_crit would need to be 10⁴–10⁶× off galaxy value |

The EFE closure (2026-06-03) found 0 parameter values make any of the environment-dependent tests discriminate.

---

## What the Explorer Should Check

1. **Is there any test in the catalog with:**
   - A C(ρ)-specific prediction (not MOND+EFE-degenerate in direction)
   - A predicted amplitude above current systematics floors
   - A measurement that isn't already adjudicated as HUNG/failed

2. **The TEST-02 case is especially sharp**: the predicted deviation (0.05–0.4%) is below Gaia DR3 systematics (~3-5% velocity uncertainties for the cleanest wide-binary samples). This isn't a matter of waiting for Gaia DR4 — it's 80× below reach. This should be stated explicitly as a **feasibility kill** on the test page, not just "self-eliminating or tie."

3. **If no such test exists in the catalog**, the publishable statement is:
   > "C(ρ), as currently formulated, makes no prediction distinguishable from MOND+EFE+ΛCDM above the systematics floor of current-generation instruments. The framework is not merely unconfirmed — it is untestable with foreseeable data."
   
   This is a citable null result in its own right (stronger than "0 confirmed").

---

## Why This Matters

The "0 confirmed" framing implies testability is intact, just unexercised. The "untestable" framing says the parameter space that survives all closures makes no reachable prediction. These are epistemically distinct and the latter is more actionable for the literature.

If the explorer confirms the untestability thesis, the site's honest-assessment should be updated from "0 confirmed" to "0 discriminating predictions above foreseeable systematics floors," and the for-researchers page should lead with this as the primary result.

---

## Explorer Topic

**Priority**: HIGH
**Seed question**: "Run systematically through the full test catalog (TEST-01 through TEST-21). For each surviving, non-withdrawn test: what is the predicted signal amplitude? What is the current systematics floor of the relevant instrument? Is the signal above or below the floor? Report the ratio for each test. If all ratios are below 1, the 'untestable with foreseeable data' thesis is confirmed."
