# Session 631: BTFR Exponent and α² Audit — Site-Visitor Flag Resolved

**Date**: 2026-04-23
**Type**: Targeted archive audit responding to 2026-04-23 back-annotation
**Grade**: A- (specific, actionable, closes two claims)

---

## Trigger

A site visitor (researcher persona, Pass 4) flagged two issues in the public-facing Synchronism site:

1. **TEST-09 BTFR exponent n ≈ 2.2** contradicts Lelli, McGaugh & Schombert 2019 (n = 3.85 ± 0.09). |Δn| = 1.65, which exceeds the stated kill criterion |Δn| > 0.3 by 5.5×.
2. **A = 4π/(α²GR₀²)** uses α = fine-structure constant (1/137)? If so, it is a unique claim connecting EM coupling to galactic dynamics and warrants a Lagrangian.

The back-annotation proposal (`Research/proposals/btfr_exponent_falsification_and_alpha_coupling.md`) requested an archive search for the derivations. This session does that search and reports verbatim.

## Findings

### Claim 1: BTFR n ≈ 2.2 — sourced to Session #48, explicitly labeled "Not rigorous"

`Research/Session48_Beta_Derivation_Resolution.md` Track B contains a table comparing derivation strategies:

| Approach | Predicted n | Status |
|----------|-------------|--------|
| Virial theorem alone | Cannot derive | Needs R(v) |
| MRH boundary sets size | n = 0 | Wrong |
| Coherence threshold | **n = 3 − B/2 ≈ 2.2** | Too low |
| MOND limit | n = 4 | **Not rigorous** |

The formula n = 3 − B/2 with empirical B = 1.62 yields n ≈ 2.19. Session #48 line 91 explicitly labels the MOND n = 4 derivation "Not rigorous." The Synchronism n ≈ 2.2 is labeled "Too low" — i.e., Session #48 itself acknowledged the discrepancy with observation.

The escape hatch offered in Session #48 (lines 109–110): "The observed n ≈ 4 may include dark matter contribution. n ≈ 2.2 represents the baryonic component only." This asserts a baryonic-only interpretation without a specific-sample derivation. **The observed value n = 3.85 ± 0.09 from Lelli+2019 IS the baryonic BTFR** (it uses M_baryon = M_star + M_gas, not dynamical mass). The escape hatch does not apply. The prediction is falsified.

Session #48's own conclusion (line 114): "BTFR (M ∝ v⁴) is NOT yet derivable from Synchronism alone, but the relationship n = 3 − B/2 provides a theoretical connection." The word "yet" has been doing load-bearing work for ~500 sessions.

**Verdict**: n ≈ 2.2 is semi-empirical (n = 3 − B/2 with fitted B), and the prediction is refuted by observation. Should be moved to the honest-assessment failure catalog alongside YBCO T_c and the Bullet Cluster viscosity sign error.

### Claim 2: A = 4π/(α²GR₀²) — α is NOT the fine-structure constant

`Research/Session66_All_Parameters_Derived.md` line 61: `α = 1.0 (fiducial)`.

The formula's numerical agreement comes from α = 1.0, not α = 1/137. If α were the fine-structure constant (~7.3 × 10⁻³), then α² ≈ 5.3 × 10⁻⁵, and the formula would give A off by a factor of ~1.9 × 10⁴ from the empirical A = 0.028.

The 4π has defensible motivation in Session #66 (Jeans mass criterion + 4πR² surface area + spherical solid-angle averaging). The α² is a dimensional tuning factor with fiducial value 1, where the symbol "α" was chosen but the value 1/137 was never used.

The **public site's labeling** of this parameter as "α = fine-structure constant ≈ 1/137" in `A = 4π/(α²GR₀²)` is misleading. Two possible readings:

- If the site intends α = 1/137: numerics fail by ~4 orders of magnitude. Claim is wrong.
- If the site intends α = 1.0 (fiducial): "fine-structure constant" labeling is incorrect notation. Claim is numerology.

Either reading fails. The site needs correction.

**Verdict**: α² in A is a fitting parameter with fiducial value 1, not the QED coupling. Any derivation connecting galactic dynamics to fine-structure would require a Lagrangian, which no Synchronism session produces. The site's parameter-derivations page should be updated from "Validated/Derived" to "Dimensional Fit | 4π justified, α² is a fiducial normalization."

## What This Adds to the Demolition Arc

S617–628 focused on foundational dynamics (transfer rule, EOS, coupling, computation). S631 adds a *site-facing* falsification that was invisible from the internal-physics angle — because it concerns what the public-facing framework CLAIMS, not what the internal mathematics can DO.

This is a distinct failure mode from S617–628:
- **Internal-physics failures** (S617–628): the math doesn't do what the framework claims.
- **Public-claim failures** (S631): the framework's site presents derivations the archive shows are not rigorous, and uses notation ("α") that invites a stronger reading than the archive supports.

The Publisher flagged "synthesis-sufficient" as the current state. S631 shows the synthesis needs one more edit: the site must be aligned with the archive. The archive says the BTFR prediction is not rigorous and the α² is a tuning parameter; the site should say the same.

## Actions Taken in This Session

- Verified both claims directly against Session #48 and Session #66 source files.
- Confirmed the subagent archive-search findings.
- This document.

## Actions Recommended (Not Taken — Would Modify Site)

The operator's CLAUDE.md states `synchronism-site/` is reference-only for worker sessions. The following updates belong to the operator or a dedicated site-editing session:

1. **TEST-09 → move to honest-assessment failure catalog** with label "BTFR exponent n ≈ 2.2 (full-sample) refuted by Lelli+2019 observation n = 3.85 ± 0.09. Derivation n = 3 − B/2 is semi-empirical (Session #48 explicit acknowledgment: 'not rigorous')."
2. **A = 4π/(α²GR₀²) → relabel** on parameter-derivations page: "Dimensional fit. 4π motivated by Jeans/solid-angle; α = 1.0 is fiducial normalization (Session #66), not the fine-structure constant. Labeling α invites misreading."
3. **Visitor feedback channel worked correctly.** The Pass 4 researcher persona caught what 630 internal sessions did not — because the issue was a disconnect between public claim and archive content, not an internal-physics failure. This is a reproducible methodology: external readers auditing public claims against archive derivations.

## So What?

The demolition arc closed the internal-physics investigation. The stopping protocol held silence correctly for 12 firings. This one broke silence because genuine new content arrived (a visitor-triggered proposal), and the required action was a bounded archive audit with specific targets.

Two site claims now have clean dispositions: BTFR n ≈ 2.2 is falsified, α² is fiducial. The site needs an edit; the archive already says what it needs to say.

The productive response wasn't a new simulation or a new frame question. It was reading the two source files and reporting what they say. The surprising part: Session #48 itself labeled the derivation "not rigorous," and Session #66 itself wrote α = 1.0 (fiducial) — the falsifications were already in the archive, six hundred sessions ago. No one brought them forward to the public site.

## Files

- `Research/Session631_BTFR_and_Alpha_Audit.md` (this document)
- Updates to SESSION_FOCUS.md Validation State table
