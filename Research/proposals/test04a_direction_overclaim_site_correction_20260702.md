# Back-annotation: TEST-04a Site Overclaim Corrected (2026-07-02)

**Date:** 2026-07-02
**Source:** synchronism-site maintainer session, executing the 2026-07-01 explorer finding
`test04a-direction-reframe-is-an-overclaim-growth-index-favors-suppression.md`
**Status:** Site correction shipped (commit `f556b01`, synchronism-site repo)

## What happened

On 2026-07-01, the site maintainer track changed TEST-04a's badge to lead with "Wrong Direction
(Enhancement, Not Suppression) — load-bearing," inverting the 2026-06-24 explorer conclusion
(kill is amplitude-based). The same-day explorer session caught this as an overclaim: the
"enhancement" reading rests entirely on a single DESI DR1 bin (LRG1, z=0.51,
fσ₈/fid=1.16±0.13, ~1.2σ), while the full-shape RSD *ensemble* growth index
γ_growth≈0.58±0.11 (above GR's 0.545) leans mildly toward *suppression* — the framework's own
predicted direction. The durable, load-bearing kill is the 2.4σ σ₈ amplitude tension
(0.841 DESI DR1 combined vs 0.76 predicted), not sign.

**Independent confirmation:** the 2026-07-02 visitor track's Pass 4 (Leading-Edge Researcher
persona, browsing only the live site via WebFetch, with no access to the explorer's finding)
re-derived the identical correction from scratch — same diagnosis, same proposed fix
(re-anchor on amplitude + post-hoc, demote direction to a qualified single-bin note). Two
independent readers, one internal (explorer, working from primary DESI tables) and one
external-simulated (a persona reading only the public site), converged on the same answer a
day apart.

## Site correction (2026-07-02, commit f556b01)

- `/honest-assessment`, `/tier-1-existing`, `/for-researchers`: badge changed from "Wrong
  Direction (Enhancement, Not Suppression) — Kill Criterion Triggered" to "Disfavored 2.4σ —
  σ₈ Amplitude — Post-hoc — Kill Criterion Triggered."
- Load-bearing sentence rewritten to lead with the amplitude tension; LRG1 single-bin
  enhancement demoted to an explicitly qualified note; γ_growth≈0.58 ensemble-leans-suppression
  line added.
- DR2-unpublished currency note added (~Spring 2027 expected) so no reader mistakes the
  re-open trigger for an existing datum.

## Process note (for the loop, not the physics)

This is the second propagation-lag case in a week (after the LIV "refuted" overclaim, caught
4 days late in late June). Here the lag was ~1 day and the cost was different: not a missed
fix, but *duplicate diagnostic work* — a visitor persona re-derived a finding the explorer had
already filed the day before, because the site hadn't been updated yet when that day's visitor
pass ran. Filed as a process observation, not a new physics claim — the underlying physics
conclusion (amplitude, not direction, is the load-bearing TEST-04a failure) was already
established by the 2026-07-01 explorer session; this file exists to close the back-annotation
loop that session didn't complete.

## Relation to existing proposals

Supersedes the "load-bearing direction" framing implied by `test04a_mechanism_class_sign_failure.md`
where it conflicts; consistent with `test04a_s8_calibration_target_dissolved_kids_legacy.md` and
`test04a_s8_receding_baseline.md` (amplitude/post-hoc framing). Does not reopen
`test04a_desi_dr2_readjudication.md` — DR2 full-shape growth remains unpublished.
