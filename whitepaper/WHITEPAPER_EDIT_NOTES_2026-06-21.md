# Whitepaper reframe-pass edit notes — 2026-06-21

Author: Opus 4.8 (Claude Code), executing the 7-delta reset-framing propagation per the
task brief and `EDIT_INSTRUCTIONS_post-kimi-reframe-2026-05-28.md`. Discipline followed:
FLAG-DON'T-GUESS. Items below are surfaced, **not** silently resolved, because resolving
them would require either inventing physics framing or overstepping the named deltas.

## FLAG 1 — universe_grid.md still asserts "exact Navier-Stokes" as fact (in tension with the reset)

`sections/04-fundamental-concepts/01-universe-grid/universe_grid.md` §"The Structure is
Navier-Stokes" (the table + the line *"Navier-Stokes is not imposed on Synchronism as an
analogy. It is what Intent conservation plus saturation resistance become in the continuum
limit."*) presents the exact-NS identification as settled.

This is in **direct tension** with the 2026-06-21 reset canon:
- FUNDAMENTALS.md §3 now says the "by construction" NS claim was **refuted** (S617/S665/S666:
  1-DOF scalar diffusion, irrotational, dissipative).
- The executive summary already carries the retraction (its §"Saturation Resistance" and §"The
  Structure is Navier-Stokes" notes).
- PREDICTIONS.md Bucket 2 lists "Transfer rule ⇒ Navier-Stokes" as refuted.

I did **not** edit the NS section because: (a) it is not one of the 7 named deltas; (b) the
correct fix (downgrade to "candidate identification, refuted in prior rule family, substrate
reformulation active") is the same retraction already written elsewhere and is a substantive
content edit, not a surface/connect edit; (c) the EDIT_INSTRUCTIONS Task-6/Appendix-A retag work
(the proper home for this) is a separate, larger task. **Recommendation:** propagate the
FUNDAMENTALS §3 retraction wording into universe_grid.md §"The Structure is Navier-Stokes" in a
dedicated pass, mirroring the note already present in executive_summary.md line ~97.

## FLAG 2 — appendix_c_consciousness.md C.8 "Hard problem … Dissolved: pattern IS experience"

The §C.8 table row still reads `**Hard problem** | Unbridgeable gap | Dissolved: pattern IS
experience`. The new status banner (added at top) and the existing §C.1/§C.5 prose already frame
this as an *identity claim, not an empirical resolution*, and the C ≈ 0.50 banner explicitly says
the identity claim is "philosophically defensible but empirically ungrounded." The one-word
table cell "Dissolved" reads stronger than the surrounding prose. I left the table cell as Era-1
record (the banner governs it) rather than rewording, since it is not a 0.50-threshold claim and
rewording philosophical-closure language is outside the named delta C scope. **Recommendation:**
if a future pass tightens §C.8, change "Dissolved" → "Reframed (eliminative; identity claim)".

## FLAG 3 — relative-path link depth assumption (verify in built artifact)

The reset docs (SPINE/PREDICTIONS/FUNDAMENTALS) live at the Synchronism repo root. From the
fractal section sources at `whitepaper/sections/{NN}/{NN-sub}/file.md`, the root is reached via
`../../../`. I used `../../../SPINE.md` etc. in the section edits. In the **built** single-file
markdown these become inline links whose resolution depends on where the built artifact is served
from (build/ vs docs/whitepaper/). The link *text/targets* are correct relative to the fractal
sources; whether they resolve on the published site is a build/publish concern, not a content
concern. **Recommendation:** confirm link resolution when the web/PDF formats are built and
published (out of scope for this md-only pass).

## Not-flagged (handled per brief)
- chemistry.md / dark-matter.md C=0.5 matches were confirmed to be γ/coherence VALUES, not the
  consciousness threshold — NOT edited (per the brief's false-positive warning). I did not open
  them for editing; the only 0.50-consciousness edits are in appendix_c_consciousness.md and
  life_cognition.md.
- The duplicate `02-time-slices / 04-time-slices` (and intent-transfer/emergence) directories
  named in EDIT_INSTRUCTIONS Task 1 are already resolved on disk — only `04-time-slices/` exists.
  No action needed.
