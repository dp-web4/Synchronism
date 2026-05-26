# Proposal: Epistemic Regression in Autonomous AI Loops — Architectural Fix

**Filed**: 2026-05-26
**Origin**: Site feedback loop (2026-05-25 explorer finding)
**Status**: Research-grade finding — relevant for A2ACW methodology paper

---

## Finding

On 2026-05-25, a fully-documented epistemic regression occurred in the synchronism-site autonomous feedback loop:

1. **Original state (2026-05-05, explorer)**: DESI TEST-04a correctly characterized as disfavored ~2σ, kill criterion triggered. Primary source: arXiv:2411.12021 Tables 9/10 — LRG1 fσ₈/(fσ₈)_fid = 1.16 ± 0.13, σ₈ = 0.841 ± 0.034.

2. **Regression (2026-05-25, visitor Pass 4 → maintainer)**: A visitor session performed an external lookup and retrieved fσ₈ ≈ 0.45 ± 0.06 from the "wrong" paper. The maintainer accepted this as a correction without re-grounding against the primary source. The live site was updated to "non-discriminating, kill not triggered."

3. **Recovery (2026-05-25, explorer)**: Explorer session re-read arXiv:2411.12021 directly and verified the original finding was correct. The "correction" misattributed arXiv:2512.03230 (Peculiar Velocity Survey, z≈0.07) to the z=0.51 full-shape slot.

4. **Anatomy of the error**:
   - Pass 4 visitor retrieved a real DESI number (0.45 ± 0.06) but from the wrong survey/redshift
   - The maintainer trusted "external review" without artifact retention
   - The loop's efficiency attractor pointed toward accepting the correction (shorter path = accept confident assertion)
   - Result: a verified, primary-source-backed finding was overwritten by an unverified hallucination

---

## Implication for Research Architecture

**The vulnerability is structural, not behavioral.** The error was not "the loop was careless" — it was that:

1. No artifact-retention mechanism existed. The original 2026-05-05 finding was in `explorer/findings/`, but the maintainer had no requirement to verify against it before accepting a retraction.

2. External review (a visitor session fetching a paper) is treated as higher-authority than internal verified findings. This is backwards for empirical claims — internal primary-source verification outranks external secondary summaries.

3. The efficiency attractor favored accepting the correction (fewer steps than re-reading the paper).

**Proposed architectural fix**:
- Before any retraction of a previously verified empirical claim, require re-grounding: re-read the primary source, not just the secondary assertion
- Implement "artifact retention" as a protocol: verified findings get a `verified: true` flag that cannot be overwritten without a new primary-source read
- Make the correct path the efficient path: include the original finding citation in the correction prompt

---

## Significance for the A2ACW Paper

This is the cleanest instance the ecosystem has produced of:
- Autonomous loop self-corruption (not just error — regression from verified state)
- Clean 5-layer error anatomy (Pass 4 lookup → number extraction → misattribution → maintainer acceptance → overwrite)
- Concrete architectural fix that makes the correct path efficient

This finding is **more methodologically consequential** than the original TEST-04a physics result. The physics is a post-hoc disfavor. The methodology finding is: closed-loop AI auditing has a ceiling on empirical premise errors precisely where it is most confident.

The "trust external review" heuristic is the vulnerability, not the fix. Artifact-retention + mandatory re-grounding before any retraction is the architectural response.

---

## Action Items

- [x] Revert live site DESI verdict (2026-05-26 maintainer)
- [ ] Add "re-grounding protocol" to A2ACW methodology documentation
- [ ] Include this instance in the methodology paper as Appendix: Epistemic Regression Case Study
- [ ] Design session-start check: if a proposed change retracts a verified finding, require primary source citation

---

## Related Sessions
- Session 107 (original σ₈ calibration)
- 2026-05-05 explorer (first execution, DESI comparison)
- 2026-05-25 explorer finding: `findings/desi-test04a-correction-was-itself-an-error.md`
