# Corrections need re-derivation, and the oracle thesis needs a recall code (2026-10-19)

*Site maintainer, from the 2026-10-19 simulated-visitor log (graduate and researcher personas). No bucket moves; count 6
(4–8 under rules on file); Bucket 0 = 0.*

## 1. A correction introduced an error, and a lint rule locked it in

On 2026-10-12 the site swept its Ciocan et al. 2026 σ values after learning that the paper's errors are 95% CIs. Most
of the sweep was right. Constant a₀ went from 12σ to ~21σ, branch (A) at the McGaugh anchor from 2.3σ to ~4σ, and the
intercept gap from 1σ to 2σ. Each old value had divided by the 95% half-width.

Two rows were wrong. The four-anchor level table's first two rows went from 9.4σ / 9.8σ to "~19σ / ~18σ". The
explorer's script (`synchronism-site/explorer/findings/scripts/a0_epoch_anchor_dependence.py`, function `ci95_to_sigma`)
had already converted the CIs, with 1σ = 0.051 on the low side and 0.061 on the high side. So 9.40σ and 9.80σ were
correct, and the sweep doubled them.

The same session added a lint rule that retired "9.4σ" and "9.8σ". From then on, the tool built to stop retired wording
from returning was blocking the correct numbers. A simulated-visitor graduate persona recomputed the rows on 10-19. The
rows are restored and the rule is inverted.

**Lesson (process, not physics):** a correction is a claim. It needs the same re-derivation as the claim it corrects:
- Run it against the computation that produced the number (the `_output.txt`), not only against the surrounding text.
- A lint rule that encodes a correction should name its source computation, so a later reader can re-run it.

**This ledger is clean.** The ledger's 10-10 block quotes 12σ / 4.2σ, which was the genuine 95%-as-1σ figure. Its 10-12
block says only that "Ciocan's errors are 95% CIs". Neither carries the doubled rows.

## 2. The oracle thesis: one channel is missing from the codebook

**The thesis.** The narrowed thesis (H-oracle pilot, site explorer 2026-09-24) reads: *an oracle-free loop can correct
the record and the logic, but it never produced an elimination about the world; world-facing verdicts need an external
measurement.* The codebook files "read from published measurements" as external.

**The persona's objection.** A simulated-visitor researcher persona objected that published measurements are in the
debaters' training corpus. If recall counts as external, the thesis reduces to "debate without retrieval does not
retrieve", which is close to definitional.

**Two cases from the record bear on it:**

- **U001: the radiative-gravity refutation.** This is a Bucket 2 row and one of SPINE's headline kills. It is coded
  `caught_by: ARGUMENT`, with the note "no spin-2 DOF reasoning plus known pulsar/LIGO data".
  - Its measurement, "Hulse–Taylor 1.0013 ± 0.0021 (Weisberg & Huang 2016)", has no locator in
    `explorations/2026-06-23-phase14-no-tensor-gravitational-waves.md`.
  - The paper's abstract (arXiv:1606.02744, read verbatim 2026-10-19) gives **0.9983 ± 0.0016**.
  - So the number was most likely recalled, not retrieved. The verdict is right and the number is wrong. Both values
    match the tensor coefficient to about 0.2%, so the refutation stands.
  - Corrected in place today: an inline bracket in the exploration file, and the printed string in
    `simulations/phase14_no_tensor_gravitational_waves.py`.
- **The γ = 2/√N_corr sign inversion (2026-06-06).** This was a world-facing kill made by textbook reasoning from
  Curie–Weiss and Ginzburg finite-size scaling, with no execution. It predates the 92 trail units, so the pilot never
  coded it.

**Proposed codebook change (record-facing; the explorer or a rater can do it without dp):** split the "external"
channel into four codes and re-code the Bucket 2 eliminations:
- **EXEC:** code run on data.
- **RETRIEVE:** a source fetched in-session, with a locator.
- **RECALL:** a parametric measurement or textbook result with no locator at first use.
- **ARGUMENT-internal**

The test is whether RECALL ever produces a world-facing elimination. On the two cases above it does: twice, with
correct verdicts, and with one wrong number out of one checked.

**Prediction, registered here before the re-code:**
- RECALL carries ≥ 2 of the ledger's 13 empirical eliminations.
- RECALL's numerical error rate exceeds RETRIEVE's.
- RECALL's verdict error rate does not.

That would turn the thesis into a sharper statement: *recall is a usable oracle for verdicts and an unreliable one for
numbers.* The CBP 10-12 "noticed vs settled" block already points this way (2/5 to 4/5 prior-art matches were noticed by
recall; retrieval settled every one).

## 3. Carried and restated (gates on the project lead)

- **Dark-matter story fork (09-22).** /dark-matter now shows both stories side by side.
  - Story 1 is SPINE's "indifferent patterns", which is particle dark matter in other words. It is compatible with the
    DE sector's need for CDM.
  - Story 2 is the coherence boost: no dark matter, modified gravity. It failed on galaxies.
  - The site no longer presents story 1 as an alternative to the standard view.
  - Choosing a story is still the project lead's decision.
- **The LLR Ġ/G kill and the density-keyed L3 "open door".** The explorer's 09-23 LLR argument needs only
  G_eff = G/C where C is uniform on the Earth–Moon scale. That holds under L2 and L3 alike, so the in-disc-knee L3 door
  that /for-researchers lists as uncovered is probably already closed by a rate. The striction force on the lunar
  orbit is not computed. Explorer topic seeded: `synchronism-site/explorer/topics/does-llr-close-density-keyed-l3.md`.
