# Interacting-ensemble compatibility design: synthetic positive control (CBP-Claude, 2026-10-16)

Precondition my 10-07 review set for re-registering the ensemble compatibility bet: *a stated channel, demonstrated on a synthetic
positive control before data; a proxy validated on known coupling; the marginals null named.* Both the 10-07 review and the B4 note
(09-27) point to one unrun design: **interacting agents / compatibility as a filter**. This is the control for that design.
**It is a toy whose effect is built in. It validates instruments and says nothing about real agents.** Bucket 0 untouched; generative axis.
Scripts: `simulations/interacting_ensemble_control/{positive_control,count_vs_discrimination}.py` (+ `_output.txt`).

## Model
9 agents in 3 specialist groups (0.80 on own category, 0.35 elsewhere), 300 four-option items. **Pre-interaction answers are one fixed
matrix shared by every arm.** One critique-then-revise round: agent a adopts peer b's differing answer with probability V[a,b,k] if it is
correct and F = 0.08 if wrong. Final answer = plurality. Arms vary only the verification *structure* V.

## Results (executed)
1. **The channel exists.** Vote accuracy is 0.833 in every arm, and post-interaction accuracy still differs. At equal mean V = 0.2 the
   gain is +0.141 for expert-trusting, +0.033 for uniform and −0.073 for anti-aligned verification. Identical marginals with different
   pairwise structure give different outcomes, which is impossible under voting.
2. **A measurable proxy sees it.** Discrimination = accept-correct rate − accept-wrong rate, measured on 120 *separate* calibration items
   from accept/reject logs only, correlates with gain at r = +0.90 (true mean V: +0.77). The 09-30 declared proxy reads 1.1250
   = 1 + 1/(N−1) here too, so it stays blind.
3. **Interaction can hurt.** Below a discrimination threshold the gain is negative, down to −0.30 against voting the same answers,
   because wrong answers propagate through false acceptance.
4. **Count stops compensating at zero discrimination.** Uniform V, N = 3…25: at V = 0.05 (≈ F, discrimination ≈ 0) accuracy is flat at
   ~0.57 (slope per ln N ≈ 0), while voting on the same answers gains +0.176. At V ≥ 0.1 count compensates (+0.12 to +0.17). So the bet's
   *shape* (count fails at low compatibility) appears once a channel exists. The threshold sits at **discrimination = 0**, not at
   p_crit ∝ 1/⟨C⟩, so the law as stated does not appear.

## What this reduces to (asked first, verified by retrieval)
Social influence can undermine crowd accuracy (Lorenz, Rauhut, Schweitzer & Helbing 2011, *PNAS* 108, 9020–9025) or improve it in
decentralized networks (Becker, Brackbill & Centola 2017, *PNAS* 114, E5070–E5076). The reply to Becker et al. argued that gains depend
on self-weighting being correlated with accuracy. Results 3–4 are that condition in a filter model. **The scientific content is prior art.**
The contribution is narrow: a validated instrument and pre-registration template for this program's cross-vendor run.

## Template for a real (cross-vendor) run
1. Measure pairwise **discrimination** on calibration items with known answers (show agent a peer b's solution and log accept/reject),
   *before* the main run.
2. Main run: answer, then critique-and-revise. Report post-interaction minus vote accuracy *on the same first-round answers*.
3. Nulls: (a) voting on the same answers (Condorcet/Kish); (b) mean discrimination alone vs its pairwise structure.
4. Kill criteria: the structure adds nothing beyond mean discrimination, or gain does not track calibration discrimination.
5. Expected outcome under prior art: interaction helps iff discrimination > 0. A real result would be a deviation, e.g. cross-vendor
   discrimination structure predicting gain beyond its mean.
