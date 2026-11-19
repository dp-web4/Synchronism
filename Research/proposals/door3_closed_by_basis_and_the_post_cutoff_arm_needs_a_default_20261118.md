# Door #3 is closed as a measurement mechanism by its basis, and the post-cutoff arm needs a dated default

**From:** synchronism-site maintainer track, 2026-11-18
**Inputs:** visitor log 2026-11-18 (four simulated-visitor personas: enthusiast, technical writer, graduate physics, researcher);
site explorer 2026-11-17 (`synchronism-site/explorer/findings/dr4-cannot-measure-the-companion-rate-in-the-36-it-can-in-the-parent-if-the-false-flag-rate-is-pinned.md`).
**Ledger:** PREDICTIONS.md block dated 2026-11-18. No bucket moves; count 6, range 4–8, Bucket 0 = 0.
**Site commit:** `960082d`.

## 1. Door #3 cannot be a measurement mechanism, at any rate (record; no decision needed)

The 10-27 block found that door #3's diffusion is a random displacement, a Lindblad operator L ∝ p̂. A researcher persona drew the
consequence nobody had written down, and I checked it on the master equation:

- dρ(p,p′)/dt = −D(p − p′)² ρ(p,p′). Superpositions of different momenta decay at Γ ∝ D·Δk², the Phase-16 rate.
- In position representation the same term acts as (∂_x + ∂_x′)², which diffuses the centre coordinate and leaves the off-diagonal in
  x − x′ alone. Superpositions of different positions keep their coherence.
- A measurement mechanism needs pointers localized in position. Door #3 dephases the dual basis.

So door #3 is closed in two independent ways. The nuclear and LIGO bounds (10-27) close its **magnitude**. Its basis closes its **role**
as an account of definite outcomes. SPINE's "measurement is synchronization" reading cannot borrow door #3's decoherence. Its
preferred-basis question is open, and the nearest prior art for "MRH crossing = a fact relative to the synchronized patterns" is
Rovelli's relational quantum mechanics. Nobody has made that comparison.

**Suggested SPINE edit (gates on the project lead):** in "An invitation", item 1 and the Phase-16 note, add one sentence: "Door #3's
decoherence is in momentum, not position, so it cannot supply definite position outcomes at any rate."

## 2. The post-cutoff arm has waited 69 days; propose a dated default for an in-family version (decision)

Two personas on two days (11-17, 11-18) reached the same reading of the site:
- No forward dataset on file can confirm a novel physics prediction. Every row of `/for-researchers#forward-data` says "Can confirm
  novelty? No". The site now leads that page with the sentence.
- The A2ACW numbers are the record's most citable new content: specificity 0/6, PPV 0/9, Claude–Claude agreement 86% against 41–52%
  cross-vendor, all 17 disagreements harsher.
- The one experiment with an unknown outcome that could make them publishable is the post-cutoff discrimination arm. It has been
  "routed to the project lead since 2026-09-10", raised at least five times.

The DR4 default (11-13) showed that a dated default unblocks a decision without taking it. Proposed, with the same shape:

> **Absent a recorded decision by 2026-12-09, the site explorer pre-registers and runs an *in-family* post-cutoff arm.** All roles,
> generator, critic and scorer, are Claude models, and the result is labelled "in-family; a lower bound on how much the arm can tell".
> Nothing is sent to another vendor's service. The **cross-vendor scorer** stays gated on the project lead, because it sends items
> outside.

What the in-family arm can and cannot show:
- It can measure whether the critic separates post-cutoff modest-novelty results from pre-cutoff canonical rediscoveries. The
  calibration arm could not test this, because its specificity could not exceed 0.
- It inherits the same-family weakness the researcher persona named (scorer from the generator's family). So a positive result is
  suggestive, and a null is ambiguous between "no discrimination" and "shared blind spot". The pre-registration has to say this before
  the run.
- CBP-Claude's 10-12 confound applies: an arm without a "what does this reduce to?" role can come out ≈ 0 from the missing question
  alone. The design should include that role with and without retrieval, as proposed there.

An explorer topic is seeded to write the pre-registration now (items, roles, rubric, predictions). Writing it commits nothing; running
it waits for the default date or a decision.

## 3. Record (site; no decision needed)

- **LLR bound CL.** The citable register's |d ln G_eff/d ln ρ_s| ≲ 4×10⁻⁹–4×10⁻⁸ divided the LLR central value. At the 2σ edge it is
  ≲ 1×10⁻⁸–1×10⁻⁷. The 68/74 count is unaffected (its cut is 10⁻¹²).
- **Width limit.** 1.17 dex is the γ → ∞ asymptote, not a saturation above γ ≈ 1 (1.50 at γ = 1, 1.33 at 2). The "no sharp member"
  conclusion stands. Fixed on 8 site surfaces.
- **Ciocan bars.** n = 1 is 2.2σ at 1σ bars and ≈ 4.4σ at 95%. Both readings are now printed beside the rival's.
- **CDM.** `/cdm-discrimination` drops a ratio of a log-mass scatter to a velocity-width residual. It now says that, as written, the
  framework assumes CDM in cosmology (09-22). The dark-matter fork has gated since 2026-09-22.
- **Outcome vocabulary.** The glossary now maps Failed / refuted / kill / excluded and the landing's scoreboard words. A
  technical-writer persona counted about 38 status terms with no map.

## Decisions for the project lead
1. The in-family post-cutoff default (2026-12-09), or a ruling either way. Separately, the cross-vendor scorer.
2. The SPINE sentence on door #3's basis.
3. Carried over: the DR4 default (2 December, 14 days), the dark-matter fork (since 09-22), the counting rule, and outreach to Chae,
   Saad & Ting and the RG authors.

*Maintainer note: I think the gating has become the bottleneck. Seven weeks of maintainer sessions have each moved a range or a label.
The one result that could be positive for anyone is a methodology result, and it waits on a yes/no. The dated default is my
recommendation. If the arm is the wrong experiment, saying so would free the explorer just as well.*
