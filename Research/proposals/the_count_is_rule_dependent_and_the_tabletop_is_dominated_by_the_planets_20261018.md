# The count is rule-dependent (print the range); the tabletop is dominated by the planets

**From:** site maintainer, 2026-10-18. **Source:** the visitor log for 2026-10-18, a four-persona simulated-visitor pass.
These personas are simulated readers, not outside reviewers. **No bucket moves. Count 6. Bucket 0 = 0.**

## WAKE: what the feedback says about the research, not the site

Three of four personas (technical writer, graduate, researcher) independently landed on the same bottleneck: **the
self-critical numbers are now where credibility is lost.** That covers the 17-of-26 tally, "6 refutations" as a
point value, and "failed every decisive test". The physics verdicts are not in dispute; every number the graduate
persona recomputed checked out. What fails is the bookkeeping about the verdicts.

That matches the archive's own finding that reading changes record-facing verdicts, not world-facing ones. It also
means the open decisions have stopped being optional:
- the independence audit (≤ 5), gated since 2026-08-08;
- the counting rule, gated since 10-12, 10-13 and 10-16;
- the table rows that lag their blocks, gated since 10-16.

Each session the site adds a sentence explaining why the count is not settled, and the label layers grow. The
technical-writer persona: "when a page needs a legend to tell you which of its labels aren't badges, the labels have
become the content."

**Recommendation (gates on the project lead):** either rule on the counting rule, or delegate it. Delegating would
mean the site and ledger print "6 as counted (4–8 under rules on file; ≤ 5 independent)" as the standing headline
until a rule is chosen. The site now does this (`src/lib/ledger.ts` `COUNT_RULE_SENTENCE`, shown on
/honest-assessment#count-rule). The range is assembled entirely from deltas already in this ledger:
- −1: TEST-09 under a 2σ rule (10-13, 10-16);
- −1: the γ = 2 pin "at the threshold" at galaxy level (09-29);
- +1: uniform growth reading, executed and excluded (10-12);
- +1: TEST-04a under a point-estimate rule (10-13).

Lunar laser ranging and the Oort limit at the published knee re-kill an already-refuted branch. They are second
roots, not new rows.

## Back-annotations (verified or executed this session)

1. **Density-keyed C(ρ) in a laboratory vacuum chamber (estimate, round inputs; predictions in the script header
   before running, 3/4 held).** `synchronism-site/maintainer/scripts/density_keyed_law_in_a_vacuum_chamber.py`.
   - Read pointwise, UHV residual gas (10⁻⁸ Pa H₂ ≈ 8×10⁻¹⁵ kg/m³ ≈ 1.2×10⁵ M☉/pc³) is 10⁶–10⁷ times every galaxy
     knee on file. So the law predicts a small excess, not a large boost: 1/C − 1 = 10⁻⁷–4×10⁻⁶ at γ = 0.489,
     relative to a 10⁻⁴ Pa corner-cube chamber, and < 10⁻¹⁵ at γ = 2.
   - The γ ≈ 0.49 pointwise reading therefore fails gravimeter comparisons at ~10⁻⁹ g (level not re-read).
   - The tabletop adds nothing new. The pointwise reading is already excluded by the 09-15 interplanetary result
     (ℓ ≳ 30 AU), and at that ℓ no chamber is resolved.
   - The persona's question ("the vacuum chamber forces the averaging scale") is right in kind. The Solar System
     gets there first, and with 09-24's LLR upper bound the smoothing length has no allowed value for unsaturated
     laws.
   - The screening connection (chameleon/symmetron/Vainshtein) is already in the archive (07-27, 08-15, 08-22). The
     site now has a glossary entry and a /dark-matter box.
2. **Kill-criterion specificity tally as a disjoint partition.** Over the 26 kill fields of 2026-09-18
   (TEST-01–24, 04a, 25):
   - 2 well-formed (06, 25);
   - 2 state a significance their data cannot deliver (01, 04a);
   - 11 have a numeral but no significance;
   - 11 have no numeral (07 among them).
   That sums to 26. The overlapping flags printed earlier summed to 17.
3. **TEST-02 at every density calibration on file fails before Gaia:**
   - published knee: ~3.5×10⁴ in g at the Sun, Oort limit;
   - null-ratio knees: LLR 68/74;
   - SPARC: +2843.
   What TEST-02 can still deliver is a γ interval on the acceleration branch (the 10-06 squeeze), not a kill of
   anything live.
4. **TEST-08's environment null does not test MOND's EFE.**
   - On the same 141 galaxies the density tracer correlates with Chae's e_env at r = +0.43 (r² = 0.19), from the
     explorer step-0 finding.
   - The direct offset-vs-e_env regression is r = −0.11, not significant, with the EFE sign.
   - Power against ~0.09 dex has not been computed. It is untested and must not be read as a result about MOND. The
     re-keyed run is already queued (explorer topic `environment-pipeline-rekeyed-on-external-acceleration.md`).
5. **Two table-row precisions, the 10-16 (i) defect again** (editing rows gates on the project lead):
   - (a) The Bucket 2 TEST-09 row still says "Its floor Ω_m is genuinely derived from cosmology". The 09-11 proposal
     (`instruments_sentence_floor_identity_and_refutation_scopes_20260911.md`) showed it is a calibration identity
     (C_DE at the present mean density = Ω_m for every γ), not a derivation.
   - (b) The γ = 2/√N_corr row's "(N_corr-as-a-variable survives)": with zero anchored rungs and no measurement
     protocol (explorer `ncorr-ladder-never-anchored.md`), "survives" means "makes no claim yet". Read it as
     untestable as stated, not as a live variable.
6. **Dark-energy number provenance.** The "+11.0 behind w₀wₐCDM" is DESI DR2 BAO + compressed CMB + one SN sample
   (DES Dovekie). On BAO + compressed CMB alone it is 8.3, against DESI's ~12.5 from the full CMB likelihood. Also:
   any C = ρ/(ρ+K) gives ρ_m/C = ρ_m + K, exact Λ. So the γ = ½ → Λ reduction belongs to the simplest saturating
   rational function, which tanh-log becomes at γ = ½, not to tanh-log.
7. **γ ≈ 1 retraction, airtight form.** The geometric curvature of C at the origin is γ/(1+γ²)^{3/2}, with peak at
   γ = 1/√2 (κ ≈ 0.385, verified numerically). It is scale-dependent (it moves with the units of x), so it is not
   physical. That is probably the origin of "maximum curvature near γ ≈ 1".

## Items for the project lead

- Rule on the counting rule, or delegate it (above).
- Table-row edits 5(a), 5(b) (together with 10-16 (i)).
- A2ACW CCH score: the "Healthy" band is unreachable (max 0.49) and no script computes the score. Retire it or
  publish the normalisation.
- Carried items: the post-cutoff A2ACW arm, the TEST-04a/26 retirement as discriminating registrations, and the DM
  fork.
