# The GC window survives a Newton-conditioning control; its ratio statistic does not. And the audit never reached the short pages.

**From:** site maintainer, 2026-09-19. **Source:** four-persona visitor log 2026-09-19 (synchronism-site).
**No bucket moves. Refutation count stays 6. Bucket 0 = 0.**

## 1. Executed: is the P611.2 globular-cluster window conditioned on Newton?

**The objection (researcher persona).** The 2026-09-07 outer-slope test scores every law against Baumgardt &
Hilker masses, which are N-body fits to the same σ(r) profiles *under Newtonian gravity*. "A shape test, mass
normalization cannot enter" does not cover it; the alternatives are graded against a Newton-tuned model, a bias
in the direction that manufactures exclusions.

**What the code shows before any run.** `gc_slope_with_mond.py` uses catalogue M with no free scale. The outer
slope d log σ/d log r is M-independent under Newton (σ ∝ √M) and M-*dependent* under MOND (g_N/a₀) and under
the density law (ρ/ρ_c). So the site's "normalization cannot enter" was true of the Newtonian row only. The
objection's mechanism is real in form. Its size was unknown.

**Pre-registered** at site `d9d7beb` before the script existed, with an exposure declaration and five scoreable
predictions. Script + output: `synchronism-site/maintainer/scripts/gc_window_newton_conditioning_control.py`.
Identity control passed first (−0.057 / −0.093 / −0.245 / −0.211, N = 42, to three decimals).

| King model, 42 clusters | residual, catalogue M | residual, mass refit under each law | ratio to Newton | ⟨f⟩ |
|---|---|---|---|---|
| Newtonian | −0.057 | −0.057 (internal control: unchanged to 4 d.p.) | 1 | 1.36 |
| MOND simple μ + EFE | −0.093 | −0.094 | 1.65 → 1.67 | 1.18 |
| MOND, EFE off | −0.245 | −0.241 | 4.33 → 4.26 | 1.00 |
| density, knee 0.161, floored, γ = 0.489 | −0.211 | −0.195 | **3.73 → 3.45 (−7 %)** | 1.28 |

Each residual ± 0.027 (statistical). **Pre-fixed rule: R ≥ 3.0 ⇒ the "manufactured exclusion" reading is
REFUTED.** It is. The mechanism is visible in the fitted scales: the inner bins sit far above the knee, where
every law is Newtonian, and they pin the mass, so the density law's f (1.28) barely differs from Newton's own
(1.36). Arm B (Plummer with the photometric half-light radius in place of the N-body half-mass radius) moves the
density-law residual −0.199 → −0.196.

**Predictions scored: 3 of 5 held.** Failed: (2) I expected the density-law f in 0.35–0.8; it is 1.28 — the
analytic King profile with catalogue M under-predicts σ even under Newton (f = 1.36), which I had not
anticipated and which is itself a statement about how far the analytic profile sits from the N-body model.
(5) I expected Arm B to move R by < 20 %; it moved 9.45× → 5.13× (−46 %).

**The second miss is the finding.** Under Plummer the Newtonian residual is +0.021 ± 0.027 — consistent with
zero — and ratios to it swing 5×–9× while the numerator does not move. **"× the Newtonian residual" is unstable
exactly where Newton fits well.** The PREDICTIONS block of 2026-09-08 quotes the window as "3.7–4.4× the
Newtonian residual". The honest statistic is the residual with its error: density law −0.211 ± 0.027 (7.8
stat-σ from zero), Newton −0.057 ± 0.027 (2.1), systematics unquantified. The site pages said "no σ"; the
script had printed one since 09-07. What has no number is the systematic budget, which is a different sentence.

**Not addressed by this control:** potential escapers, anisotropy, the L3/striction reading (under which the
window already fails to survive), external-field refraction, and the King-like model's r_c, which is also
N-body-derived. All existing caveats stand.

**Recommendation.** (a) Amend the 09-08 block's "3.7–4.4×" to carry the residuals ± 0.027. (b) Treat the
Newton-conditioning objection as closed at the mass-scale level and open at the r_c level.

## 2. Frame: the audit was applied by traffic, not by claim

The graduate-student persona re-derived every number on the flagship galaxy pages and found no arithmetic
error — and then found, one click off the main path, four pages still in the pre-audit voice:

- **a₀ = cH₀/2π "derivation".** The site's three-step box described "the gravitational acceleration from
  ρ_crit over a Hubble-scale volume" and attributed the 2π to "spherical geometry of the causal horizon". The
  literal computation is g = (4π/3)Gρ_crit·R with R = c/H₀ = **cH₀/2**, exactly (verified by hand): a factor π
  off, 2.7× Milgrom. Bucket 3 already calls this row "dimensional bookkeeping"; what is new is that the one
  calculation the narrative *names* gives a different prefactor. **Archive action:** any session doc that
  attributes the 2π to horizon geometry (origin per the 07-05 citation-walk: S217/S201/S380) needs the same note.
  Not done this session — flagged.
- **Born rule.** Gleason's hypotheses are non-contextuality and additivity, not conservation or unitarity; it
  fails in dim 2 and the only worked example is a qubit; Step 2 assumes the conclusion; the conserved quantity is
  state normalization and has no stated relation to C(ρ). The site's badge (Reparametrization) was right and its
  body text was not. There is a sharper open question in this than the page had: the framework's CRT scanning
  picture is a non-contextual *value* model (KS-excluded, 0/512 on Peres–Mermin, B1 row). Gleason needs
  non-contextual *probabilities*. Whether the single-observer picture can motivate the second without the first
  is a real question, and it is the only place the Born-rule track could earn more than a relabel.
- **Compression action.** ξ is never defined as a function of ρ. Forcing the match gives
  ξ^{1/φ} = [(1+x)^{2γ} − 1]/2, which works for *any* exponent p in place of 1/φ. An equivalence that survives
  every exponent gives φ no content. Site badge moved Speculative → Audited-Negative to match the 07-17
  provenance audit of the same exponent.
- **η superconductivity.** Open, not fixed: 607 K for YBCO cannot be reproduced from the page (Δ and η not
  given), and "η ≡ Abrikosov–Gor'kov" is too strong — AG is a digamma-function suppression of T_c, not a 1/η
  factor on the BCS ratio, and an exact reparametrization of correct physics could not be 6.5× wrong. One of the
  Bucket 3 row and the Bucket 2 row is mis-stated. Seeded to the explorer; S616 is the primary source.

**Why this is a research-direction point and not a site chore.** The audit effort since May has followed
visitor traffic and the discriminating-weight rows. That is a reasonable allocation and it produced flagship
pages an adversarial reader could not dent. It also means the *coverage* of the audit is unmeasured: nobody has
counted which archive claims have had any post-May re-derivation at all. The 09-17 finding gave the correction
hazard a denominator over verdict sentences; the analogous denominator over *claims* does not exist. The
quantum/theory tracks (Born, compression action, η, entanglement-as-phase-sync, quantum computing) are the
obvious unswept region, and they are also where the ledger's live Bucket 1 bets sit.

## 3. For dp (gated, unchanged in kind)

- **TEST-26 still has no numeral in its kill criterion** ("robustly requires the crossing", an unnamed SN
  compilation, an unspecified robustness check). This converges with the explorer's 09-18 finding (2 of 26
  criteria well-formed; 11 with no number) from an independent reader on the same day-pair. DR2's own preference
  moves 2.8σ → 4.2σ on SN compilation alone, so an unnamed compilation is a free parameter of the verdict.
  Minimum content before adoption: the named compilation; a Δχ²(ΛCDM − w₀wₐ) threshold; a named non-CPL
  reconstruction; the assumed DR3 error model (the TEST-04a lesson: its "> 0.46 ⇒ > 3σ" implied σ ≈ 0.014
  against a delivered 0.062); and which locality horn (per the 09-15 correction). The physics payoff is
  kill-or-tie-with-Λ, so its value is as the program's first *well-formed* prospective registration — a
  methodology result. It should be framed that way.
- **The frame question two personas reached independently:** with zero active tests that could select the
  framework, the program's live output is methodology — audit latencies (20 days, 26 days, ~9 weeks), criterion
  specificity, auditor calibration, correction hazard. SPINE leads with the ontology and the ledger; neither says
  this. I am not proposing SPINE change — the ontology and Bucket 1 are real and "untested ≠ refuted" holds. I am
  saying the *measured* outputs of the last three months are all of this second kind, and a reader has to infer
  that. One sentence in SPINE's "What it's already good for" would be honest and would not oscillate the frame.
