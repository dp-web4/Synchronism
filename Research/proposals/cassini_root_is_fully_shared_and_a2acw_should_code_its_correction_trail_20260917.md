# Proposal: the Cassini root is fully shared with MOND, and A2ACW should code its own correction trail

**From:** site maintainer track, 2026-09-17 (WAKE). **Count stays 6; Bucket 0 = 0.** Two ledger-governance items gate on dp.

## Why this is a research item, not only a site fix
Two of today's visitor personas (graduate physics, leading-edge researcher) independently reached the conclusion the
09-16 proposal reached: **the credibility risk has moved to the refutation side.** The arguments about *why* things
fail now draw more referee objections than the claims do. The 09-16 proposal argued this in general terms. Today
there are two concrete cases, one resolved by execution and one about what the program should measure next.

## 1. TEST-25 (Cassini) is inherited in the literal sense. The asymmetry argument is withdrawn.

**Context.** The site gave three incompatible accounts of the Solar-System failure:
- /galaxy-plotter: McGaugh's RAR ν = 1/(1 − e^−√y) is "Cassini-safe".
- /galaxy-rotation: the empty intersection is a tail-shape mismatch (power law vs e^−707 at Saturn).
- /honest-assessment: "MOND picks a different μ and survives. Synchronism cannot."

Cassini bounds the external-field-induced quadrupole Q₂. Q₂ is set near √(GM☉/a₀) ≈ 7000 AU, where g ~ a₀, so it
probes the interpolating function (IF) in its transition, not its tail at Saturn. McGaugh's ν is the δ = 1 member of
the δ-family that Desmond, Hees & Famaey 2024 put at 8.7σ.

**Executed.** Pre-registered at site commit `43a66a3` before computing. The instrument is TEST-25's own
`simulations/sparc_cassini_q2.py`. Script and output:
`synchronism-site/maintainer/scripts/cassini_q2_mond_interpolating_functions{.py,_PREREG.md,_output.txt}`.

Controls:
- C1: the compander at γ = 0.489, a₀ = 5.33265e−11 reproduces z = +17.95.
- C2: δ = 1 matches `nu_rar` exactly.

Results on the current Cassini interval, grid a₀ ∈ {1.128, 1.20}×10⁻¹⁰ × g_ext ∈ {2.00, 2.32, 2.48}×10⁻¹⁰:

| IF | z range | inside 95% |
|---|---|---|
| McGaugh RAR ν (δ = 1) | +15.9 … +20.9 | 0/6 |
| Milgrom simple (n = 1) | +15.3 … +20.1 | 0/6 |
| standard (n = 2) | +9.5 … +10.8 | 0/6 |
| δ = 2.5 | +5.5 … +7.7 | 0/6 |
| δ = 3 | +3.2 … +5.8 | 0/6 |
| δ = 4 | +0.9 … +3.3 | 3/6 |
| compander γ = 0.489 (TEST-25 point) | +17.95 | no |
| compander γ = 1.5 / 2 / 3 (**post-hoc**, g_ext 2.32) | +1.1…+7.9 / −0.4…+3.2 / −0.9…−0.1 | a₀-dependent / a₀-dependent / yes |

Predictions: P1 (δ = 1 fails at |z| ≳ 5) **held**. P2 (simple fails harder than δ = 1) **refuted**: it fails
slightly less. P3 (a passing δ near 2.5) **held only at δ = 4**, on half the grid, on this unmarginalized instrument.

**Reading.** Both families have the same structure. Their RAR-preferred members fail Cassini at comparable z, and
their sharp-transition members pass: MOND's at δ ≈ 4, the compander's at γ ≳ 1.5–2. For the compander, SPARC
excludes the passing members (γ = 2 at ΔBIC +184; the retained interval ends at 0.600). For MOND, Desmond+ report
that the RAR prefers δ ≈ 1 and that the tension persists across the families they tested. The SPARC cost of δ = 4 was
not computed here. So "MOND picks a different μ and survives" is unsupported by the program's own instrument and by
its own citation. **TEST-25 stays booked "inherited from MOND", and the inheritance is literal.** The direct-tail
argument (e^−707) is a real but separate effect and is not what TEST-25 computed.

**Asks (dp-gated):**
- (a) Confirm TEST-25's booking stays "inherited" with the asymmetry clause withdrawn from the archive's own framing.
- (b) If the preprint cites TEST-25, cite it as a QUMOND-family result reproduced on the program's instrument, not as
  a framework-specific kill.

**Next (unexecuted, for explorer):** run the frozen SPARC likelihood at δ ∈ {3, 4}. If δ = 4 costs ΔBIC ≫ 10, the
Cassini root is a QUMOND-class closure: no single universal IF fits both. That is a transferable null with the
program's instrument benchmarked to Desmond+ at 0.76%.

## 2. A2ACW: the better experiment is the correction trail already in git

**Context.** /a2acw was still stating the superseded null ("6 externally-audited", "human audit", "sensitivity 6/6",
"citable null"). That is fixed on the site today. The researcher persona then questioned the frame of the proposed
post-cutoff control:
- Both benchmark arms are **famous** items, so they test recall rather than novelty judgment.
- There is **no human-referee arm**, so "LLM audit maps novelty onto prior art" has no reference rate.
- The program already holds a better dataset: **hundreds of dated corrections**, each with a direction, the track
  that caught it, what caught it, and a latency.

**I agree, and I think it changes what A2ACW should measure.** The archive's A2ACW question (H1: nothing novel vs H2:
an LLM auditor rewarded for prior art over-maps) is one slice of a broader, answerable question: *what is the error
profile of LLM research agents that have an executable oracle?* The correction trail answers it without a post-cutoff
control, because the ground truth for most corrections is an execution (a re-run script, a control, a paper read in
full), not another LLM's judgment. Recent instances show both directions:
- 09-16: compander 2.10× withdrawn by controls.
- 09-15: the "no instability" claim missed the Jeans term.
- 09-17: the Cassini asymmetry was withdrawn by execution, and its opposite, "Cassini-safe", was also refuted.

**Proposed design (for explorer; pre-register before coding):**
- **Unit.** A correction event: a dated "corrected / withdrawn / retracted / restated" marker in PREDICTIONS.md, site
  `src/`, or the maintainer/explorer logs, deduplicated by object.
- **Codes:**
  - direction: over-claim, over-refutation, or neutral/provenance
  - caught-by: execution/control, primary-source read, re-derivation, or persona objection
  - catching track
  - latency in days from first appearance (git blame) to correction
  - whether the correction itself was later corrected
- **Inter-rater.** Two independent codings (different model vendors if available). **dp codes a random 10% as the
  human arm.** That gives the human-referee reference rate the benchmark lacks.
- **Registered predictions to test:**
  - over-refutation's share rises after 2026-08-01
  - execution-caught corrections have lower re-correction rates than read-caught ones
  - latency is longer for over-refutations
- **Why it matters beyond this program.** It is a measured error taxonomy for autonomous AI science with an oracle,
  and "0 of 9" cannot be interpreted without one.

The existing topic `correction-palimpsest-rate.md` (text share of correction notes vs contradiction reports) is
compatible and should share the event extraction.

## 3. Open gap surfaced, not executed: the matter budget
The dark-energy sector fits DESI + Planck with Ω_m ≈ 0.315 of *clustering* matter and sets C₀ = Ω_m. The galaxy sector
is meant to replace dark matter, and its most permissive boost normalization, Ω_m/Ω_b ≈ 6.4, is the ΛCDM baryon
fraction. Does the framework's cosmology contain non-baryonic clustering matter?
- If yes, the galaxy sector double-counts.
- If no, C₀ = Ω_m has no referent and the DE fit is a fit of a different model.

No archive document reconciles the two. This is seeded to the explorer.

## Also back-annotated today
- `simulations/efe_locality_vs_phi_dependence.py` docstring: "EFE = 5.6e−13" certifies superposition. The refracted
  host field was subtracted (explorer 09-16).
- PREDICTIONS.md: a 📌 block covering the three items above, plus the GC window as an L2 object.
