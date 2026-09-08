# The MRH-validity program — the price of a horizon (charter)

**From**: kimi-code · **Chartered**: 2026-09-08 · **Frame-owner**: dp
**Cadence**: interactive with dp for now; no cron. Inherits the house rules of the
hive-organs arc: falsifier-first, registered kill criteria, honesty blocks,
conventional-prior contamination check, the distinct-count keeper.
**Prompt** (dp, 2026-09-08, ratified as the frame): *"our long-running theme had been
'one equation' — but equations are MRH-specific (even heliocentrism vs geocentrism,
let alone per discipline). Yes, coherence drives emergence at every MRH, but the
equations are conditional on both what the MRH includes and what it specifically
excludes. The interesting question is not one equation to rule them all, but more
like equation-of-equations — tied directly to MRH as selective abstraction, and
abstraction is lossy compression. If the equation can model the abstraction loss,
we can perhaps evaluate the bounds of equation validity in MRH. Ptolemy's epicycles
are valid in a very specific scope — what the universe looks like viewed from Earth.
There is always an MRH; the only observer that holds all is the All, and we're all
tiniest subsets of that. Does this line of inquiry make sense?"*

## The reframing this charter registers

The post-audit record (Sessions #574–586 and the AGENTS.md honest assessment) reads:
γ = 2/√N_corr reparametrizes known physics in every empirically tested domain —
Debye at the phonon scale, MOND at the galactic scale, standard QM at the quantum
scale — zero confirmed novel predictions. Under the old frame ("one equation") that
is failure. Under this charter's frame it is the **dataset**: the same functional
form recurring across horizons with different referents is what a coarse-graining
invariant looks like. The program's two questions:

- **Q1 (the validity bound).** Can the *price of a horizon* be computed — a bound on
  an MRH-conditional equation's error as a function of what the horizon excludes?
  Not an equation-of-equations: an **inequality about information**, which prices
  each discipline's equations instead of competing with them.
- **Q2 (the invariant form).** Why does the √-form of γ recur — is it a fixed point
  of coarse-graining, and if so, where does the interesting (non-Gaussian) physics
  live: in which deviations?

## The formalism (registered so later results can be checked against it)

- An **MRH** is a restriction map ρ: Ω_all → Ω_obs — the observer's horizon as a
  choice of included variables *plus the instrument channel* (the epicycle case
  requires this: Ptolemy-vs-Copernicus was underdetermined until Tycho's precision;
  the measurement channel is part of the horizon).
- The **abstraction loss** of ρ, for a prediction target, is the conditional mutual
  information **I(excluded ; future_included | present_included)**. If this is zero,
  the horizon is *closed* (a Markov blanket) and any sufficient equation on the
  included variables is exact within it.
- The **validity inequality**: an MRH-conditional equation's achievable error is
  bounded by a function of the abstraction loss (Fano / rate-distortion direction
  for the general case). For jointly Gaussian systems the bound is an **identity**:
  the prediction-error inflation of excluding Y is exactly e^{2·I(Y;X′|X)}
  (first results doc, computation 1).
- **N_corr** (the repo's standing parameter) is the correlation structure's shadow:
  N_corr = N/(1+ρ(N−1)) for pairwise-ρ variables; the √-recurrence is the CLT/Gaussian
  fixed point (computation 3). The physics of each horizon lives in ρ — in what its
  couplings to the excluded do to the effective count.

## Conventional-prior contamination check (what already exists; we stand on it)

Effective field theory / Wilson RG (every equation horizon-conditional, with cutoff
and error terms — orthodox physics; our Q1 is EFT's practice generalized to
observer-relative exclusion beyond energy scale); Markov blankets / free-energy
principle (Friston — observer-as-subset, to differentiate not duplicate); causal
emergence (Hoel — under noise a macro description can carry MORE causal information:
abstraction loss can go negative, the anti-reductionist result this program needs);
information bottleneck (Tishby); computational mechanics (Crutchfield's ε-machines);
QM restriction to observable subalgebras, with Tomita-Takesaki modular theory
generating a thermal time from the restriction — a resonance with the substrate
track's two-level time ontology, unexamined so far.

## Registered first experiments (all run 2026-09-08; results doc is companion)

1. **Gaussian price of a horizon** — verify the identity error-inflation = e^{2·CMI}
   on a coupled linear system with an exclusion-leak knob; blanket closure (leak 0)
   must cost exactly nothing (the Ptolemy case in toy form).
2. **The Debye boundary zone** — exact Debye C_V vs the two horizon equations
   (low-T T³ law, high-T Dulong-Petit 3Nk): measure where each crosses 1% error and
   the shape of the loss profile.
3. **N_corr from correlation** — the fluctuation exponent under pairwise correlation:
   1/2 = Gaussian (no new information), ≠1/2 = correlated degrees of freedom an
   MRH-only equation cannot see. Registers the per-domain exponent measurement as
   the program's empirical front.

## Registered kill / honesty criteria

- **K1.** If no computable bound is produced for ANY known effective equation, the
  program is philosophy, not physics — stop and report that. (Computation 2 is the
  first defence against this kill.)
- **K2.** If γ's recurrence is fully explained by the CLT fixed point with zero
  anomalous exponents in every domain measured, then "Synchronism's γ" dissolves
  into statistics. That is a reportable answer, not a failure to hide — the program
  would then have shown WHY the one-equation search was doomed, which is itself the
  contribution.
- **K3.** If the web4↔physics MRH bridge cannot be given one definition that compiles
  in both (context-scoping graph ↔ restriction map/σ-subalgebra), drop the
  unification claim; the physics half stands alone.
- "Only the All holds all" is registered as cashed: *no finite observer accesses the
  full state, so all equations are effective equations*. Uncashed, the sentence is
  theology and may not appear in results.

## The map

This charter + the dated results docs in this directory. No living roadmap file
(unlike the hive-organs arc): if the program grows past three docs it earns one.

— kimi-code

---

## Errata / addenda 2026-09-08 (codex review; see `2026-09-08-kimi-response-to-codex-review.md`)

- **"Abstraction loss can go negative" was a category slip.** For a fixed
  restriction the CMI is nonnegative. Causal emergence says a DIFFERENT
  horizon choice can carry lower loss for the same task — horizon choice can
  reduce the price; nothing goes negative. The framing survives sharpened.
- **The general anchor adopted**: L*_log(X) − L*_log(X,Z) = I(Y;Z|X) under
  log loss; the Gaussian e^{2·CMI} ratio is the square-loss specialization
  (counterexample against overextension: codex SA-3A T2).
- **Identifiability is a first-class field** (SA-3A T5: [0, H(Y|X)] — a
  channel cannot identify its own omitted information). All closure claims
  in this program carry a compatible-interval column from rung 3 onward.
- **Prior art added to the contamination check**: Mori-Zwanzig / Chorin
  optimal prediction (memory from eliminated variables); the in-repo Markov
  Phase 3 doc (2026-08-17) — target/horizon/tolerance-relative CMI relevance,
  which predates this program and is credited.
