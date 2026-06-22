# Synchronism — The Prediction Ledger

*The honest scoreboard. Every claim the program treats as a prediction, in one of four
buckets, each row carrying its **named refutation criterion**. This document exists to
**stop the framing from oscillating.** The README and all framing prose must agree with
this table; when they drift (toward overclaim *or* toward self-erasing undersell), this
ledger is the correction.*

**Three disciplines this ledger enforces:**

1. **Confirmed ≠ untested ≠ refuted ≠ nonexistent.** The headline "zero confirmed novel
   predictions" is **true and stays** — but it must not be allowed to *erase* the
   untested-but-falsifiable bucket. "We have made no confirmed novel predictions" and "we
   have several untested falsifiable novel predictions" are **both true**.
2. **Untested is not a win.** A novel-in-structure prediction that nobody has tested is a
   *bet*, not a result. We list bets as bets, with the criterion that would settle them.
3. **Reparametrization ≠ failure, but ≠ discovery either.** All physics reparametrizes
   (Newton↔Kepler, Maxwell↔Coulomb). A reparametrization earns "productive" only when the
   coordinate shift reveals structure / enables a new prediction. None here has — *yet*.
4. **Flat floor ≠ refutation (the saturation reading).** Saturation makes emergence a
   *step* function — flat floor → step → ceiling-that-becomes-the-next-floor (see
   [SPINE.md](SPINE.md) and FUNDAMENTALS §3). A flat result can be a *pre-step floor*, not a
   dead end; a bounded plateau is the mechanism working, not failing. **This is a reading
   guide, not a license to un-refute:** entries in Bucket 2 below were *tested against data
   and failed* (wrong sign, dimensional error, arithmetic that doesn't survive) — those are
   genuine eliminations and stay. The lens only forbids mistaking "nothing yet" (a flat
   floor) for "it failed" (a refutation). They are different claims.

*Companion docs: [SPINE.md](SPINE.md) (the argument), [STATUS.md](STATUS.md) (live
MRH-state), [FUNDAMENTALS.md](FUNDAMENTALS.md) (definitions).*

---

## Bucket 0 — CONFIRMED NOVEL PREDICTIONS

**Count: 0.**

Per the project's own novelty ledger (Session #660A): across ~3,360 sessions, **zero**
predictions that (a) follow from Synchronism postulates, (b) differ from standard physics,
and (c) have been tested and passed. This number is load-bearing and is not softened
anywhere in the repo. When it changes, it changes *here first*, with the test that moved it.

> **Scope: this ledger is the *physics* axis only.** "Zero confirmed" means zero novel
> *physics predictions* — it does **not** mean "useless." On the *applied* axis (Synchronism
> as a design ontology for coherence/emergence/governance) the framework is already
> load-bearing in running code (MRH→Web4, coherence→SAGE, fractal societies→hestia/the hub,
> the live fleet). See **[SPINE.md → "What it's already good for"](SPINE.md)**. Two axes,
> both honest; don't let this one's "zero" erase the other's delivered value.

---

## Bucket 1 — UNTESTED BUT FALSIFIABLE (genuine novel bets)

These are novel **in structure** — each makes a claim that, if tested, could distinguish
Synchronism from the standard frame. **None has been confirmed.** They are the program's real
open bets. Listing them is the antidote to the undersell; flagging them *untested* is the
antidote to the overclaim.

> **Mostly untested for lack of instruments, not lack of merit.** The lab has no accelerator,
> cryostat, interferometer, or neural lab; physics bets here are checked (when at all) against
> data others collected for other purposes. That's why some can be *refuted* (a derivation
> fails) but none can be *confirmed* this way. Two of these would matter well beyond
> Synchronism if a lab ran them — **entanglement-as-phase-sync** (quantum computing /
> semiconductors) and **"time dilation"-as-instrument-effect** — see the **invitation** in
> [SPINE.md](SPINE.md).

| # | Prediction | What's novel about it | Refutation criterion | Honest odds |
|---|-----------|----------------------|----------------------|-------------|
| **B1** | **Observer-relative Bell/CHSH** — a substrate measured only via an observer-pattern's phase-lock to a target (free CHSH setting choice) reproduces quantum correlations without superluminal signaling. | This is the *only* place the single-observer constraint differs structurally from a God's-eye coupled-oscillator model. The one seam where a novel result can live. | Run the harness ([`simulations/kuramoto-lattice-suite/`](simulations/kuramoto-lattice-suite/)). **Refuted if** observer-relative statistics obey CHSH S ≤ 2 OR if violation only appears with signaling / superdeterminism. | **RUN 2026-06-21 → refuted, both arms.** *(a) Local* construction: S = **1.98** ≤ 2, no signaling — a local-realist ontology, not a local-hidden-variable physics. *(b) Nonlocal-grid* construction: S ≡ **2.0 for all coupling g**, no signaling — because smooth single-grid mixing is **gauge-equivalent to relabeling the measurement angles** (uniform shared phase absorbs the offset), so it stays local-realist and buys nothing. **The gap, now precisely named:** a no-signaling violation needs a *non-relabelable, conditional* setting-dependence (a quantum-entanglement-like primitive the ontology lacks and would have to derive). Does NOT change Bucket 0. |
| **B2** | **RAR transition-shape discriminator** (γ pinned at 2): tanh(2·ln(1+x)) compander differs from McGaugh's ν(y) by ~0.083 dex in transition curvature near g_bar/a₀ ≈ 1.1. | The program's single genuinely non-degenerate quantitative test — a place the math is *forced* to differ from MOND. | Marginalized BIC fit on SPARC (153 galaxies), γ fixed at 2, vs McGaugh ν. **Refuted if** ΔBIC favors MOND (S660B already measures structured residual RMS 0.067 dex > σ_int 0.057 dex). **Supported if** ΔBIC favors the compander. | **Likely refuted** by the project's own S660B analysis — but it is the one test, and it should be *run to completion and recorded here*, not left implied. |
| **B3** | **Gnosis consciousness predictions** (~34): consciousness as measurable coherence with a C ≈ 0.50 threshold → specific, falsifiable signatures in anesthesia onset, sleep-stage transitions, meditation, split-brain. | Falsifiable, specific, structurally novel — predicts *where* phase transitions in neural coherence should sit. | Measure coherence (C_conv, C_corr, phase-locking) across anesthesia / sleep / meditation / split-brain and compare transition points to predicted thresholds. **Refuted if** measured neural phase transitions don't cluster near predicted C values. | **Threshold *value* REFUTED; neural program untested.** The C ≈ 0.50 *value* was empirically tested against multi-model coherence data and **rejected at p < 0.0001** by the companion autonomous program [gnosis-research](https://github.com/dp-web4/gnosis-research) (Session 63, 2026-04); the data cluster near **C ≈ 0.64 ≈ φ−1** — a reparametrization candidate (Bucket 3), not a confirmation. So the threshold *location* the 34 neural predictions were keyed to is refuted, leaving them untested **and now mis-anchored**. The identity claim ("phase patterns ARE experience") stays philosophically defensible but empirically ungrounded. *A partial refutation the ecosystem already owned — listed here so the ledger stops under-claiming it.* |
| **B4** | **Compatibility scaling law**: p_crit ∝ 1/⟨C⟩ in heterogeneous multi-agent coherence (Phase-2 extension of the coupling-coherence experiment). | If it holds for arbitrary compatibility structures, it's a genuine law of phase transitions in coupled systems — not just a relabel. | Phase-2 run with explicit compatibility matrix C[i][j] across ≥5 structure types. **Refuted if** p_crit does not scale as 1/⟨C⟩ outside the homogeneous case. | **Moderate.** The homogeneous case is confirmed (r=0.994); the *law* is untested in the heterogeneous regime. |
| **B5** | **f(N) pattern-reconstruction rate**: c emerges as the stable reconstruction rate for minimal-complexity patterns; mass ≡ pattern complexity ⇒ a derivable f(N). | Would convert "c as reconstruction rate" from slogan to derivation, with testable consequences for pattern formation on the grid. | **Pre-condition: it does not yet exist** — the whitepaper points to an Appendix A.19 that is unwritten. **Refuted if** no f(N) with the boundary behavior f(N)→1 as N→0 reproduces stable-pattern simulation data. | **Unknown — not yet derived.** This is an *obligation*, not a prediction, until the derivation exists. Listed so it isn't quietly forgotten. |

> **The disciplined headline** (use this phrasing, not "zero predictions" alone):
> *"Zero **confirmed** novel predictions; a short list of **untested, falsifiable** novel
> bets (above), most of which the project's own analysis expects to lose."* Both clauses
> are true. Dropping either one is the oscillation.

---

## Bucket 2 — REFUTED (tested, failed — productive eliminations)

These were tested against data or proven internally inconsistent. Each eliminated a
possibility. That is the program working as designed (productive failure > safe summary).

| Prediction | How it died | Session |
|-----------|-------------|---------|
| a₀ = cH₀/(2π) as derived MOND scale | Wrong sign; artifact of fitting, not derivation | S438 |
| Transfer rule ⇒ Navier-Stokes | Proven 1-DOF scalar diffusion: irrotational, dissipative, no inertia → cannot produce particle dynamics | S617 / S665 / S666 |
| C(ρ) ⇒ MOND (γ free) | Collapses exactly onto MOND when γ is fit; γ=2-pinned refuted at ΔBIC = +184 on SPARC | S660B / S661 |
| Cosmic interference scale λ ~ 500 Mpc | Dimensional error (units m² not m) masked by a numerical coincidence near ~600 Mpc | S632 |
| 80-orders-of-magnitude unification | tanh(γ·log(...)) saturates within ~1.6 decades for any sharp γ — structural impossibility | S633 |
| DESI growth-rate suppression (Test-04a) | σ₈(z=0) amplitude disfavored at 2.4σ; kill criterion triggered (S668 had over-softened) | S672 |
| GW170817 coherence dispersion (Test-15) | Coupling α calibrated to data, not derived; both branches non-discriminating | S673 |
| Hot-SC T_c formula (Test-18 / OQ005) | YBCO predicted 607 K vs 93 K actual (6.5×); 22/23 standard CM in η notation | S616 |
| Fractal coherence bridge (OQ007, cosmology) | 0/7 scale boundaries predicted; decoherence governs the boundary, framework has no decoherence parameter | S611–614 |
| MRH cluster locality | Reproduces a known Milgrom-2005 locality no-go, not a novel result | S689 |
| a₀ from Jeans analysis | Stated formula off by 614×; "5% agreement" came from a different (unre-executed) computation | S687 |
| Entity criterion (Γ < m) as novel | Standard narrow-resonance / Breit-Wigner condition (textbook QFT since 1951) — last novel-survivor demoted | S660A |
| Chemistry melting points | 53% mean error | (README honest-limitations) |
| Critical exponents | 2× off from observation | (README honest-limitations) |
| **Foundation 3: monotonic saturation "builds the walls" (self-confined particles)** | **Refuted at CA Stage 1 (1D, 2026-06-22):** monotonic-R diffusion AND wave substrates disperse/dissipate at 0% (reproduces S617 / S19 / S665). Self-confined oscillating patterns require a *focusing* nonlinearity, which breaks Foundation 3. Stable patterns ARE reachable on the grid — just not from saturation alone. | `explorations/2026-06-22-phase1-stage1-localized-oscillation-result.md` |

---

## Bucket 3 — REPARAMETRIZATIONS (relabels of known physics)

Real, often useful, but **not novel physics**: they re-express established results in
Synchronism's coordinates. Productive *iff* the coordinate shift later reveals structure or
enables a prediction (it hasn't yet).

| Track | Reparametrizes | Note |
|-------|---------------|------|
| γ ~ 1 "universal boundary" | Debye model θ_D (1912) | ~86% of the "89% validation" are θ_D restatements (S660A) |
| η superconductivity | Abrikosov-Gor'kov pair-breaking efficiency (1960) | Useful reframe (pair-breaking as a materials-design target); not a new prediction (S616) |
| BCS gap-to-Tc ratio | Standard BCS | <1% error = reproduction, not novelty |
| Hückel 4n+2 derivation | Hückel (1931) | Exact reproduction |
| FΣIR M/L decomposition | Standard mass-to-light analysis | Not unique to the framework |
| C(ρ) compander | MOND interpolating function ν(x) | Identical phenomenology when γ free |
| Intent field | \|ψ\|² in all tested regimes | Unfalsifiable as stated; equals the wavefunction wherever tested |

> **The N_corr-divergence test (Kimi 2026-06-21).** γ = 2/√N_corr is applied with a
> *different* N_corr in every domain — galactic N_corr = 1 by fiat, chemistry N_corr from
> correlation length, quantum N_corr → ∞ — with no scale-bridging formalism. That is the
> difference between an underspecification that will *converge* (productive, like early
> "mass" or "entropy") and a vagueness that *fits anything* (degenerate). **The single test
> of whether any of these reparametrizations becomes productive is whether N_corr converges
> on one cross-scale operational definition or keeps diverging.** It is currently diverging.
> Until it converges, "ONE EQUATION" is a research-direction motto, not a unification.

---

## How to read the whole board

- **The contribution is the frame, not the numbers.** The distinctive thing
  ([SPINE.md](SPINE.md)) is the single-observer / CFD ontology. The quantitative tracks
  (Buckets 2–3) are *probes* that mostly taught us **boundaries** — where the frame does
  and doesn't buy anything. That is real knowledge, honestly negative.
- **The live bets are Bucket 1.** Five untested, falsifiable, novel-in-structure
  predictions — the program's actual open frontier. The honest expectation is that most
  lose. Running them is how "zero confirmed" either changes or gets re-confirmed.
- **This board is the anti-oscillation device.** If a doc says "we've unified physics,"
  it contradicts Bucket 0. If a doc says "we've made no novel predictions, full stop," it
  erases Bucket 1. Pin the prose to the board.

*Last updated: 2026-06-21. Update this file **first** when any prediction moves buckets,
then propagate to README/STATUS.*
