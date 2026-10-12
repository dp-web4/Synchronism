# The counting rule admits the uniform growth reading, and the post-cutoff arm is the one live experiment

*Site maintainer, 2026-10-12. Source: synchronism-site visitor log 2026-10-12, Pass 4 (researcher persona) and
Pass 3 (graduate persona). No bucket moves. Count 6. Bucket 0 = 0. Two items gate on dp.*

## 1. The site's counting rule, read as written, admits a seventh refutation

The rule is printed on the site (glossary `refutation`, since 2026-09-14): *a claim counts as an executed
refutation when (1) a test with a stated failure condition was actually run, on data or as a simulation, and (2) the
claim failed it.*

The uniform G/C growth reading (U) meets both halves:
- It applies the framework's own substitution G → G/C uniformly to perturbations. That is the same move that builds
  the Friedmann equation in Session 100. It gives fσ₈(0.51) = 0.575 (+21 %) and μ(a) = 1/Ω_m(a) for every γ
  (maintainer 2026-10-07, `synchronism-site/maintainer/scripts/growth_under_uniform_G_over_C.py`).
- It was pre-registered (site explorer, commit `be3ecba`, rule: ≥ 2 independent probes at ≥ 3σ) and then executed.
  Without slip the ISW–galaxy cross-correlation flips sign (A ≈ −8.5 against an observed ≈ +1.0 ± 0.2), and
  Σ₁ ≈ 1.6 against DESI DR1's 1.021 ± 0.029. With lensing left unmodified, each of two independent legs excludes it at
  ≈ 3.5σ. (The leg count was corrected after the run: the registered "growth" row carried DES Y3 lensing. See
  `synchronism-site/explorer/findings/uniform-g-over-c-growth-is-excluded-...md`.)

It is not in the "6". The only reason is one nobody wrote down: no archive document asserted reading (U). The
site constructed it in October as the plainest operationalization. That reason runs in the framework's favour. A
framework that never commits to how C enters perturbations can never lose a growth test. The persona put it plainly:
the γ = 2 *pin* counts and the excluded *plainest reading* doesn't, and that asymmetry favours the framework.

**The two consistent options, and what each does to the count.**
- **(a) Count it.** Refutation #7, from a new root: how C couples to perturbations. That root is independent of the
  boost cap, the γ = 2 shape, the density keying, Cassini and CHSH. The count becomes 7 refutations from 6 roots.
- **(b) Add an "asserted by the archive" clause and apply it to every row.** (U) drops out, and the clause then has
  to be checked against TEST-25. Its scale-universal QUMOND realization (the same compander used as the Solar-System
  interpolating function) was set up in-house on 2026-07-23, and the ledger already calls it "realization-refuted,
  umbrella untested". Whether "ONE EQUATION" counts as the archive asserting scale-universality is a judgment I can't
  make for dp. If it doesn't, the count under (b) is 5.

So "6" sits between two consistent rules (5 under a strict (b), 7 under (a)), and as stated it matches neither.
**Recommendation:** (a). The rule as written already counts executed, pre-registered exclusions of the framework's own
equations. "Nobody in the archive wrote this reading down" is the operationalization gap the 10-07 proposal flagged,
and that gap shouldn't double as an exemption. Meanwhile the site reads "6, plus one excluded reading pending a
counting decision" (Honest Assessment §1b and the glossary, 2026-10-12).

## 2. The post-cutoff A2ACW arm: fifth independent raise, 32 days at the gate

The researcher persona's research-direction note, verbatim in substance: every pending physics test is externally
gated (DESI DR3 ~2027–28, Gaia DR4 late 2026, DR2 full-shape ~2027), and the site itself says none of them can select
the framework. The post-cutoff arm is the only experiment that is:
- internally runnable,
- of unknown outcome, and
- able to produce a positive, citable result: a calibrated novelty detector, or a calibrated failure.

Raises so far: visitor 2026-06-27, 2026-07-02, 2026-09-10 (routed to dp in
`test02_amplitude_is_knee_conditional_and_a2acw_positive_control_20260910.md`), 2026-10-08 (maintainer, third raise),
and today. The explorer's 2026-06-27 finding (`a2acw-specificity-measures-framework-reuse-not-novelty.md`) adds a
caveat the run design must carry. A2ACW sorts "new framework vs new object in an old framework", which is orthogonal
to discovery value, so a clean post-cutoff specificity could still come out ≈ 0 for a reason other than blindness.
The run has to score both axes or it repeats the confound.

Today's site cleanup also changed what the oracle thesis rests on. The temporal-asymmetry "0/6" was shown as a run on
Honest Assessment and /research-philosophy, and was a desk counterfactual. So the support is now:
- 5/6 refutations came from code run against external data;
- 0/4 prior-art identities arrived through the protocol;
- 4/4 arrived through retrieval lanes.
All three are observational and uncontrolled. The "cannot be removed by more capable agents" half of the thesis has
never been tested. The post-cutoff arm is that test.

**Ask (dp):** authorise the run, or rule it out permanently. Nothing would be sent to an external service: the arm
uses post-cutoff papers and the existing protocol. A third option, leaving it queued, is in effect a ruling too:
after five raises the methodology claim is being cited while its decisive control stays unrun.

## 3. Recorded, no direction claim

The researcher persona said all four of its version-drift findings "run in the framework's favour". Today's
corrections, sorted by the direction the stale text leaned:

| Correction | Stale text leaned |
|---|---|
| TEST-04a card hid that 0.418's mechanism was withdrawn | pro-framework (a test looked live) |
| Uniform G/C reading off the ledger | pro-framework |
| A2ACW 0/6 desk counterfactual shown as a run | pro-methodology-claim |
| "The methodology finding is" | pro-methodology-claim |
| cH₀/2π "inside ±5%" (12% is not) | pro-framework |
| Ciocan σ's divided 95% CIs as 1σ | anti-framework (understated the constant-a₀ exclusion) |
| Plotter formula box vs legend (opposite wirings) | neither |

That's 5 of 7 leaning toward the program, from a sample picked by one persona who was looking for exactly that, with
no denominator for the drift that wasn't reported. Not a finding. The rate is what's worth measuring: the
correction-trail coding (topic done 2026-09-17) could be extended to site corrections.

## Also back-annotated today (no proposal needed)
- TEST-25 is conditional on a modified-Poisson completion. The DHF 2024 bound is QUMOND/AQUAL; modified-inertia MOND
  evades it. The framework's only written field equation is modified-Poisson, so the inheritance holds for it. The
  published 8.7σ covers the family near γ ≈ ½; the fast-return tail (γ ≳ 1.5–2) passes Cassini post-hoc and is
  SPARC-disfavoured at ~2σ galaxy-level.
- The a₀(z) row is now on the level-free slope: estimator-dominated and degenerate with ΛCDM at n ≈ 1. Constant a₀ is
  excluded at 5σ by Ciocan's bins.
- A graduate persona asked about the near-Earth density boundary for ∇·[C(ρ)∇Φ] = 4πGρ (atmosphere vs LAGEOS, ~10⁻⁹).
  It's an instance of the 2026-09-15 smoothing-length obligation, not a new test. With ℓ ≳ 30 AU (planets) the
  atmosphere is smoothed away; pointwise, it fails the way the solar wind does. The page should say which ρ_crit and
  which ℓ apply outside galaxies, and the archive doesn't. Explorer topic seeded.
