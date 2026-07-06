# Proposal: Dim-4 LIV Exclusion as a General No-Go for Absolute-Time Substrates

**Date:** 2026-06-30  
**Source:** Pass 4 researcher (leading-edge LIV/quantum-foundations specialist) site review  
**Upstream computation:** Phase-13 (2026-06-23), Phase-16 (2026-06-24), explorer triage 2026-06-26

---

## The Finding (as-is, in Synchronism terms)

The framework's dim-4 LIV computation established:

- Tree-level c_μν = 0 (single-substrate universality — all species share one substrate, so no species-dependent k² coefficient at tree level)
- Radiative c_μν ~ α/π ~ 10⁻² (Collins–Perez–Sudarsky–Gambini–Pullin, PRL 93, 191301, 2004 — one-loop, UV-dominated, Planck-cutoff-independent)
- Existing bounds: |c_μν| ≲ 10⁻¹⁸ (cavity Michelson–Morley) to 10⁻²⁹–10⁻³⁰ (nucleon comagnetometer)
- Fine-tuning gap: 16–28 orders of magnitude
- Root cause: absolute time / global preferred frame removes boost invariance, the only custodial symmetry that could forbid dim-4 LIV radiative generation

## The General Form (beyond Synchronism)

Pass 4's question: **does this exclusion generalize to any discrete absolute-time substrate, or is it specific to Synchronism's parameter choices?**

The Collins–Perez–Sudarsky result is a theorem about *any* Lorentz-violating UV completion at a discrete cutoff, not about Synchronism specifically. Its content:

> If a theory has a preferred frame (breaking boost invariance) and couples to the Standard Model at the Planck scale, dimension-4 LIV operators are generated radiatively at O(α/π) unless a symmetry forbids them. The only symmetry that can forbid them is boost invariance — which a preferred-frame theory has renounced by definition.

This creates a general no-go:

**Theorem (CPSU 2004, restated):** Any quantum field theory with (a) a discrete Planck-scale substrate, (b) a preferred frame (absolute time), and (c) SM-perturbative gauge couplings, generates c_μν ~ α/π at one loop, refuted by existing bounds by 16–28 OOM.

## Why This Is Broader Than Synchronism

The three conditions (a), (b), (c) apply to a wide class of proposals in quantum gravity / digital physics:

- Lattice QFT with a preferred time slicing
- Causal dynamical triangulations (CDT) variants with preferred foliation
- Horava–Lifshitz gravity in its non-projectable form (broken Lorentz at the UV)
- Any "ether" or "granular spacetime" model that assigns physical significance to a global clock

The Synchronism-specific step is identifying that the framework's ontological commitment to absolute substrate ticks (Level-0 time) directly instantiates condition (b). But the refutation doesn't depend on the rest of Synchronism being true — it depends only on (a)+(b)+(c).

## What Would Make This a Publishable Negative Result

The path to a citable note (PRL Comment, LRR review section, or arXiv:hep-th):

1. **State the general theorem** (CPSU 2004) clearly, with the three conditions
2. **Apply it as a triage criterion**: any proposal that commits to (a)+(b)+(c) is refuted at natural radiative value — no fitting needed
3. **Scope the escape routes honestly**:
   - Add boost invariance → gives up (b) → preferred frame is gone
   - Non-perturbative strong-coupling UV completion → requires exhibiting the mechanism; constrained by Bednik et al. 2013 (log running closes only ~1.6 OOM)
   - SUSY → custodial, but independent commitment
4. **Benchmark against Synchronism as the worked example**: tree=0, radiative=10⁻², bounds=10⁻¹⁸–10⁻³⁰ → 16–28 OOM

## Relation to Existing Literature

This is NOT a new theorem. CPSU 2004 already establishes it. What this program contributes:

- The explicit instantiation in a worked example (Synchronism's single-substrate universality → tree=0, then radiative overwhelms)
- The identification that "single-substrate universality" creates a *distinct* tree-level suppression that doesn't help at one loop
- A clear triage criterion applicable to contemporary digital-physics proposals

Honest novelty: we are writing a **quantified instance + application** of an existing result, not a new no-go. The publishable framing is as a "worked application of CPSU 2004 to discrete absolute-time substrates, with a triage criterion for the contemporary emergent-gravity literature."

## Action

- **Decision gate for dp**: is this worth writing up as a short arXiv note (hep-th, ~4 pages)? The CPSU citation is solid, the computation is clean, and the scope is precise.
- **If yes**: the writeup should cite CPSU 2004 as the primary result and Synchronism as the worked example, not the primary claim.
- **If no**: at minimum, update `/for-researchers` to frame the dim-4 result as a transferable no-go rather than a framework defeat — which makes it more credible to arriving experts, not less.

**Status: this proposal is the research-direction signal from the 2026-06-30 maintainer session WAKE. Decision belongs to dp.**

---

## Amendment 2026-06-30 (explorer novelty/feasibility audit)

Explorer ran the writeup-feasibility audit against the external literature (`synchronism-site/explorer/findings/dim4-liv-nogo-feasibility-escapes-missed.md`). Two corrections to this proposal:

### 1. The escape-route triage is incomplete and partly contradicts itself

- **Internal contradiction:** §"The Finding" line "boost invariance [is] the only custodial symmetry that could forbid dim-4 LIV" is contradicted by this proposal's own escape bullet "SUSY → custodial." SUSY **is** a custodial symmetry for exactly this problem (Groot Nibbelink & Pospelov, *PRL* 94, 081601, 2005: SUSY forbids dim-3/4 LV operators, forcing them to dim ≥ 5, UV-suppressed, with no high-energy dispersion modification). Boost invariance is therefore **not** "the only" custodial symmetry — delete that phrasing.

- **Missing the decisive escape — anisotropic scale-hierarchy (Pospelov & Shang, *PRD* 85, 105001, 2012).** This is the most damaging omission because it is demonstrated for **Hořava–Lifshitz gravity — the canonical preferred-foliation (absolute-time-class) theory**, i.e. condition (b) of the restated theorem. Mechanism: the LV sector acquires Lifshitz scaling above Λ_HL and couples to the SM only through M_pl-suppressed operators; the percolated SM-sector LV is then protected by Λ_HL²/M_pl², and a wide separation Λ_HL ≪ M_pl removes the tuning **without restoring boost invariance and without SUSY, perturbatively.** This directly refutes the site's stronger claim (`/for-researchers`) that any survivor must be "boost-symmetry-independent AND non-perturbative." Scale-separation is boost-breaking-compatible and perturbative.

- **Causal sets (Bombelli–Lee–Meyer–Sorkin 1987; Bombelli–Henson–Sorkin, "discreteness without symmetry breaking"):** random Poisson sprinkling gives a discrete substrate with **no preferred frame** — Lorentz-invariant on average, no modified dispersion. This shows the *discrete* class as a whole is not excluded; only the *regular-lattice + preferred-frame* subclass (which Synchronism's absolute time keeps) inherits the problem.

### 2. Reclassify: naturalness problem, not refutation

A naturalness problem the field **routinely solves** (SUSY, scale-hierarchy) is not a "decisive falsification." The honest status of the dim-4 face is:

> The *minimal* framework (no added custodial structure) predicts c_μν ~ α/π, refuted by 16–28 OOM — the standard LV-naturalness problem (CPSU 2004). It is **not doubly-obstructed**: SUSY and anisotropic scale-hierarchy are standard, perturbative, boost-breaking-compatible custodial escapes. Synchronism adopts none (single-substrate universality is the obstruction to a Λ_HL ≪ M_pl hierarchy; SUSY/random-sprinkling are unmotivated additions or require giving up absolute time), so the problem is **unresolved — open custodial-mechanism gap, not closed refutation.**

Keep "minimal framework refuted at natural value" (true, conditional). Drop "doubly-obstructed / most decisive falsification / only-non-perturbative-escape." The genuinely decisive negatives remain the **locality no-go** and the **DESI sign failure** — those have no comparable standard escape.

### 3. Feasibility verdict: NOT a standalone paper

CPSU 2004 is one of the most-absorbed results in LV phenomenology (reviewed in Mattingly 2005; Liberati 2013; Collins et al. hep-th/0603002; mature Polchinski-era debate, arXiv:1106.1417). The community does not ignore it — it *solves* it (Hořava: scale-hierarchy; causal sets: random sprinkling; LQG: known tension). A referee would see the dim-4 c_μν argument as CPSU applied to one more substrate, and would flag the "only custodial symmetry" claim as wrong. **Standalone arXiv note: not feasible.** The only viable form is a short *critical comment* on the digital-physics / absolute-time-substrate class (Synchronism, and — if they likewise ignore CPSU — 't Hooft cellular-automaton QM, Wolfram model), noting each inherits the naturalness problem and specifying which custodial mechanism it must adopt — **and only after the escape-route correction above.** This is the **third** "citable result" to reduce to *known-theorem-instance + contemporary-triage* under novelty audit (cf. locality = Milgrom instance; A2ACW = assembled prior art); it ranks *with* them, not above them.
