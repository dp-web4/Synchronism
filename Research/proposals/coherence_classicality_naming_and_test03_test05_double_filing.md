# Research Proposal: Two Epistemic Consistency Problems Surfaced by Site Audit

**Date**: 2026-05-27  
**Source**: Synchronism-site maintainer session (daily visitor feedback, four-persona pass)  
**Status**: Proposal — open questions, no archive changes yet

---

## Problem 1: "Coherence" Means Its Opposite in Condensed-Matter Physics

### Observation
C(ρ) is named "coherence" and its axis runs C ≈ 0 (quantum/independent) → C ≈ 1 (classical/ordered). But in condensed-matter physics, "coherence" means quantum phase coherence — a property *maximized* in BEC condensates and BCS superconductors. These are the most quantum-coherent systems known to physics.

Under Synchronism's C axis, BEC (N_corr ~ 10²³) and BCS (N_corr ~ 10⁶–10⁹) sit at very low γ (γ ≈ 0.002–0.02) because their massive N_corr places them in the "Collective / correlated" regime. Their C values at any given density also tend to be low relative to uncorrelated systems. So Synchronism assigns "low coherence" to systems that condensed-matter physicists call "maximally coherent."

This means the central variable is named the opposite of its meaning for the audience most likely to find it relevant (physicists working on BEC, superconductivity, quantum computing).

### The Ambiguity
C(ρ) conflates two distinct axes:
- **Classical ordering** (Ising m: 0 = disordered, 1 = ordered) — what C actually measures
- **Quantum phase coherence** (off-diagonal long-range order: maximized in BEC/BCS) — what "coherence" conventionally means in CM

A BEC is simultaneously:
- **High quantum coherence** (off-diagonal long-range order, massive macroscopic wavefunction)
- **Low C** in Synchronism (because N_corr is huge, γ is small, C output depends on density)

These two axes are orthogonal — macroscopic quantum states are simultaneously quantum *and* collective.

### Research Question
Should C be renamed to clarify what it actually measures? Options:
1. **Rename to "classicality"** — directly describes C ≈ 1 = classical/ordered, C ≈ 0 = quantum/independent
2. **Rename to "decoherence fraction"** — maps to the quantum decoherence program
3. **Keep "coherence" but add explicit disambiguation** — state on every physics-facing page that "coherence" in Synchronism is anti-correlated with quantum phase coherence in the CM sense
4. **Introduce a distinction** — "Synchronism coherence" (classical ordering) vs "quantum coherence" (off-diagonal long-range order)

Option 3 is the minimal site fix. Option 1 is the honest rename if the framework commits to the compander/classicality interpretation. Option 2 aligns with the decoherence literature but has its own implications.

### Why This Is Research-Level, Not Just Site-Level
The naming is not cosmetic — it determines which community can engage with the framework:
- Quantum information / QEC community: would read "coherence" as quantum coherence, find BEC/BCS at C≈0 baffling
- Statistical mechanics community: would recognize "classicality" or "order parameter" immediately
- Decoherence / open systems community: would recognize "decoherence fraction"

The site currently uses the word that misleads the most technically relevant audience. The resolution should be a research-level commitment, not just a terminology patch.

---

## Problem 2: TEST-03 and TEST-05 Double-File the Same Regression Result to Opposite Conclusions

### Observation
The Tier 1 test catalog currently has:
- **TEST-03** (ALFALFA-SDSS TFR Scatter): "Kill criterion TRIGGERED — R² = 0.14 environmental term is below the <20% threshold — presumptively FAILED"
- **TEST-05** (RAR Environment Partition): "RAR scatter shows NP2 environment dependence (p = 5×10⁻⁶)" — no kill badge; no failed status

These are the **same underlying hypothesis** (does the RAR/rotation-curve residual carry an environmental term?) analyzed on overlapping datasets. The p=5×10⁻⁶ figure and the R²=0.14 figure are from the same regression — a small-but-significant effect that passes a threshold on significance but fails a threshold on effect size.

### The Epistemic Problem
Filing the same result as both "Failed" (TEST-03) and implicitly "Untested / open" (TEST-05, which has no kill verdict) allows a reader to reach opposite conclusions depending on which card they look at:
- TEST-03 reader: "This failed, R² = 0.14 < 20% kill criterion"
- TEST-05 reader: "This looks like a near-miraculous discovery, p = 5×10⁻⁶"

This is not honest reporting. p=5×10⁻⁶ + R²=0.14 is one signal: statistically significant but explains only 14% of variance, below the pre-registered kill threshold. The pre-registration said: "environment dependence must explain ≥20% of scatter to constitute confirmation." It does not. That is a failure under the pre-registered criterion, regardless of p-value.

The p-value can be small and the effect can still fail the kill criterion — this is the statistical distinction between *significance* (probability of seeing a result at least this extreme under the null) and *effect size* (how much it matters). The kill criterion was wisely set on effect size, not p-value.

### Resolution
TEST-05 should be labeled consistent with TEST-03: "Environment-RAR dependence detected at p=5×10⁻⁶ but R²=0.14 below kill criterion (<20%); hypothesis fails by effect size even if not by significance." This should be stated explicitly, not left for a reader to infer. The two tests should reference each other explicitly and give the same verdict.

Alternatively, consolidate into one entry: "Environment-RAR dependence: detected (p=5×10⁻⁶, TEST-05/SPARC; TEST-03/ALFALFA-SDSS), but R²=0.14 below pre-registered kill criterion. Hypothesis fails by effect size."

### Why This Is Research-Level
The distinction between "p-value significant → open question" and "R² below kill → failed" is a core research decision about what kind of evidence the framework requires for confirmation. If the kill criterion is R²>20%, then R²=0.14 is a failure regardless of p-value. The current double-filing implicitly lowers the bar (p-value is sufficient) while the pre-registration raised it (effect size required). This inconsistency should be resolved at the research level and reflected consistently on the site.

---

## Summary

| Problem | Type | Site Fix | Research Fix |
|---------|------|----------|--------------|
| "Coherence" = classicality inversion | Naming/framing | Add disambiguation note to physics-facing pages | Commit to rename or explicit anti-correlation statement in framework documents |
| TEST-03/05 double-filing | Epistemic consistency | Add kill-verdict to TEST-05; cross-reference | Clarify in test catalog that R²<20% = failed by effect size |

Both issues originate in the same pattern: a term or result was introduced to describe one thing, and then applied to a context where it means something different, without flagging the shift.

---

## Explorer adjudication (2026-05-27) — Problem 1: the rename options above are wrong

**This proposal's Problem 1 conflicts with `gamma_definitional_collision_regime_label_inversion.md` (2026-05-04), and that earlier proposal is closer to right.** Reconciliation (full analysis: site finding `coherence-naming-three-axis-adjudication.md`):

C = tanh(γ·ln(ρ/ρ_crit+1)) is a function of **density ρ and collectivity N_corr only**. It is high when (dense AND single-particle), low when (sparse OR collective). It contains no ℏ, temperature, action, or decoherence rate — i.e. **zero quantum-vs-classical content**. "Classicality" and "decoherence fraction" (Options 1 and 2 above) are quantum/classical-axis names; attaching them to C does not add the missing axis, it converts a confusing label into an asserted falsehood (it would *assert* that a lone electron, γ=2, high C, is "maximally classical," and a BEC, low C, is "maximally quantum"). **Recommend against Options 1 and 2.**

The 2026-05-04 proposal's Case 1 is correct: the γ-axis is a **collectivity** axis (single-particle ↔ collective), explicitly not quantum/classical. The tools (coherence-explorer, phase-boundary-visualizer) already implement that relabel. What remains unfixed:
1. The relabel never propagated to front-of-site copy (landing page, why-synchronism, coherence-function still say C = "how quantum or classical").
2. Even the fixed tools left **C itself** labeled "0=quantum, 1=classical" — that label must also go.
3. The **defining metaphor is inverted**: landing page / glossary / coherence-explorer say "lockstep, e.g. electrons in a superconductor = high coherence," but a superconductor (large N_corr → small γ) is pinned at C≈0 by the equation. The marching band the framework uses to *define* high coherence is a low-C system by the framework's own equation. (terms.ts even contradicts itself: the γ entry says "γ≪1 = quantum" while the N_corr entry calls the crystal at γ≈10⁻¹² "classical.")

**The actual decision is a scope decision, not a vocabulary one:**
- **Option I (scope restriction):** C(ρ) is a density→saturation map valid at the single-particle reference (N_corr=1, γ=2) — the galaxy regime where it was calibrated. Drop the BEC-to-cosmos "quantum→classical" framing. Consistent with the audited state (compander, 0 discriminators) and with the gamma-dual-role conclusion that C degenerates for N_corr≫1.
- **Option II (keep multi-scale, state the truth):** "C is a density-driven saturation index normalized by collectivity; high for dense weakly-correlated matter, low for sparse *or* strongly-correlated matter; it does NOT measure quantum-vs-classical character — macroscopic quantum systems (BEC/BCS) are strongly correlated and therefore low-C though maximally quantum-coherent."

Deepest point for a framework named *Synchronism*: C doesn't measure synchronization either — a marching band and a BEC are both maximally synchronized (collective) and both low-C. A genuine "lockstep" order parameter (Kuramoto r, ODLRO) is a different object from C(ρ). The naming inversion is the lossy multi-axis→one-scalar projection (compander class, γ dual-role, three-C) made visible at the front door.
