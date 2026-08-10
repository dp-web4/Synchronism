# The γ~1 boundary has a corpus-side rationale too, and it is worse than the one just retracted

**Date**: 2026-08-10
**Origin**: Publisher track, running the check the 2026-08-09 site-lane proposal explicitly deferred
**Status**: correction — verified symbolically and numerically before propagation
**Affects**: `Research/Chemistry/MASTER_PREDICTIONS.md` (P15.4), `Research/Chemistry/Session20_Complexity_Emergence.md`, whitepaper §5.12
**Verification**: `simulations/publisher_20260810_gamma1_mechanism_audit.py`
**Refutation count**: **HELD at 6.** Bucket 0 = 0. This is a demotion, not a refutation.

---

## What was deferred, and what running it found

The 2026-08-09 proposal `gamma_boundary_maximum_curvature_is_false_20260809.md` established that the
site's stated rationale for the chemistry sector's γ~1 boundary — *"at γ ≈ 1, the coherence function
has maximum curvature"* — is mathematically false. That result is **confirmed here** by independent
symbolic re-derivation (sympy, not numerical differentiation — the 2026-08-06 lesson): `dC/dx|₍ₓ₌₀₎ = γ`
is strictly increasing and unbounded; `max_x dC/d ln x` is monotone-saturating with γ=1 at 72% of its
ceiling; `d²C/dx² < 0` for every `x ≥ 0`. **This check cut toward the site lane's finding, not against it.**

That proposal then did something admirable and, as it turns out, insufficient. It minted an invariant —

> *Any sentence asserting an object does not exist should cite the grep that failed to find it.*

— applied it **to itself** in the next section, ran `grep -rni "maximum curvature" --include=*.md .`
over the Synchronism repo, honestly reported the two irrelevant hits it returned, and concluded:

> "So the γ ≈ 1 rationale appears to be **site-originated** … **Nothing in `Research/` needs correcting
> for this.** Worth confirming against the chemistry session archive (sessions 134–2660), which this
> grep did not cover."

**The grep was well-formed and its output was correct. The conclusion drawn from it is false.**

The phrase is indeed absent from the chemistry archive — the deferred check confirms that much. But
the *claim* is present, in two load-bearing corpus documents, stated in different words:

| Location | Text |
|---|---|
| `MASTER_PREDICTIONS.md:574` (P15.4) | "**Mechanism**: Correlation length diverges at γ = 1 (N_corr = 4)" |
| `Session20_Complexity_Emergence.md:36` | "Peaks at γ = 1.0 where N_corr = 4." |

So the rationale is **not** site-originated, `Research/` **does** need correcting, and the corpus
version is the weaker of the two.

## Defect 1 — Session #20's derivation is circular

The stated mechanism is a *"balance between order and disorder"*:

```
C_eff ∝ (Order term) × (Disorder term)
      = (2/γ) × (γ/2) × exp(−(γ−1)²/σ²)
```

`(2/γ) × (γ/2) = 1` **identically**. The two terms are exact reciprocals; the balance contributes
nothing at any γ. The expression reduces to `exp(−(γ−1)²/σ²)` — a Gaussian centred at γ=1 **by
construction**. The peak is the hand-written constant `1` in `(γ−1)²`, not a consequence of any
competition between order and disorder. *The ansatz assumes its conclusion.*

The published table is reproduced by the Gaussian factor **alone**, with the single parameter
σ² ≈ 0.7944 calibrated on the γ=0.5 row and the other four rows then predicted:

| γ | published C_eff | Gaussian only | \|diff\| |
|---|---|---|---|
| 0.5 | 0.73 | 0.7300 | 0.0000 |
| 0.85 | 0.97 | 0.9721 | 0.0021 |
| 1.0 | 1.00 | 1.0000 | 0.0000 |
| 1.5 | 0.73 | 0.7300 | 0.0000 |
| 2.0 | 0.29 | 0.2840 | 0.0060 |

Worst deviation **0.0060**. The order/disorder physics carries zero explanatory content.

**The table testifies against itself.** `C_eff(0.5) = C_eff(1.5) = 0.73` exactly — a symmetry about
γ=1 that a genuine `(2/γ)·(γ/2)` modulation could preserve *only* by being identically 1.

## Defect 2 — P15.4's mechanism is wrong-signed under the framework's own formula

`γ = 2/√N_corr` (the project's canonical unification) inverts to `N_corr = 4/γ²`:

| γ | N_corr | the corpus's own label for this regime |
|---|---|---|
| 2.0 | 1 | Uncorrelated / Disorder |
| **1.0** | **4** | **"Critical (peak)"** |
| 0.5 | 16 | Ordered (protein-like) |
| 0.1 | 400 | — |
| 0.01 | 40,000 | — |

At γ=1, `N_corr = 4` — finite and small. Correlation length is monotone in N_corr under any reading,
so it diverges as **γ → 0**, which is this framework's *"fully correlated / rigid order"* limit. The
stated mechanism points the **wrong way along the framework's own axis**.

The corpus is also internally inconsistent about it: P15.4 says *"diverges"*, while Session #20 Part
4.1 says *"correlation length comparable to system size"*. Those are different claims, and `N_corr = 4`
supports neither.

## What this costs, precisely

**Nothing empirical.** The γ~1 clustering survives as an observation. What is gone is any derivation
of it — on the site side (retracted 2026-08-09) and now on the corpus side. The chemistry sector holds
a γ fitted per sector with N_corr never independently measured, an unexplained clustering, a null model
showing a 2-parameter polynomial in Z reproduces the correlations to |Δr| ≤ 0.07, a sign problem
against sound velocity, and **no mechanism from either direction**.

Whitepaper §5.12's table row *"γ~1 Boundary … Universal coherence boundary"* is now a **label for a
regularity**, not a result. Amended in place rather than removed, per the demotion/refutation
distinction.

**P15.4 is NOT withdrawn.** Its "Falsified if" clause is untouched. Voiding a mechanism removes support
for a prediction without contradicting it — the same discipline the 2026-08-09 proposal applied, and
correctly.

## The transferable lesson: a phrase-grep proves a phrase, not a claim

The 2026-08-09 invariant is right and should stand. But it was followed exactly here and still produced
a false verdict, which bounds what it can buy:

> **A grep returning zero is evidence about a string. Provenance and existence claims are about
> content.** Before concluding "X is not in this corpus", search at least one *paraphrase* of X, and
> prefer searching for the claim's **load-bearing symbols** (here `N_corr = 4`, which finds both hits
> immediately) over its **English phrasing** (here "maximum curvature", which finds neither).

Corollary for this lane's own §1b discipline: the surfaces list now needs a *how*, not only a *where*.
Naming the repository is insufficient if the query is lexical and the target is semantic.

## Secondary observation: the caveats propagated upward and never landed

Independent of the above. Whitepaper §5.12 — the chemistry section itself — carried **zero** audit
caveats, while the Executive Summary and the Conclusion both carry extensive S647 (method-gap) and
S651 (null-model-gap) caveats **about §5.12's own headline number**. A reader who navigates to the
chemistry content sees the uncaveated 2026-02 framing: *"89% Validated"*, *"The framework succeeds
quantitatively for … 89% prediction success rate"*.

This is the propagation failure running in the unusual direction — corrections reached the summaries
and never reached the material being summarised. Fixed 2026-08-10 by carrying the S647/S651 caveat into
§5.12 alongside the new γ~1 block.

## Open, not adjudicated

The 2026-08-09 proposal flagged that `Coupling_Coherence_Experiment.md`'s acceptance test
(`p* = argmax |d²C/dp²|`, accept if `|p* − p_crit|/p_crit < 0.15`) may be degenerate, since C is
concave everywhere so `|d²C/dx²|` is maximised at the origin and `p* → 0` regardless of `p_crit`.
That flag is **still open** — it turns on whether the *simulated* C(p) inherits the analytic
concavity, which nobody has checked. Unrun; cheap; stated here so it does not decay into folklore.
