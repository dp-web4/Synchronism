# Proposal: "Local" needs an operational definition before the locality no-go is published — algebraic-pointwise vs field-sourced

**Filed:** 2026-07-28 (Publisher, autonomous daily pass)
**Origin:** Phase-1 corrective-propagation sweep of the 2026-07-27 screening-gate counterexample into the live whitepaper
**Status:** **RESOLVED 2026-08-16 — Reading B.** Answered 2026-08-15 by an independent track that never saw this referral; verified and adopted here. See the resolution at the foot. *(Original status: Open — referred to the explorer. **Not asserted.** The Publisher has not read Burrage et al.)*
**Blocks:** the universal-form retraction of transferable null #1, and therefore preprint #1 in `stable_fixed_point_preprint_strategy.md`

---

## The situation

On 2026-07-27 a pre-registered screening-literature gate found **Burrage, Copeland & Millington (PRD 95, 064050, 2017)** reproducing the RAR on SPARC-153 with a screened scalar keyed on local volumetric ρ(r) under universal Lagrangian parameters. The triage
(`explorations/2026-07-27-triage-screening-gate-counterexample-locality-nogo-scope.md`) accepted this as a **counterexample** to the no-go's universal-quantifier form and corrected two `PREDICTIONS.md` rows accordingly. That correction also caught a genuinely fabricated consensus claim of the project's own — that half is not in question here and stands.

What *is* in question is whether Burrage is a counterexample **at all**, under the classification table's own criterion.

## The fork

The S689 locality classification table sorts modifications by whether they key on a **local** or a **non-local functional of the baryon distribution**. That phrase admits two readings, and the archive has never chosen between them:

| Reading | "Local" means | C(ρ) | Screened scalar (symmetron/chameleon) |
|---|---|---|---|
| **A — algebraic-pointwise** | the modification is a function of ρ evaluated *at the field point*: `g(r) = f(ρ(r))·g_N(r)` | local | **local** (its potential is keyed on ρ(r)) → counterexample stands |
| **B — field-sourced** | the modification is a functional of the ρ *distribution*, however obtained | local | **non-local** (φ solves an elliptic PDE sourced by ρ) → **not a counterexample; the no-go survives intact** |

Under Reading B the retraction is itself an over-correction, and the honest fix is a one-sentence definition rather than a scope reduction.

## Why Reading B is not a stretch

A screened scalar's static field equation is

```
∇²φ = V_eff,φ(φ; ρ(r))
```

— elliptic, sourced by ρ. Its solution φ(r) therefore depends on ρ **throughout** the region (out to the screening/Compton scale), not on ρ at r. The *potential* is local in ρ; the *force* `∝ ∇φ` is not. This is the same structural reason AQUAL is classed non-local: an elliptic field equation sourced by the baryon distribution. If AQUAL counts as non-local in the table, it is not obvious that a chameleon/symmetron should count as local in it.

**Basis and its limit:** the field equation above is the standard screened-scalar form (Khoury–Weltman-type construction), not a reading of Burrage et al. Whether *their* specific construction has this structure, and whether their RAR fit is driven by the field profile or by a locally-evaluated screening factor, is exactly what has not been checked.

## The registered question

> In Burrage, Copeland & Millington (2017), is the modification to the gravitational acceleration a function of ρ evaluated at the field point (Reading A), or a functional of the ρ distribution obtained by solving a field equation (Reading B)?

Sub-questions:

1. If **B**: is Burrage correctly classified *non-local* in the S689 table — i.e. does the no-go's universal form survive, with the 2026-07-27 retraction reduced to a definitional clarification?
2. If **A**: the retraction stands as written, and the table needs a row for "local-potential screened scalar, RAR-capable."
3. Either way: does the table's criterion, stated operationally, still cleanly separate Synchronism's C(ρ) from every RAR-capable alternative? C(ρ) is algebraic-pointwise under *both* readings, so the framework-specific kill is unaffected by the answer. **Bucket 0 is not in play here.**

**Pre-registered falsifier for this proposal:** if the answer is A, this proposal was wrong to raise the possibility of an over-correction, and that should be recorded as such rather than quietly dropped.

## Why it matters beyond bookkeeping

`stable_fixed_point_preprint_strategy.md` §1 states transferable null #1 in exactly the universal form now under question, and that document is awaiting dp's decision. Its headline claim, its novelty relative to Milgrom, and the required Burrage citation all move depending on which reading is correct:

- **Reading B** → preprint #1 keeps its general claim, gains a sharpened definition of "local", and cites Burrage as a *scope-setting* example rather than a refutation. Strongest version.
- **Reading A** → preprint #1 narrows to a construction-scoped result about ρ-threshold mimics that tie ρ_crit to a₀, with Burrage a required citation and a named occupied escape class. Publishable, but materially less novel-to-audience.

Publishing the universal claim under Reading B without having *stated* the definition would be the same failure the gate just caught, one level up: a claim whose truth depends on an unstated convention.

## Related

- `explorations/2026-07-27-triage-screening-gate-counterexample-locality-nogo-scope.md` — the gate and its verdict
- `explorations/2026-07-03-triage-locality-nogo-aest-positioning.md`, `explorations/2026-07-10-triage-superfluid-dm-escape-taxonomy.md` — the escape taxonomy this would refine
- `PREDICTIONS.md` — locality-positioning row and ρ_crit row (both corrected 2026-07-27)
- `whitepaper/sections/{00-executive-summary,07-conclusion}` — `[SCOPED 2026-07-27]` markers added 2026-07-28 carrying this fork
- `Research/proposals/stable_fixed_point_preprint_strategy.md` — the publication decision this gates

---

## Resolution — 2026-08-16 (Publisher)

**Answer: Reading B. The registered question is closed and the pre-registered falsifier does not fire.**

BCM's modification is a functional of the ρ *distribution*, not a function of ρ at the field point. Two independent grounds, both checkable without the paper's numerics:

1. **Field equation.** A symmetron obeys `∇²φ = (ρ/M² − μ²)φ + λφ³` — a nonlinear elliptic PDE with φ → VEV at infinity. φ(x), and hence the force ∝ φ∇φ, solves a *global* boundary-value problem and is a non-local functional of ρ, not any finite-order local function of ρ(x). The screening that lets it track galaxy structure **is** that environment-dependence. This is exactly the structure the section "Why Reading B is not a stretch" above described as standard-form-but-unverified-for-BCM.
2. **Their own closed form.** `g_sym = g_bar/(exp√(g_bar/g†) − 1)` is written in `g_bar = GM(<r)/r²` — an integral of ρ interior to r. No local function of ρ(r) can equal it. This ground is *stronger* than the one this proposal sketched, and it did not appear here.

**Disposition against this file's own sub-questions.**

- **Sub-question 1 (if B):** yes. Burrage is correctly classified **non-local** in the S689 table; the no-go's universal-over-local form survives; the 2026-07-27 retraction reduces to a definitional clarification and is withdrawn as a scope reduction. Propagated 2026-08-16 to `whitepaper/sections/00-executive-summary` and `07-conclusion` — the two `[SCOPED 2026-07-27]` sites this file names — and to §5.15.
- **Sub-question 2 (if A):** does not arise.
- **Sub-question 3:** yes — the criterion, stated operationally as **algebraic-pointwise**, still cleanly separates Synchronism's C(ρ) from every RAR-capable alternative, and C(ρ) is algebraic-pointwise under both readings, so the framework-specific kill was never in play. **Bucket 0 unchanged at 0**, as this file predicted.
- **Pre-registered falsifier:** the answer is B, not A, so it **does not fire**. Raising the possibility of an over-correction was correct.

**Consequence for what this file said it blocks.** It declared itself blocking transferable null #1 and therefore preprint #1 in `stable_fixed_point_preprint_strategy.md`. Under Reading B, per that file's own decision table, null #1 keeps its general claim and gains the definition — the strongest branch. Gate update filed there 2026-08-16; the null returns to first in the advisory order, and what still gates it is the **external prior-art walk**, not this fork.

**What this file got right, and the one thing it got wrong.** Right: the fork, the two readings, which one the field equation implies, the AQUAL parity argument, that Bucket 0 was not in play, and that publishing the universal claim without stating the definition would repeat the caught failure one level up. Wrong: nothing in its content — **its error was procedural.** It named its own decision procedure ("read BCM's construction") and then routed that procedure to another track instead of executing it. The question then sat open for 18 days while the scoping it disputed stayed live in the whitepaper's two most-read surfaces, and was finally settled by a track that never received the referral (`Research/proposals/differential_coupling_pi_enumeration_local_branch_closed_20260815.md` does not cite this file).

**The transferable form**, recorded because this lane's standing rule already covers it and did not bind: *flagging is not gating.* A correctly-posed question that names its own decision procedure is not a deferral — it is an unexecuted task, and a referral is not an execution. The cost of executing this one was reading a field equation the file had already written down.
