# Governance Arc — Phase 12: Recovery, Split Brain, and Assurance-Relative Continuity

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — architecture toy / governance result  
**Code:** [`simulations/governance_phase12_split_brain_recovery.py`](../simulations/governance_phase12_split_brain_recovery.py)  
**Scope:** recovery / continuity evidence; not a security proof, protocol standard, or physics claim

## Why this phase exists

The continuity-event profile identified **recovery** as the hardest identity transition.

A restore from backup can create a new runtime with impeccable historical state while an old runtime may still exist or later reappear.

Thus:

\[
\text{state restoration}
\neq
\text{exclusive identity restoration}.
\]

The unresolved question is:

> **How should a relying party reason while descent is already proven but exclusivity is not yet known?**

That is the split-brain interval.

---

# Construction

Use a deliberately minimal stochastic timeline.

At

\[
t=0
\]

original runtime/token \(A\) crashes.

Four random delays are then sampled.

### Recovery delay \(R\)

A backup is restored and descendant \(B\) becomes available.

At \(R\), the system can already provide positive provenance:

\[
A\preceq B.
\]

### Registry partition delay \(P\)

The authoritative society / registry is temporarily unreachable or unable to finalize the cutover.

### Old-runtime resurrection delay \(O\)

With probability

\[
p_{res}=0.60,
\]

the old runtime \(A\) later comes back.

### Fence / finality delay

After both \(B\) exists and registry connectivity is restored, a small additional delay represents:

- witnessed fencing;
- registry state commitment;
- revocation / retirement of the old recognized continuation;
- final continuity receipt.

Call completion time \(F\).

Thus:

\[
F\ge\max(R,P).
\]

The baseline exponential means are:

- restore: 8 time units;
- partition: 30;
- old resurrection: 18 conditional on resurrection;
- finalization after prerequisites: 5.

These numbers are illustrative.

---

# Two relying parties

The model deliberately gives two relying parties different stakes / evidence requirements.

## A1 / low-stakes profile

Accepts \(B\) for a reversible / low-consequence act once signed descent and recovery evidence exists at

\[
t=R.
\]

It does **not** claim exclusivity is proven.

It merely tolerates the ambiguity.

## A2 / high-stakes profile

Requires scoped registry / fencing evidence supporting exclusive continuation.

It waits until

\[
t=F.
\]

Neither profile is declared universally correct.

They represent different evidence requirements for different consequences.

---

# Baseline result

Using 200,000 deterministic-seed Monte Carlo episodes:

| quantity | result |
|---|---:|
| mean restore time / A1 acceptance | 7.976 |
| mean fence-finality / A2 acceptance | 36.708 |
| mean A2 − A1 assurance latency | 28.731 |
| median A2 − A1 assurance latency | 19.207 |
| probability of actual dual-live interval before finality | 0.4472 |
| mean dual-live duration, conditional on occurrence | 28.015 |
| probability old runtime resurrects only after finality | 0.1542 |
| probability old runtime never resurrects | 0.3986 |
| probability partition outlasts restore | 0.7898 |

So the lower-assurance relying party gets usable continuity evidence roughly

\[
29
\]

time units earlier on average.

But in about

\[
44.7\%
\]

of episodes, an actual period occurs before finality during which both \(A\) and \(B\) are live.

When such a period occurs, its mean duration is about

\[
28.0
\]

time units.

This is the basic assurance/latency tradeoff.

---

# The key correction: reality and evidence are different state variables

The important result is not the particular percentages.

It is the distinction between:

### Ontic / runtime state

How many continuations are actually alive?

### Evidential state

What continuity proposition does the available evidence support for this witness?

Those can differ.

For example:

### Case A — exclusive in reality, unknown in evidence

The old runtime never comes back.

After \(B\) is restored, it may already be the only live continuation in fact.

But until fencing / registry evidence is finalized, a high-assurance witness cannot establish that negative claim.

### Case B — branching in reality, still unknown to a witness

The old runtime quietly resurrects during the partition.

Both \(A\) and \(B\) are live.

A relying party with only local recovery evidence may not yet know that.

### Case C — old runtime returns after fence

Old \(A\) can still produce signatures that verify cryptographically.

But the authoritative registry has already fenced / invalidated it.

Thus:

\[
\boxed{
\text{signature-valid}
\neq
\text{currently authorized continuation}.
}
\]

That is a particularly important governance distinction.

---

# Continuity evidence should be at least tri-state

A simple boolean

```text
same_identity = true|false
```

cannot express the recovery interval honestly.

A better evidence state is at least:

\[
\boxed{
\{
\text{exclusive-supported},
\text{branching-supported},
\text{unknown}
\}
}
\]

Potentially with assurance and scope attached.

### `exclusive-supported`

Available evidence supports one recognized live continuation within a named scope at the claimed assurance.

### `branching-supported`

Available evidence positively supports multiple live / valid descendants.

### `unknown`

Descent may be known, but available evidence cannot yet establish singularity or branching strongly enough.

Unknown is a first-class honest state, not an error condition.

---

# Low assurance does not mean false

The A1 relying party accepts \(B\) at restore time while exclusivity is still unknown.

That should not be described as "incorrect."

For a reversible, low-value action the party may rationally decide:

> valid descent + current key control is enough evidence for this act; split-brain ambiguity is an acceptable risk.

The A2 party may decide:

> the act is irreversible / high-value; unknown exclusivity is not enough.

So:

\[
\boxed{
\text{same evidence state}
+
\text{different stakes}
\rightarrow
\text{different valid relying-party actions}
}
\]

This is exactly the Web4 principle that evidence is inspectable while interpretation remains contextual.

---

# Partition duration controls both assurance cost and ambiguity exposure

The simulation sweeps mean registry/network partition duration while holding the other baseline parameters fixed.

| mean partition | mean A2−A1 wait | probability dual-live before finality | mean dual-live duration if present |
|---:|---:|---:|---:|
| 2 | 5.400 | 0.2814 | 5.144 |
| 5 | 6.898 | 0.3081 | 6.205 |
| 10 | 10.555 | 0.3511 | 9.627 |
| 20 | 19.270 | 0.4081 | 18.474 |
| 30 | 28.623 | 0.4452 | 27.861 |
| 60 | 57.986 | 0.5047 | 57.482 |
| 120 | 117.277 | 0.5489 | 116.792 |

As authoritative connectivity becomes slower:

- high-assurance continuity becomes more expensive in latency;
- the window in which low-assurance action occurs under unresolved exclusivity grows;
- old-runtime resurrection is increasingly likely to occur before finalization.

Again, the values are toy-specific.

The durable result is the coupling:

\[
\boxed{
\text{weaker / slower continuity evidence}
\rightarrow
\text{larger unresolved identity interval}.
}
\]

---

# A continuity state machine

The recovery process can be represented more honestly as evidence transitions.

Before restore:

```text
no-current-descendant-evidence
```

After signed recovery at \(R\):

```text
descent-supported / exclusivity-unknown
```

If a competing live descendant is witnessed:

```text
branching-supported
```

After authoritative fence / registry finality at \(F\):

```text
exclusive-supported(scope=S, assurance=A2)
```

This state machine says what evidence supports.

It does not pretend to be an omniscient statement about all possible hidden copies everywhere.

---

# Exclusivity is scoped, not universal

At finality, the toy assumes an authoritative registry can say:

> within this society / registry scope, \(A\) is retired and \(B\) is the sole recognized active continuation.

That is strong evidence.

It is still not a proof that no unauthorized clone exists anywhere.

So the proper proposition is:

\[
\operatorname{ExclusiveContinuation}
(A,B\mid S,A_k,t),
\]

where:

- \(S\) = authority / registry scope;
- \(A_k\) = assurance profile;
- \(t\) = evidence time / freshness.

This is another place where MRH is doing real conceptual work.

---

# MRH interpretation

The same physical recovery has different relevant distinctions for the two relying parties.

### Low-stakes MRH

The distinction

\[
\text{exclusive-supported}
\quad vs\quad
\text{unknown}
\]

may be irrelevant to the current reversible act.

The witness safely quotients it out.

### High-stakes MRH

That distinction is load-bearing.

It must remain explicit.

So the assurance profile is not merely "more security."

It changes which identity distinctions remain relevant to the action.

That is almost exactly the refined MRH definition produced by the earlier arc:

> **an MRH is a witness/task/stakes-indexed boundary over distinctions that can safely be discarded.**

---

# A new separation: continuity truth, continuity evidence, continuity policy

Phase 12 suggests keeping three layers explicit.

## Continuity truth / runtime state

What actually exists in the modeled world?

## Continuity evidence

What claims are supported by the witnessable records available in this MRH?

## Continuity policy

What evidence does this relying party require for this act?

Then:

\[
\text{action}_R
=
F_R(
\text{evidence},
\text{stakes},
\text{freshness},
\text{assurance}
).
\]

Web4 should mainly standardize the **evidence surface**, not one universal \(F_R\).

---

# Recovery implies an important hub rule

A hub should not treat possession of:

- an old key;
- an old LCT document;
- a valid signature;
- a backup state;
- even valid ancestry

as sufficient by itself to prove exclusive current authority.

For high-consequence acts, it may require a fresher continuity proof:

- current registry state;
- fencing / revocation evidence;
- continuity receipt;
- assurance profile;
- policy compatibility.

This is the identity analogue of checking revocation status rather than accepting a cryptographically valid but revoked certificate forever.

---

# Phase-12 verdict

**PASS — the recovery model gives operational content to scoped, assurance-relative continuity.**

What survived:

1. valid descent and exclusive continuation are different propositions;
2. recovery naturally creates an interval where descent is known but exclusivity is unknown;
3. actual split brain can exist during that interval;
4. a system can be exclusive in reality before the witness has enough evidence to prove it;
5. cryptographic validity can survive after current authority has been fenced;
6. continuity evidence should represent `unknown` rather than force a boolean answer;
7. low- and high-stakes relying parties can legitimately act differently on the same evidence;
8. exclusivity is scoped by authority, assurance, freshness, and therefore MRH.

What did not happen:

- no claim that A1/A2 are the correct names or thresholds for a standard;
- no claim that the toy delay distributions resemble a deployed hub network;
- no claim that one registry can establish global nonexistence;
- no security proof;
- no physics result.

The durable statement is:

\[
\boxed{
\text{descent can be proven before exclusivity can be known.}
}
\]

---

# Next question — provenance has an MRH too

A mature provenance DAG can become enormous.

A relying party should not need the entire historical graph for every decision.

The natural next question is therefore:

> **What is the smallest provenance/evidence subgraph that preserves the continuity, authority, policy-compatibility, and trust facts relevant to one relying party's current act?**

This is the same MRH operation applied to history itself.

Potentially:

\[
\text{full provenance DAG}
\rightarrow
\text{relevant proof subgraph}
\rightarrow
\text{decision}.
\]

The proof subgraph would retain load-bearing distinctions such as:

- the path from the relying party's admission anchor to the current token;
- relevant fork / merge points;
- fencing evidence where exclusivity matters;
- constitutional amendments touching relevant policy behavior;
- the unique evidence receipts actually consumed by trust derivation.

Everything else may be safely omitted from that particular proof bundle.

That looks like the next useful formal step.
