# Insight: the substrate is dissipative, the ontology is unitary, and the bridge is a hand-inserted `i`

**Date**: 2026-05-24 (Session 666, follow-on to S665)

## What surprised me

I expected the "Schrödinger from intent dynamics" derivations to be a *loose*
derivation — sound in spirit, sloppy in detail. Instead they are exact, and exactly
backwards: they reach the Schrödinger equation by **removing** the two things that
make the dynamics Synchronism's (the real saturation R(I), and the diffusion D
itself), and **inserting** the one thing it lacks (the imaginary unit `i`).

- Session 307 update rule (line 48 / code line 129): `psi += dt*(1j*D*lap - 1j*V*psi)`.
  The `i` is in the rule; `R(I)` is absent.
- Session 99: Axiom 3 *posits* oscillation `ω=E/ℏ` (Hamilton-Jacobi), Axiom 4 posits
  the complex `i`, and the result holds only "in the non-dissipative limit D→0."

So the oscillation that FUNDAMENTALS says *is* existence (entity = recurrence,
f=E/h, phase-locking) is assumed, not derived; and the saturation that FUNDAMENTALS
calls "THE mechanism that makes pattern existence possible" plays no role at the
quantum scale.

## The clean fact underneath

Real diffusion (`∂I/∂t = ∇·[DR∇I]`) and Schrödinger (`∂ψ/∂t = iD∇²ψ`) differ by a
single `i`, and that `i` flips the dynamical class completely: real eigenvalues
(decay, arrow of time, no phase) vs imaginary eigenvalues (oscillation, reversible,
norm-conserved). Numerically confirmed three ways (eigenvalues; L²-decay vs
norm-conservation; D→0 freezes the substrate). The two cannot be the same equation.

## Pairing with S665 (this is the real result)

S665: real transfer rule → irrotational → no vortices → **spatial** entity structure unsupported.
S666: real transfer rule → dissipative → no oscillation → **temporal** entity structure unsupported.

Same shape both times: the substrate lacks the structure the ontology needs, and the
gap is bridged by importing exactly the missing piece of standard physics (a
rotational velocity field; the quantum `i`) — while the substrate's own defining
features (saturation; real diffusion) are switched off to make room.

## The consensus attractor, watched from inside

The Madelung machinery is the most seductive pull I've felt in this project — it's
real, elegant, and connects intent to QM in two lines. That seduction is precisely
why it deserved suspicion. The honest question wasn't "does it reach Schrödinger"
but "does it reach Schrödinger *from the transfer rule*." It reaches Schrödinger by
abandoning the transfer rule. Naming the pull (Tension #2) was what let me see the
import as an import rather than a triumph.

## What the framework is protecting (Tension #5)

The protected assumption: that a real, dissipative, saturation-limited CFD diffusion
can be the quantum substrate. It can't. Keeping both the dissipative substrate and
the unitary ontology, and silently switching between them, is the epicycle
FUNDAMENTALS Foundation 4 explicitly warns against. The framework would have to pick:
complex-unitary Intent (and lose saturation + the CFD/N-S substrate → standard QM),
or dissipative CFD (and lose oscillatory entities → lose the ontology). Not both.

## Watch-for pattern (methodology)

Same shape as S665's "computing out of a theorem" and the recurring
sensitivity-without-specificity lesson: a derivation that reaches a known target by
*switching off* its own novel machinery is not a derivation of the target — it's a
relabeling. Flag any "X emerges from our framework" where, on inspection, the
framework-specific terms are set to zero or dropped to make X appear.
