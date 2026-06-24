# Phase-13 (door #2/LIV) — the custodial escape is closed: absolute time ⊥ Standard-Model couplings without dim-4 LIV (2026-06-23)

**Status:** `[ACTIVE-MRH]` — repays the Phase-12 target ("compute the radiatively-induced dim-4 LIV,
or exhibit the custodial symmetry that suppresses it"). **Result: the only custodial symmetry that
can forbid dim-4 LIV is Lorentz *boost invariance* — and the framework's defining commitment
(absolute time / a single global clock) is *precisely the absence* of boost invariance. Phase-6's
"emergent Lorentz" is an external-leg, tree-level, soft-pattern (N≫1) effect, **not a symmetry of
the action**, so it cannot constrain loop corrections (which sample N~1 cutoff modes). A toy
computation confirms the boundary: radiative dim-4 LIV *decouples only when the loop is UV-finite*
(super-renormalizable), and goes *marginal / non-decoupling* exactly at D=4 with marginal couplings
— the dimension of real physics. So a substrate that (a) keeps absolute time and (b) contains the
Standard Model's marginal couplings generically produces dim-4 LIV at O(g²/16π²), bounded at
~10⁻¹⁸–10⁻²² — i.e. refuted-or-fine-tuned. Door #2 is the framework's *sharpest* falsification
risk, and the risk is tied to its most defining commitment, not a tunable parameter.**
**Sim:** [`simulations/phase13_radiative_dim4_liv.py`](../simulations/phase13_radiative_dim4_liv.py) · result: `simulations/results/phase13_radiative_dim4_liv_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The question Phase-12 left

Phase-12 showed door #2's *time-of-flight* channel is untestable, but the **dim-4** channel
(species-dependent limiting speed / c-anisotropy) is even-in-k (so the even-k symmetry doesn't
forbid it) and *reachable* (~10⁻¹⁸–10⁻²²). It left one escape: maybe a custodial symmetry suppresses
dim-4 LIV below reach — candidate, Phase-6's emergent Lorentz invariance (the PN barrier of an
extended kink vanishes as exp(−N)). This phase tests that escape.

## Why Phase-6's mechanism is the wrong kind of suppression

Phase-6 measured the **Peierls–Nabarro barrier of a static, extended kink** — a *free-field,
on-shell* quantity. It vanishes as exp(−N) because a pattern spanning **N≫1** cells **averages
over** the lattice; the suppression lives in the *external pattern size N*. Likewise Phase-5's
moving-soliton boost test is an on-shell, propagating-pattern (N≫1) statement. These establish
emergent Lorentz invariance **for soft external patterns** — they are *not* a symmetry of the
microscopic action.

Radiative dim-4 LIV is a **loop** effect. The loop integral runs over virtual modes at *all*
scales, and the dim-4 coefficient is sourced where the integral has its UV weight — at the cutoff,
where a mode spans **N~1** cells. There is no large-N averaging there. An external-leg/soft-pattern
restoration of Lorentz invariance (exp(−N), N≫1) **cannot constrain a loop sourced at N~1.** So
Phase-6 is structurally incapable of being the custodial symmetry. (This is the Collins–Perez–
Sudarsky–Gambini–Pullin point, *PRL* **93**, 191301, 2004: tree-level/emergent Lorentz invariance
does not survive radiative corrections without a *symmetry of the action* forbidding the LIV
operators.)

## The toy: decoupling ⟺ UV-finiteness; marginal LIV ⟺ marginal couplings

Minimal boost-violating setup faithful to the framework: continuous time ω (absolute clock) + a
**spatial lattice** (spacing a), scalar propagator `G=1/(ω²+K(k)²+m²)`, `K(k)=(2/a)sin(ka/2)`, one
φ³ interaction, 1-loop self-energy bubble. The renormalized inverse propagator is
`Γ ≈ const + A·p0² + B·p²`; Lorentz invariance needs `A=B`, so `c_LIV ≡ B/A−1` is the dim-4 LIV.

The honest first result was a **surprise that had to be diagnosed**: in the bare 2D φ³ toy, `c_LIV`
*decouples* (`|c_LIV| ∝ (m a)^1.70 → 0`). Taken at face value that would *support* "untestable." But
2D φ³ is **super-renormalizable** — its bubble is UV-finite — so it *cannot* exhibit the Collins
percolation by construction; the decoupling is an artifact of UV-finiteness, not protection. Tuning
the loop's UV sensitivity (radial phase-space weight `|k|^{d_s−1}`, a power-counting proxy for going
to D=4) makes the boundary explicit:

| effective spatial dim `d_s` | spacetime D | loop UV behaviour | decoupling power `p` (`|c_LIV|∝(ma)^p`) | verdict |
|---|---|---|---|---|
| 1 | 2 | finite (q⁻³) | **1.70** | decouples (safe) |
| 2 | 3 | finite (q⁻²) | 0.54 | marginal |
| 3 | **4** | **log-divergent** | **0.10 ≈ 0** | **non-decoupling (Collins)** |

**Decoupling ⟺ UV-finiteness (super-renormalizable). Marginal, non-decoupling dim-4 LIV appears
exactly at D=4 with marginal couplings** — the dimension of real physics. (The proxy is not a full
4D lattice loop; it is a power-counting knob that *reproduces the known Collins result*, used to make
the structural boundary visible. The decisive content is the symmetry argument above; the scan
illustrates it.)

## The tension, named

The framework needs **marginal couplings** to contain the Standard Model (gauge α~1/137, strong
g_s~1, Yukawas, the Higgs quartic — all dimensionless). On the **absolute-time** substrate, those
couplings sit in the D=4 / marginal class where dim-4 LIV does *not* decouple. The only thing that
could forbid it — Lorentz **boost invariance** as a symmetry of the action — is *exactly what
absolute time gives up*. So:

> **Absolute time (the framework's defining "one move") and a Standard-Model-containing interacting
> completion cannot coexist with sub-10⁻²² Lorentz invariance, absent a custodial symmetry the
> framework has not identified and whose natural form (boost invariance) its ontology forbids.**

Both horns are load-bearing and both are taken seriously by the program, so this is a genuine
foundational tension, not a parameter choice:
- Keep absolute time + SM couplings ⇒ generic dim-4 LIV at O(g²/16π²) ⇒ **refuted-or-fine-tuned**
  by table-top cavity/Hughes–Drever/clock nulls (10⁻¹⁸–10⁻²²).
- Suppress dim-4 LIV by a real custodial symmetry ⇒ you need boost invariance **in the action** ⇒
  you have given up absolute time, i.e. the framework's distinctive ontology.
- Suppress it by fine-tuning ⇒ the part-in-10²² tuning is exactly the "additional fine-tuning
  problem" Collins et al. named — the program would have to own it explicitly.

This is the deepest answer the arc has reached to "what is the framework protecting?" (the standing
tension #5): **it is protecting absolute time**, and absolute time is simultaneously its identity and
its sharpest exposure to refutation by existing precision data.

## What this does to the ledger / the doors

- Door #1 (galactic) = MOND cage (Phase-11). Door #2 time-of-flight = untestable (Phase-12, leak
  sealed). Door #2 **dim-4 = falsifiable, and the custodial escape is now closed** (Phase-13):
  not "untestable," but "the framework's sharpest falsification risk, structurally tied to absolute
  time."
- This is **not** a Bucket-0 confirmation and **not** (yet) a Bucket-2 refutation — the framework
  hasn't specified its interactions, so the coefficient isn't computed. It is a **named tension and a
  sharp obligation**: any interacting completion must either exhibit the custodial mechanism
  (and explain how it coexists with absolute time) or accept refutation. That obligation is the real
  loan-repayment target now — far sharper than "test B7's Umklapp."

## Honesty / caveats

- The symmetry argument (boost invariance is the custodial symmetry; absolute time forbids it;
  Phase-6 is external-leg/tree-level not a loop symmetry) is the rigorous core and needs no numerics.
- The toy is 1+1D φ³ with a power-counting proxy for higher D; it *reproduces* the established
  Collins decoupling-vs-marginal boundary, it does not independently prove it. A full 3-lattice +
  continuous-time 4D loop (or a genuine marginal-coupling theory) would be the next-order rigor, but
  the conclusion is a textbook EFT result, not a novel numerical claim.
- "Refuted-or-fine-tuned" is conditional on the framework's (unspecified) interactions containing
  marginal couplings — which it must, to contain the SM. I claim a structural obligation, not a
  completed refutation. **Bucket 0 unchanged (0).**
- Exotic escapes not in the framework (supersymmetry can relate/forbid some LIV coefficients; a
  bespoke discrete custodial symmetry) are not excluded — but none is part of Synchronism, and
  invoking one is itself the unpaid obligation.

## So what

The arc's last open door looked like a safe exit ("untestable, symmetry-protected"). Walking through
it, the symmetry that was supposed to protect it turns out to be the very thing the framework's
identity gives up. Door #2 is not where Synchronism hides from falsification — it is where its
**deepest commitment (absolute time) collides with the precision frontier**, and the collision is
structural. That is the most falsifiable, most Synchronism-specific statement the whole arc has
produced — and it points at refutation, which is a real result, not a comfortable one.
