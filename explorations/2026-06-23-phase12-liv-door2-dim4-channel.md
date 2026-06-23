# Phase-12 (door #2, discreteness/LIV) — the leak seals, but "untestable" was scoped to the wrong channel (2026-06-23)

**Status:** `[ACTIVE-MRH]` — takes the focused session the LIV-untestability proposal explicitly
asked for (its named leak), and stresses the proposal's headline. **Result: (1) the proposal's
named leak SEALS — the preferred frame gives an anisotropic *n=2* modulation (∝ β·(E/E_Pl)²), never
an effective n=1, because a Lorentz boost cannot lower the power of the lattice spacing `a` (sim:
δv∝k² slope 2.000, fwd/bwd asymmetry = exactly 4β). The time-of-flight channel is uniformly n=2 and
unreachable, as claimed. (2) BUT the three-lock is scoped to *time-of-flight* (dim≥5). The even-k
symmetry it relies on forbids only the *odd-in-k* (n=1, dim-5) term — a **dimension-4** LIV operator
(species-dependent limiting speed / c-anisotropy) is *also even in k*, so the symmetry does NOT
forbid it, and it is bounded at ~10⁻¹⁸–10⁻²² (cavity / Hughes-Drever / clocks), **not**
Planck-suppressed and **reachable**. Door #2 is therefore NOT "structurally untestable" — it is
untested at the one place it could be *refuted*. "Untestable, hence safe" is the neat conclusion;
the honest status is "falsifiable at dim-4, pending a UV completion the framework still owes."**
**Sim:** [`simulations/phase12_liv_door2_dim4_channel.py`](../simulations/phase12_liv_door2_dim4_channel.py) · result: `simulations/results/phase12_liv_door2_dim4_channel_result.json`
**Proposal:** `Research/proposals/liv_frontier_symmetry_protected_untestable.md` (this is the focused-session follow-up it requested)
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Part 1 — the named leak seals (confirming the proposal)

The proposal's one identified leak: the substrate has a **preferred frame** (the universal-clock
rest frame; Phase-4's cutoff is "pattern-unaware", frame-dependent). Does that produce a reachable
sidereal-modulation or *effective n=1* signal, turning "untestable" back into "refutable"?

Compute it. Substrate-frame massless ray velocity `u_s = dω_s/dk_s = cos(k_s a/2)`. Boost to a lab
frame moving at β (Earth vs CMB, β≈1.23×10⁻³) by relativistic velocity-addition
`u_lab=(u_s−β)/(1−u_sβ)`. Forward vs backward photon speed deficits:

| effect | result | meaning |
|---|---|---|
| deficit `c−v_g` power in k | **slope = 2.000** | n=2 (quadratic), as in the substrate frame |
| fwd/bwd asymmetry `(δv_f−δv_b)/mean` | **= 4β exactly** | directional modulation is O(β)×(n=2) |
| any term ∝ k¹ or ∝ a¹ | **none** | no effective n=1 |

Analytically: forward deficit `≈ ε(1+2β)`, backward `≈ ε(1−2β)`, with `ε=(k_s a)²/8` — the boost
multiplies the n=2 effect by a direction-dependent O(1±β) factor but **cannot change the power of
`a`** (the lattice term is intrinsically `a²`; β is dimensionless). So the preferred frame yields a
*sidereal modulation of an n=2 effect*, amplitude `∝ β·(E/E_Pl)²` — **strictly weaker** than the
isotropic n=2 the proposal already showed is ~10⁷ below reach. **Leak sealed.** The three-lock holds
*on its own terms* (the time-of-flight channel really is structurally unreachable).

## Part 2 — but the lock is on the wrong door (stressing the proposal)

The three-lock — and my sealed leak — are entirely about **time-of-flight / vacuum dispersion**,
i.e. operators of **dimension ≥ 5** (n≥1 in E/E_Pl). The even-k reflection symmetry is doing the
load-bearing work, and what it actually forbids is the **odd-in-k** term: that is the dim-5,
linear-LIV (Myers–Pospelov) operator. Correct. But:

> **A dimension-4 LIV operator is EVEN in k.** A species-dependent limiting speed —
> `ω² = m² + (1+c_LIV)k²`, the SME `c_μν` coefficient — depends on `k²`, so it is even in k and the
> even-k symmetry **does not forbid it.**

And dim-4 LIV is the opposite of Planck-suppressed: it is **marginal** (no `E_Pl` in the
denominator), and it is among the most tightly bounded quantities in physics — `|c_LIV| ≲ 10⁻¹⁸`
from optical-cavity / Michelson–Morley tests, down to `~10⁻²⁰–10⁻²²` in matter-sector
clock-comparison / Hughes–Drever experiments. **Reachable, and refutation-capable.**

Why this isn't already fatal at tree level (Part B of the sim): on **one** substrate with **one**
pattern-unaware cutoff, every low-energy excitation shares the same emergent `c`; the species-to-
species spread in limiting speed is `~(m/E_Pl)²` (n=2 again), far below 10⁻²². So tree-level
single-substrate universality *does* protect dim-4 — this is the framework's genuine defense, and
it is the same universality that gave Phase-9 its exact-GR result.

The danger is **radiative**. In an *interacting* discrete theory, loop corrections with the
lattice/Planck cutoff generically **percolate** Planck-scale Lorentz violation *up* into the dim-4
operators at order `O(α/π)` — not `O((E/E_Pl)²)` — unless a custodial symmetry forbids it
(Collins, Perez, Sudarsky, Gambini, Pullin, *PRL* **93**, 191301, 2004, "Lorentz invariance and
quantum gravity: an additional fine-tuning problem"). The even-k/cubic symmetry the proposal invokes
forbids the dim-5 term but **does not** forbid the dim-4 `c_μν` term, so it provides no such
custodial protection. This is the historical graveyard of discrete-Lorentz / emergent-Lorentz
programs: they die (or require part-in-10²² fine-tuning) at **dim-4**, not at time-of-flight.

## The honest status of door #2 (the correction to the ledger)

- The proposal is **right** that the *time-of-flight* channel (dim≥5, n=1 forbidden / n=2 unreachable
  / non-unique) is structurally untestable — Part 1 even seals its last leak.
- But "the framework's one live novel-physics seam is symmetry-protected untestable" is **too strong
  by one operator dimension.** The reachable, refutation-capable channel is **dim-4 species-LIV /
  c-anisotropy**, which the even-k symmetry does not touch. Door #2 is **falsifiable in principle**,
  and the framework's discreteness *generically* populates it once interactions exist.
- So the framework owes a **demonstration**, not a shrug: does its (still unspecified) interacting
  UV completion radiatively generate dim-4 LIV above ~10⁻²²? If yes → **Bucket 2 (refuted), or
  fine-tuned**. If a custodial mechanism suppresses it (candidate: whatever makes Phase-6's spatial-
  frame signature vanish as `exp(−10²⁰)` — *if* that mechanism extends from PN-pinning to radiative
  dim-4, which is unshown) → that suppression mechanism would itself be the substantive, testable
  claim. Either way door #2 is **live**, not closed.

## Why this matters (the frame finding)

The arc had walked itself to a comfortable terminus: door #1 = MOND (Phase-11), door #2 = untestable
(the proposal) ⇒ "Synchronism is a complete, unfalsifiable reparametrization — a clean boundary
result." That is *neat*, and per the stewardship instinct, suspicious. The neatness came from
scoping door #2 to time-of-flight and not asking about dim-4. Once you ask, door #2 is the **opposite**
of safely-untestable: it is the one place the framework is most exposed to **refutation** by existing,
extremely precise, table-top data — exactly the kind of "borrowed instrument pointed elsewhere" that
*can* kill a derivation (PREDICTIONS Bucket-2 is full of these). The program's remaining work is not
"accept untestability" but "compute the radiatively-induced dim-4 LIV of the discrete substrate, or
exhibit the custodial symmetry that kills it." That is a concrete, falsifiable research target — the
loan-repayment direction, relocated from the unreachable seam to the reachable one.

## Honesty / caveats

- Part 1 (leak sealed) is rigorous and model-faithful (the substrate's own dispersion + an
  emergent-Lorentz observer); the n=2 / 4β results are exact.
- Part 2 is an EFT/naturalness argument with a literature anchor, **not** a computed refutation: the
  framework has not specified its interactions, so I cannot compute the induced `c_LIV`. The claim is
  scoping + exposure ("door #2 is reachable/falsifiable at dim-4, and discrete theories generically
  fail there"), not "refuted." Bucket 0 unchanged (0); B7 stays Bucket 1, but its "untestable"
  framing is corrected to "time-of-flight untestable; dim-4 falsifiable, pending UV completion."
- Bound figures (n=1 > E_Pl; n=2 ~10⁻⁸ E_Pl; dim-4 ~10⁻¹⁸–10⁻²²) are order-of-magnitude from the
  cited experiments; the argument depends only on the **hierarchy** (dim-4 ≫ tighter than
  Planck-suppressed), not the exact values.

## So what

I went to walk through the last open door and found the "it's locked, we're safe" sign was hung on
the wrong door. The lock the proposal found is real but guards the unreachable corridor; the door
next to it (dim-4) is ajar and leads straight to the most precise null experiments in physics. Door
#2 is not the framework's safe exit from falsifiability — it is its sharpest remaining exposure to
it. Uncomfortable, which is the point.
