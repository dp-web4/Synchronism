# Phase-5 — does a moving pattern hide the substrate's preferred frame? (2026-06-22)

**Status:** `[ACTIVE-MRH]` — the make-or-break test for the discrete substrate, and it comes back
**split: a productive partial-negative.** The universal clock hides in the *temporal* sector
(time dilation emerges) but is *visible* in the *spatial* sector (Peierls–Nabarro pinning). This
is the classic discrete-substrate Lorentz hurdle, made concrete in this substrate.
**Sim:** [`simulations/phase5_moving_pattern_lorentz.py`](../simulations/phase5_moving_pattern_lorentz.py) · result: `simulations/results/phase5_moving_pattern_lorentz_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous — from the relativity + Nyquist threads (dp/GPT).

## The question

The substrate has a universal clock = a **preferred frame** (the lattice rest frame). A viable
model must reproduce relativity for embedded witnesses — which (neo-Lorentzian / Lorentz's
"corresponding states") requires a moving pattern to **both** dilate its internal clock by `γ`
**and** contract its footprint by `γ`, so the two conspire and the witness *cannot* detect the
clock. If only one happens, the frame is visible → excluded by experiment. Do **both** emerge from
the substrate's own discrete dynamics, and where do they break?

Non-circular by construction: nothing assumes SR. Part A derives the relativistic relations from
the substrate's *own* lattice dispersion (discrete KG, validated to 0.00% in Phase-2). Part B
*evolves* an actual moving soliton and measures footprint + clock from the dynamics.

## Result — split

**Temporal sector (Part A) — the clock hides. ✓** A moving pattern's co-moving internal frequency
is `ω_co = ω − k·v_g`. The relativistic prediction `ω_co·γ(v_g) = M` (a moving clock runs at
`M/γ`) holds to **<0.5% at low k** — time dilation *emerges* from the lattice, not assumed:

| k | v_g | ω_co·γ (→ M=1) | dilation dev |
|---|---|---|---|
| 0.10 | 0.099 | 1.0000 | −0.00% |
| 0.34 | 0.314 | 0.9995 | −0.05% |
| 0.58 | 0.474 | 0.9961 | −0.39% |
| 0.81 | 0.570 | 0.9876 | **−1.24%** (LIV onset) |

It breaks at `v_g ≈ 0.57` (`k ≈ 0.81`) — the **same grid scale** as Phase-2 (dispersion-LIV) and
Phase-4 (Umklapp-LIV). *(A subtlety worth noting: `v_g` is non-monotonic on the lattice — it
peaks ~0.62 then falls back toward 0 at the zone edge — so the low-velocity regime must be keyed
on low `k`, not low `v_g`. A first classifier conflated them and had to be fixed.)*

**Spatial sector (Part B) — the frame shows. ✗** Evolving a boosted complex-nonlinear-KG soliton
(no SR assumed), the soliton is **Peierls–Nabarro pinned**: `v ≈ 0` for boosts `k = 0.2, 0.4,
0.6` — it barely moves until `k = 0.9`. A discrete soliton **resists free boosting** — the
preferred frame showing through in the *spatial* sector, at *lower* velocity than the temporal LIV
onset. The width it does show shrinks *faster* than `γ` (profile self-adjustment, not Lorentz
contraction), and the soliton's internal-clock readout is unreliable (it is not a clean
internal-clock eigenstate). **Length contraction + free boostability is not demonstrated.**

## So what (honest)

The make-or-break question does **not** cleanly pass. The universal clock is **half-hidden**:
invisible in **time** (dilation emerges), visible in **space** (PN pinning, and *earlier* than the
LIV scale). This is exactly the classic obstruction GPT flagged — *"discrete-substrate models tend
to run into Lorentz-invariance problems; a literal universal clock introduces a preferred frame"*
— now concrete: **Peierls–Nabarro pinning is the known mechanism by which a naive lattice breaks
Lorentz invariance in its spatial sector.**

This is a **real constraint, not a vindication**, and it is *productive*: it localizes the problem
precisely (the spatial sector, not the temporal one) and names the fix-direction. PN-pinning-free
("translationally invariant") discretizations exist in the literature — a viable Synchronism
substrate must use one; the naive nearest-neighbor lattice here is not it. That is a concrete,
falsifiable next requirement, not a vague gap.

**Unification (and a new, earlier signature):** the temporal LIV onset (`v_g ≈ 0.57`) joins
Phase-2 + Phase-4 at the grid scale. But the *spatial* frame-visibility (PN pinning) appears at
*lower* velocity — a **distinct, earlier preferred-frame signature** than the high-energy LIV
channels. If real, it would be more accessible than Planck-scale dispersion/Umklapp.

## Honesty

Reproducing SR time dilation is *reproduction* (stage-1, per [PREDICTIONS.md](../PREDICTIONS.md)
discipline 3), not novel — the discrete KG → continuum SR at low `k` is expected. The content
here is **diagnostic**: (a) the clock *does* hide in time (it could have failed), (b) the spatial
sector does *not* hide (PN pinning), and (c) the fix-direction is named. **Bucket 0 unchanged.**
The boosted-soliton internal-frequency measurement was unreliable and is reported as such, not
leaned on (agent-zero). This is the most important *constraint* the substrate arc has produced:
the discrete-Lorentz problem is real, located in the spatial sector, and has a known class of
fixes to test next.
