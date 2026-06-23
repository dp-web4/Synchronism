# Phase-10 (door #1, galactic) — do flat rotation curves fall out, and a Phase-8 correction (2026-06-23)

**Status:** `[ACTIVE-MRH]` — opens the galactic frontier (Phase-9's Bucket-0 door #1: the inflow
profile departs from GP via the capacity rule). **Result: it scopes door #1 to a precise,
falsifiable target — *and corrects a Phase-8 overclaim*: the GR-matching capacity rule does NOT
give flat rotation curves; it gives a slowly *rising* `r^{1/4}` curve. Flat curves require an
isothermal `ρ∝1/r²` halo (a *different* exponent), so a single rule must TRANSITION, with the
transition scale = the MOND `a0` / the RAR (bet B2).**
**Sim:** [`simulations/phase10_galactic_rotation_curves.py`](../simulations/phase10_galactic_rotation_curves.py) · result: `simulations/results/phase10_galactic_rotation_curves_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The relation and the result

Newtonian limit of the inflow/swimmer (Phase-9): a circular orbit has `v_c² = G·M(r)/r`. For a halo
`ρ ∝ r^(−α)`: `M(r) ∝ r^(3−α)`, so **`v_c ∝ r^(1−α/2)`**. Measured large-`r` slopes:

| profile | large-`r` `v_c` slope | behaviour |
|---|---|---|
| visible matter only (exponential disk) | −0.50 | **falling** (Keplerian — the dark-matter problem) |
| n=3 GR-capacity halo (`ρ∝r^(−3/2)`) | **+0.25** | **rising `r^{1/4}`** — *not flat* |
| isothermal halo (`ρ∝r^(−2)`) | −0.00 | **FLAT** ✅ |

## The correction to Phase-8 (agent-zero / honesty)

Phase-8 Faith-B said: *saturating capacity → cored centre + GR tail → "exactly what dark-matter
halos look like."* **That was an overclaim.** It conflated **cored** with **isothermal**:
- A dark-matter halo that produces *flat* curves is **isothermal**, `ρ ∝ 1/r²` (`α=2`).
- The GR-matching capacity rule's tail is `ρ ∝ r^(−3/2)` (`α=1.5`), which gives a **rising `r^{1/4}`**
  curve, not flat. A saturating centre doesn't change the tail.

So the n=3 (GR) capacity rule — and n=3+saturation — do **not** solve the rotation-curve problem.
The redirect was pointing at the right *door* but with the wrong *profile*.

## Door #1, now precisely scoped (the real requirement)

Flat curves need an **isothermal intent-halo `ρ∝1/r²`** at galactic scale — a *different* exponent
from the n=3 rule that matches GR near matter. So for one capacity rule to give **both** (GR/Keplerian
near a mass, flat curves at galactic scale) it must **transition**: `α: 3/2 → 2` as acceleration
drops. The transition acceleration **is** the MOND `a0` / the radial-acceleration-relation (the
repo's existing bet **B2**). Two requirements door #1 now makes explicit, *neither yet shown*:
1. **intent self-gravitates** — the flowing intent halo must itself source gravity, or only visible
   matter gravitates and curves stay Keplerian (no help);
2. **the capacity rule transitions** `α: 3/2 → 2` at a characteristic acceleration `a0`.

This is a sharp, falsifiable target — and the genuine loan-repayment direction, since it lives at
door #1 (profile-departure from GP), the *only* place outside discreteness where a novel prediction
can live (Phase-9).

## On "fit" (dp 2026-06-23)

`a0` would be a **fit** constant — and that's fine. Almost all physics constants are fit (`G`, `c`,
`α`, masses); "fit" means "this is the value that works, we don't yet know why," not "wrong." The
real test (productive vs degenerate) is whether **one** transition law with **one** scale `a0`
reproduces the RAR across **all** SPARC galaxies (few parameters, many data) — that would be good
physics even with `a0` fit. One knob per galaxy would be degenerate. That SPARC test is B2, the next
step.

## Honesty

Not novel physics yet; **Bucket 0 unchanged**. The value of this phase: (a) a real **correction** of
Phase-8's saturation→halo overclaim (GR tail rises, isn't flat); (b) door #1 scoped from a vague
"galactic frontier" to a precise falsifiable target (an `α: 3/2→2` capacity transition + self-
gravitating intent + `a0` = RAR); (c) the connection made explicit to the existing B2 bet. Caveats:
Newtonian-limit circular orbits, spherical-halo approximation, power-law profiles (the transition
itself not yet modeled — that *is* the open work). The next step is the actual B2/SPARC test: does a
single transitioning capacity rule reproduce the RAR.
