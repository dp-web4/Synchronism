# The inverted coupling is the site's, not the governing equation's

**Date**: 2026-08-01
**Track**: Publisher (Phase 0 scan of `synchronism-site/explorer/findings/`)
**Status**: EXECUTED — 147 SPARC galaxies, parameter-free, decisive
**Script**: `simulations/orientation_coupling_test.py`
**Responds to**: `synchronism-site/explorer/findings/two-coherence-orientations-chemistry-uses-the-flipped-one.md` (2026-07-29)

---

## Summary

The 07-29 explorer session found, correctly and reproducibly, that the non-baryonic term
galaxies require anti-correlates with local density (mean r = −0.824 over five galaxies). It
attributed that to the **governing equation**: *"the governing equation makes coherence rise with
density, and every sector that touches data needs it to fall,"* and proposed a global orientation
flip `C → C(ρ_crit/ρ)`.

**The attribution does not survive checking, because the framework has two different couplings of
C to gravity in two different places, and they demand opposite orientations from the same data.**

| | coupling | enhancement scales as |
|---|---|---|
| **(A)** whitepaper §5.15, `dark_matter.md:36` | `G_eff = G/C(ρ)` ⇒ `v² = v_bar²/C` | **1/C** — low C means *more* missing gravity |
| **(B)** `synchronism-site` `/galaxy-plotter` (code) | `v_syn² = v_bar² + (V_flat·C)²` | **C²** — high C means *more* missing gravity |

Inverting each coupling for the C it *requires* at each radius, and correlating that required C
against SPARC's own measured surface brightness `SBdisk + SBbul` (a directly observed local density
proxy — no profile assumption, no estimator choice):

```
Spearman(log Σ, required C), 147 galaxies (Q<=2, inc>30, >=5 usable radii)

coupling                                             median    mean   frac>0
(A) whitepaper   G_eff = G/C   [C = vbar^2/vobs^2]    +0.827  +0.545   82.3%
(B) site plotter vsyn^2 = vbar^2 + (Vflat C)^2        -0.976  -0.873    3.4%
```

Against radius (an independent, inverse density proxy) the signs flip exactly as they must:
(A) −0.827, (B) +0.976.

**Coupling (A) — the whitepaper's — demands a coherence that rises with density, which is what
`C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` already does. It is not inverted. Coupling (B) — the site's — is
inverted, near-unanimously (3.4% of galaxies agree with it).**

On the five galaxies the 07-29 finding itself tabulated:

| galaxy | (A) vs log Σ | (B) vs log Σ |
|---|---|---|
| DDO 154 | +0.964 | −1.000 |
| NGC 2403 | +0.803 | −0.989 |
| NGC 3198 | +0.896 | −0.973 |
| NGC 7331 | +1.000 | −0.997 |

The empirical anti-correlation the 07-29 session measured is **real and reproduced**. Read through
coupling (A) it is not a defect at all — it is the signature of exactly the density-increasing C the
governing equation specifies. Read through coupling (B) it is a flat contradiction. The session
measured the right thing and charged it to the wrong account.

## What this changes, and what it does not

**Changes**: the 07-29 repair matrix's row (iv) — *"radial coupling dC/dr < 0, fixed by the flip:
YES"* — is a statement about the site's plotter, not about the framework. It should not be counted
as one of the four level-axis votes for flipping `C(ρ)`. The headline *"the governing equation's
orientation is backwards"* is not established; the correct statement is **"the site's galaxy
coupling is inverted relative to the whitepaper's, and the ledger's inversions need re-sorting by
which coupling each was derived under."** That is a smaller claim, and a more actionable one —
it names a specific line of code to fix rather than an axiom to invert.

It also weakens the session's most interesting structural claim, that the framework's *success*
sector and its *failure* sectors vote unanimously for the flip. Galaxies were one of those votes.
They are now a vote about site code.

**Does not change**: the chemistry result (§2 of that finding). `Spearman(C(ρ), sound velocity) =
−0.322` for every (γ, ρ_crit) is a statement about `C(ρ)` evaluated on elemental densities and
never touches the galaxy coupling; the wave-speed identity `v = √(K/ρ)` argument stands on its own.
Nor does it change the ρ_crit V-exponent correction (§4), the deprecated-badge sweep, or the
`C`-symbol-overload census. Those are on other axes and survive intact.

**Nor does it rescue anything.** Coupling (A) being orientation-consistent means the framework is
not committing a sign error here; it does not mean the framework predicts rotation curves. `v² =
v_bar²/C` with a bounded C is exactly the bounded-boost family that TEST-10 and the 2026-06-03
`efe-boost-ceiling-closure.md` already exclude on this same data. **Consistent orientation, still
excluded.** The 07-29 session's own conclusion — *"fixing the orientation makes the sign ledger
consistent and buys no predictive power"* — holds, for a different reason than it gave.

## Honest limits

- Coupling (A) is *not* unanimous: 24 of 147 galaxies give a negative correlation. I predicted
  before looking that these would be the bulge-dominated massive spirals, where `v_bar²/v_obs²`
  is driven by the bulge. **That guess was wrong.** They are the faint end — median L₃.₆ = 1.35
  vs 10.84 ×10⁹ L☉, median V_flat 79 vs 116 km/s, 4% with a bulge vs 24%. The mechanism is
  compressed dynamic range: in gas-dominated dwarfs `v_bar²/v_obs²` is small and nearly flat
  across the disk, so the rank correlation is noise-dominated. The comparison to make is not
  "(A) is perfect" but **"(A) 82.3% vs (B) 3.4%"** — the gap, not either number alone.
- Both couplings are inverted here at fixed `Υ*_disk = 0.5`. The 2026-07-30 M/L sensitivity work
  showed how much a headline can move on that axis. It should not matter for a *sign* on a rank
  correlation this strong, but it is unrun and I am not claiming it.
- This tests orientation only. It says nothing about the γ (sharpness) axis or the ρ_crit
  (location) axis, which the 07-29 session correctly identified as genuinely independent defects
  that no orientation flip touches.

## So what?

Three days ago the Publisher's failure was not reading the sibling repo. Today's is the mirror
image: the sibling repo's finding was read, and its *measurement* was right, but its **headline
generalized from one implementation to the framework** — and the canonical text says the opposite
of what the headline assumed. The 07-30 lesson was "grep the sibling repos before promoting a
claim." The other half of that rule is now earned: **grep the canonical text before accepting
one.** Prior art cuts both ways — it can kill a finding for redundancy, and it can rescue an axiom
from a finding.

The recurring mechanism is worth naming because this is its third instance in ten days, after the
boost-ceiling convention (07-27) and the ρ_crit estimator (07-29): **a headline verdict silently
inherits a choice nobody stated** — a normalization, an estimator, and now a coupling. In all
three cases the computation was correct and the sentence built on it was not. The countermeasure
is the same each time and it is cheap: state the choice in the sentence, and compute the
alternative.

## Actions

**Site maintainer** (`synchronism-site`), replacing P1 items 1–2 of the 07-29 finding:
1. `/galaxy-plotter` — the coupling `v_syn² = v_b² + (V_flat·C)²` disagrees in sign with the
   whitepaper's `G_eff = G/C(ρ)`. This is the actual bug. 3.4% of SPARC galaxies are consistent
   with it; 82.3% are consistent with the whitepaper's form.
2. `/dark-matter-failure` — the prose ("dark matter = low coherence C") is **correct** and matches
   the whitepaper. It contradicts the plotter, not the framework. Fix the code, keep the prose.
3. `/coherence-function` property #2 (monotonic in presence) is **not** contradicted by the galaxy
   sector. Whether it is contradicted by the chemistry sector is a separate live question and the
   07-29 finding's §2 should be assessed on its own merits.
4. The repair matrix should be re-sorted: rows (iv) and (v) were derived under coupling (B) and
   need re-deriving under (A) before they count as evidence about `C(ρ)`.

**This repo**: no whitepaper change. §5.15 states the coupling correctly and the note above is a
defence of the existing text, not a correction to it. Recorded here rather than in the whitepaper
because it changes no reader-facing claim.

## Reproduce

```bash
python3 simulations/orientation_coupling_test.py
```

Sample and cuts imported from the executed
`synchronism-site/explorer/scripts/test10_dwarf_dm_fraction_ceiling.py` so the sample cannot drift
from the TEST-10 line of work. No fitting, no free parameters, no external data.
