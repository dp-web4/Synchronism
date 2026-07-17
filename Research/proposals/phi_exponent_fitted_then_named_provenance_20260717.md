# φ Provenance Audit: Fitted-Then-Named — S45's Own Guardrail Was Breached Silently

**Date**: 2026-07-17
**Source**: synchronism-site explorer track
(`explorer/findings/2026-07-17-phi-exponent-provenance-fitted-then-named.md` — full chain
with per-file line citations)
**Status**: Provenance audit — proposes corrections, gates on dp for whitepaper edits

## Result

The golden ratio has **zero surviving derivations** in this archive and **three mutually
independent retro-justifications** that do not cite each other:

1. **S186** "information conservation x + x² = 1" — circular: x + x² = 1 *is* the
   defining identity of 1/φ; the "quadratic delay" premise appears nowhere else in the
   archive and has no independent motivation.
2. **S219** "scale recursion THEOREM" — λ = 1 + 1/λ is inserted (nothing in the
   postulated recursion forces it) and no link between recursion factor λ and exponent β
   is computed; the Fibonacci "verification" table verifies the definition of φ.
3. **Whitepaper** (line ~4123) "optimal balance between local and non-local coherence" —
   asserted, uncomputed.

**S218 itself concedes the point**: its own Boltzmann route "gives exponent 1, not 1/φ"
(line 98), and it poses "Why 1/φ?" as an open question eight days after S186 claimed to
have closed it.

## The archive already ruled on this — and the ruling was never cited again

**S45 (2025-11-25)**, adjudicating S44's B = 1.62 ≈ φ: "0/8 other astrophysical scalings
have φ exponents … B = φ vs B = 1.62: same success rate … INTRIGUING COINCIDENCE but NOT
SIGNIFICANT … **Don't: Claim φ is fundamental**."

No session from #46–#169 mentions φ. It reappears fully formed inside C(ρ) at **S170**
with no derivation and no citation. This is the second instance of the pattern already
documented for the BTFR (S58 honest record → S193 rescue): **early honest adjudication →
un-cited override → override becomes load-bearing.**

## The empirical ledger points at other constants

| Slot | Measured | Nearest constant | φ-candidate distance | Chosen |
|---|---|---|---|---|
| S185 C(ρ) exponent | 0.66 | 2/3 (1.0%) | 1/φ: 6.8% | 1/φ |
| S217 a₀ exponent | 1.469 needed | 3/2 (2.1%) | φ: 10.1% | φ |
| S239 Gaia transition | 0.688 ± 0.10 | 2/3 (3.1%) | 1/φ: 11.3% (0.009 inside 1σ edge) | "VALIDATED" |
| Gnosis salience | 0.40 | — | 1−1/φ: 4.7% | "intriguing" |

S192 states the actual selection criterion in writing: the φ formula for a₀ is preferred
for "**symmetric structure** (coherence uses 1/φ inside, a₀ uses φ outside)" — brand
consistency over fit. S217's own numbers show Ω_m^φ ≈ 1/(2π) to 3.1%: the a₀ formula
works exactly to the extent it re-encodes Milgrom's dimensional a₀ ≈ cH₀/(2π)
(exact-match exponent: 1.591, also not φ).

## Proposed corrections (gate on dp)

1. **Whitepaper line ~4049**: "Golden ratio exponent 1/φ | VALIDATED | 1σ" — an
   edge-of-interval consistency (S239: 1σ = [0.609, 0.802], 1/φ = 0.618, 2/3 = 0.667
   near center) cannot carry VALIDATED. Downgrade to "consistent within a weak
   constraint; 2/3 fits 3× closer."
2. **PREDICTIONS.md** (if/where φ is described as derived): mark φ fitted-then-named,
   citing S45's original adjudication as the correct precedent.
3. **Proposed audit rule** (extends the citation-walk): for any parameter labeled
   "derived," walk the derivation chain; **two independent derivations that do not cite
   each other ⇒ treat as fitted-then-named** until one derivation survives scrutiny.
   Real derivations accumulate; retro-justifications get re-invented. (Ω_m's floor has
   one derivation stated identically everywhere — that's what derived looks like.)

## Consequence for TEST-09

Strengthens it. The BTFR kill was framed as "the two derived-from-cosmology ingredients
(Ω_m, φ) put it off the BTFR." Correct form: one derived ingredient (Ω_m) and one free
exponent dressed as a constant — and the 2026-07-14 no-rescue scan already showed the
kill fires for **every** exponent value at the framework's own Ω_m. A free parameter that
still cannot reach the data is a deeper failure than a derived one that cannot.
