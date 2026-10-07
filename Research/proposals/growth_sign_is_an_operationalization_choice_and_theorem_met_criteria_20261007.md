# Proposal: The growth sign is an operationalization choice; two Bucket-1 criteria are met by theorem

**Filed**: 2026-10-07 (site maintainer, WAKE phase)
**Source**: synchronism-site visitor log 2026-10-07. The researcher persona raised the growth-sign question; the
graduate-physics and researcher personas raised B1.
**Bucket moves**: none. Count 6; Bucket 0 = 0.
**Gates on dp**: items 1(c) and 2.

---

## 1. Three readings of "G/C" in the growth equation, and they disagree on the sign

The site carries two growth numbers for one observable. /dark-energy gives an fσ₈ shift of −0.22 %, and TEST-04a gives
fσ₈(0.51) = 0.418 (−12 %). A researcher persona asked: if gravity is G/C with C ≤ 1, why is growth ever *suppressed*?

The archive contains three ways to put C into the growth source. Nobody had written the third one down.

| Reading | Growth source | fσ₈(0.51) | Status |
|---|---|---|---|
| (S) Session 107 | G_local/G_global = C_cosmic/C_galactic < 1 | 0.418 | TEST-04a's number. Rests on C_cosmic ≠ C_galactic, withdrawn 2026-08-11 (08-18 block). |
| (F) DE fluid | δ_DE/δ_m = 1 + w_DE, gravity unmodified | ≈ 0.473 | ΛCDM-like. The local fluid's Jeans term pins γ = ½ to 10⁻⁵ (09-15 block). |
| (U) Uniform G/C | The same G → G/C that gives H² = 8πGρ_m/(3C) also sources δ. With C ≡ Ω_m(a) (09-14 identity), μ = 1/Ω_m(a), for every γ | **0.575 (+21 %)**, σ₈(0) = 0.916 | **New, executed 2026-10-07.** |

Executed: `synchronism-site/maintainer/scripts/growth_under_uniform_G_over_C.py` (+ `_output.txt`).
- Setup: ΛCDM background (Ω_m0 = 0.315; the substituted family sits at Λ's corner), D normalized at a = 10⁻³.
- Identity control: μ = 1 reproduces f(0.51) to 3×10⁻⁴ of Ω_m^0.55, and fσ₈(0.51) = 0.4742.
- Predictions were written in the script header before it ran. P1 (D₀ ratio > 1.5) **failed** at 1.13. P2 (fσ₈ > 0.8)
  **failed** at 0.575. P3 (identity) held. The enhancement is late-time-limited: μ = 1.6 at z = 0.51 and 3.2 today, but
  growth has mostly stopped by then. I expected a blow-out and got +21 %.

What this changes:

(a) **The sign is a modelling choice.** The framework's own G/C, applied uniformly, *enhances* growth. The 05-09/05-12
"suppressor class" diagnosis asked whether C_galactic/C_cosmic was inverted. The cleaner statement is that suppression
needs a ratio of two C's, and the framework no longer has that ratio.

(b) **(U) is not a win, and it is not prospective.** It sits +0.4σ from DESI DR1 LRG1 (0.550 ± 0.062), but that bin
was published before this reading was written down, so it is a retrodiction. Its σ₈(0) = 0.916 is 2.2σ above DESI's
0.841 ± 0.034, which is a GR-conditioned number. If lensing sees the same μ, S₈ ≈ 0.94 is far above DES Y3 and
KiDS-Legacy. Cardassian modified-gravity branches are also ISW-excessive (Koivisto+2005, imported). Its likely fate is
Bucket 2 on existing data. That is **unexecuted**: it needs a Σ (light-deflection) specification, which the archive
does not have. Explorer topic seeded: `uniform-g-over-c-growth-vs-s8-isw-and-mu0.md`.

Structurally, (U) is worth noting for the nesting table (09-27). It is parameter-free, it differs from ΛCDM, and it is
not nested in it, which is the shape a Bucket-0 candidate would need. It arrived post-hoc, though, and points toward a
probable refutation.

(c) **For dp (TEST-04a DR2 registration).** As adopted, the DR2 pre-commitment tests reading (S), a mechanism the ledger
already calls withdrawn. Reading (U) predicts the opposite direction. Recommendation: before DR2 full-shape publishes,
the registration names its reading, or it is retired as "tests a withdrawn mechanism". Without that, a DR2 fσ₈ near 0.55
would fire branch B against (S) and *agree* with (U), and the record would carry both readings unlabelled. That is the
same unrecorded-as-different pattern as Session 100 vs Session 107 (173×). This joins the 09-21 recommendation that
branch B is underpowered.

## 2. Two Bucket-1 refutation criteria are met by theorem (B1, B6)

- **B1:** "Refuted if observer-relative statistics obey CHSH S ≤ 2." For the local arm this is Bell's theorem. The
  two nonlocal-grid arms are construction nulls; they are not theorems (the site's 09-14 relabel is correct).
- **B6:** clause (a) is met by no-signaling alone (10-06 block).

Two graduate/researcher personas again (2026-10-07) read refutation #6 as "Bell's theorem on the scoreboard". The
site already splits it as "5 on external data + 1 construction check". **Recommendation (gates on dp; first routed
10-04 as "count texture"):** adopt a ledger rule that a refutation criterion satisfiable by theorem for any
construction in the class is a *construction check*, not a refutation. Then quote the headline as
"5 refutations on data + 1 construction check". The underlying facts do not change; the count's label does.

## 3. Record notes (no verdict change)

- **The V⁺² exclusion (~11σ) is already galaxy-level.** One point per galaxy, N = 129, with a 4000-draw galaxy
  bootstrap (`explorer/findings/scripts/rho_crit_exponent_is_freemans_law.py`, line 47). It is
  **estimator-conditional**. Forward OLS of Σ_c on V gives V^(−0.15 ± 0.18); orthogonal regression gives V^(+2.1).
  The 08-27 finding argues forward is the question the law poses, and the gap is intrinsic scatter (r = 0.64). The site
  now states the unit and the estimator.
- **No Reparametrization badge on the site rests on an A2ACW verdict alone** (audit 2026-10-07, all badge sites). Each
  cites an identity, a computed null, a data fit, or named prior art. Thinly supported: the archive's FΣIR row ("standard
  M/L analysis", no specific result cited) and the Intent-field row (no identity shown); on the site, two gamma-boundary
  rows, σ_int, and the Core Idea regime badges. Caveat: the prior-art *matches* were made by LLM sessions.
- **A₀ ↔ ρ_Λ.** The researcher persona's "last cross-sector bet" is a derived κ in a₀ = κ c √(Gρ_Λ), fixed by the
  function rather than by the cH₀/2π coincidence. Seeded as an explorer topic. Its first obstacle is the 09-24 (ii)
  qualifier: the galaxy γ is C_g on acceleration and the cosmology γ is C_ρ on mean density, so one knee cannot carry
  both without a bridge between the variables.

## So what

Item 1 converts a recurring visitor confusion (−0.22 % vs −12 %) into a three-row table with one new executed number.
It also exposes a live governance hazard in the program's one adopted prospective registration. Item 2 is the third
time the count label has been raised, and the fix is a rule, not a recount.
