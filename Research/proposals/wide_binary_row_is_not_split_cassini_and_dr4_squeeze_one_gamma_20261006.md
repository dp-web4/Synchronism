# The wide-binary row is not "split": Cassini and Gaia DR4 squeeze the same γ

*Site maintainer, 2026-10-06 (WAKE). Source: visitor researcher persona, 2026-10-06 log, plus an order-of-magnitude
check run this session. No bucket moves. Count 6. Bucket 0 = 0.*

## The claim being corrected

Since 2026-09-27 the site (TEST-02 card, /tier-1-existing "inherited" panel) and the 2026-10-04 stopping table
(`synchronism-site/explorer/findings/the-stopping-table-is-empty-for-novelty-not-for-difference-...md`, row 6) say the
wide-binary test is **split**. On that reading a Chae-type boost refutes the density-keyed branch and "leaves the
acceleration branch standing with MOND". A Banik-type null does the reverse. So each outcome spares one realization.

## Why that is wrong on the ledger's own terms

1. **Acceleration branch.** TEST-25 already counts the acceleration-keyed C_g at the SPARC fit (γ ≈ 0.489) as failed,
   through QUMOND + Galactic EFE (g_ext ≈ 1.8 a₀). The wide-binary boost comes from the same physics: the same ν, the
   same field equation and the same g_ext. So a Chae-type boost cannot leave standing a member that is already failed.
   The γ that passes Cassini (γ ≳ 1.5–2, PREDICTIONS 2026-09-17 block) is a different member of the family.
2. **Density branch.** Under ambient keying, TEST-02's γ_g ≡ 1 identically, because the ratio cancels. SPARC (ΔBIC
   +2843) and LLR Ġ/G (68/74 laws at ℓ = 1 pc) already exclude this branch. A Chae boost would add a third root, not
   a first one.

## The quantitative version (estimate, not pre-registered)

`synchronism-site/maintainer/scripts/wb_boost_vs_gamma_efe_bracket.py` (+ `_output.txt`). The compander is
C(g) = tanh(γ ln(1 + g/(γA))), with A = a₀′/γ = 1.08×10⁻¹⁰ m/s², the identified combination. In the EFE-dominated
quasi-1D bracket the internal boost in g lies in [ν_e(1+L_e), ν_e], with g_e = 1.6–2.2×10⁻¹⁰ m/s²:

| γ | boost bracket at g_e = 1.9×10⁻¹⁰ | Cassini (TEST-25 instrument) |
|---|---|---|
| 0.489 | 1.16 – 1.58 | fails (+17.95σ) |
| 1.0 | 0.97 – 1.30 | fails |
| 1.5 | 0.92 – 1.22 | edge |
| 2.0 | 0.90 – 1.17 | passes (SPARC-disfavoured, 1.6–2.3σ galaxy-level) |

(The lower bracket end falling below 1 is an artefact of the crude 1D form, which takes the field-parallel component
alone. Only the ordering in γ is meant.)

So **within the one family**:
- a Chae-type boost of ≈1.4 pins γ ≲ 1, where Cassini already fails;
- a Newtonian null (γ_g ≈ 1.0 ± a few %) fails every γ ≲ 3;
- a precise intermediate boost of ≈1.1–1.2 is the only outcome that lands on the Cassini-passing γ. Even there it is
  the γ = 2 member that SPARC disfavours at ~2σ.

The DR4 row is therefore **not split**. It is a **three-way squeeze on one parameter**: SPARC shape, Cassini Q₂ and the
WB amplitude each pick a γ. It has exactly one winning branch, a narrow intermediate boost at γ ≈ 1.5–2. That branch is
still MOND-class (a different interpolating function), so it cannot move Bucket 0. It can only decide which compander
member, if any, survives all three.

## Proposed ledger changes (gate on dp)

1. Stopping-table row 6: "live but split" → "one-parameter squeeze (SPARC × Cassini × WB amplitude); the only surviving
   outcome is an intermediate boost at γ ≈ 1.5–2; cannot move Bucket 0."
2. TEST-02 status: drop "each outcome spares one branch". Replace it with: "the acceleration branch at the SPARC γ is
   already failed (TEST-25); a Chae boost adds no rescue; the density branch is identically null in the ratio and
   excluded by SPARC + LLR."
3. Explorer: replace the quasi-1D bracket with a proper QUMOND wide-binary computation (Chae's γ_g statistic, projected,
   with the Galactic EFE) across γ ∈ [0.3, 3]. Then pre-register the DR4 reading **before** DR4 arrives (planned
   Dec 2026), as a γ interval per outcome.

## What this does not do

- It does not add a refutation. Nothing new has been executed against data.
- It does not show the γ ≈ 1.5–2 window is empty. That needs the proper computation and the galaxy-level SPARC interval.
- The bracket is crude (1D, EFE-dominated, no projection, no Chae selection function). The γ ordering is robust, and
  the amplitudes are ±0.1-level.
