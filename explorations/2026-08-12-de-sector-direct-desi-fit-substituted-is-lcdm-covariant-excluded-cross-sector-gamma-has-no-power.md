# DE sector, direct DESI DR2 fit: the substituted branch IS ΛCDM (measured), both covariant completions fail the fit, and the "cross-sector γ concordance" has no power to fail (2026-08-12)

**Status:** `[ACTIVE-MRH]` — gate-fired by a new *executed* result in my own DE-sector lane (not door-#1 galaxy, so the park rule does not apply; a new likelihood fit clears it). **Verdict: triaged and inscribed deflationarily. The covariant algebra I re-ran in-repo myself (PASS); the likelihood Δχ² numbers are site-executed (caveated, not re-run by me). Three findings, all deflationary: (a) the substituted branch is ΛCDM at γ=0.487, measured — the 08-11 "forced-wₐ 3.4–6.3σ" exclusion pricing DIED on execution; (b) both covariant completions FAIL the fit (A=EdS χ²≈9900, B Δχ²≥+79) — stronger than the prior quadrant-membership framing; (c) the first cross-sector γ test (γ_cosmo=0.487 vs γ_galaxy=0.489, "0.1σ") has NO POWER to fail and must not be promoted to a concordance. Bucket 0 unchanged (0); count unchanged (6); still Bucket 3; TEST-26 statistic change gates on dp.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What fired

Three 08-12 commits landed on my last HEAD (`ad8cb889`, the 08-11 both-directions guard on the DE-sector retraction), all in the DE/DESI/TEST-26 line:
- `f6de09d3` Publisher — appended to my DE block: the site-side covariant execution is re-run in-repo (`simulations/publisher_20260812_covariant_00_checks.py`), the "not derived from a covariant action" caveat discharged, the no-go hardened to class level.
- `0cc468cb` site maintainer — TEST-26 registration draft (restate class-level + projection-robust); gates on dp.
- `dabcdb39` (HEAD) site explorer — **a direct likelihood fit of the w(z;γ) family to DESI DR2 + Planck-compressed CMB + DES-SN5YR**. This is the new *executed* object.

Per the standing park discipline: door-#1 galaxy gate-fires get verify/route, not new framings — but this is the **DE/cosmology sector** and a **new executed likelihood fit**, which is exactly the "genuinely new executed result still gets triaged" exception.

## What I verified myself (in-repo)

Ran `simulations/publisher_20260812_covariant_00_checks.py` — ALL PASS:
- **Class identity** `w_DE = d(ln F)/d(ln x)` for `ρ_DE = ρ_m·F(x)`, symbolic — any F.
- **Completion A** (Appendix D §D.3 as written): Bianchi forces `ρ/C ~ a⁻³` ⇒ `H² ~ a⁻³` = **exact Einstein–de Sitter**, `q = +½`, no acceleration; vacuum floor `ρ_crit/γ`; universe **ends at a_end = 1.0372** under Session100's calibration (no FRW solution past it).
- **Completion B** (Brans-Dicke-pinned closure): `H² = 8πG ρ_m/(3 C B)`, `B = 1 − 3ε − (3ω/2)ε²`, `ε → 1` as x→0 ⇒ `B → −2` (the w=−1 attractor is destroyed). Table anchors reproduced for γ ∈ {0.2, 0.489, 0.5, 2.0}, ω=0: **every w₀ > −1 member has wₐ > 0** ⇒ no member reaches the DESI quadrant (which needs w₀>−1 AND wₐ<0).

So my 08-11 caveat "(ii) covariant route executed site-side, *not re-run by me*" is now **discharged** — twice-independently verified (Publisher + me).

## What is site-executed (caveated, NOT re-run by me)

The likelihood fit itself (`synchronism-site/explorer/findings/scripts/fit_gamma_family_to_desi_dr2.py` — site is reference, I read it, did not re-run the chains). Its numbers, taken on the pipeline's own validation (ΛCDM BAO χ²=12.6/13; w₀wₐ recovers DESI's (−0.42,−1.75); full-combo crossing Δχ²=−11.0 ≈ 2.9σ, inside the published 2.8–4.2σ band):

| Model | Δχ² vs ΛCDM | reading |
|---|---|---|
| w₀wₐCDM (+2) | −11.0 | the crossing this data prefers |
| **substituted (+1)** | **−0.3** | best γ=0.487; nests ΛCDM at γ=½; +11 *behind* w₀wₐ |
| completion B (ω∈{0,1,5,50}) | **+79 … +187** | fails the fit, every ω |
| completion A (EdS) | **χ²≈9,900** | the pre-1998 background |

## The three deflationary points

1. **The substituted branch IS ΛCDM — measured, not asserted.** My 08-11 both-directions guard said the DE sector is a tautological reparametrization (Bucket 3). The fit confirms it at likelihood level: the family sits at its Λ-corner (best member ≈ CPL (−0.993, +0.023)) and pays exactly ΛCDM's price. Critically, the **08-11 "forced-wₐ 3.4–6.3σ" exclusion pricing of the substituted branch died on execution** — the likelihood *never forces* the family to DESI's w₀; there was no exclusion, only non-novelty. (That number was site-side, not on my ledger; this is an external repricing I verify, not a self-correction. It is the 5th "no-go stands, exclusion-flavored number dies" instance — a familiar shape.)

2. **Both covariant completions fail the fit** — a strictly *stronger* statement than the "0/192 γ miss the quadrant / forced wₐ wrong-signed" framing already in the ledger (which is about quadrant *membership*). A misses by χ²≈9900 (it's EdS); B by Δχ²≥+79 for every ω. So the DE sector survives **only** non-covariantly and **only** by being ΛCDM. This is the same pattern the whole arc keeps hitting: the algebraic substitution reproduces a known result; the covariant route that would make it a *prediction* is excluded.

3. **The cross-sector γ test has NO POWER to fail — do not promote it.** This is the both-directions catch of the session. γ_cosmo = 0.487±0.02 (DESI) vs γ_galaxy = 0.489 (SPARC) "agree at 0.1σ" reads like Synchronism's first cross-sector confirmation — two independent sectors landing on the same shape parameter. It is not. **γ=½ is exactly Λ** (the Möbius/tautology point) and **γ=0.489 is exactly MOND simple-μ**; the two sectors' *standard models* sit 0.011 apart in γ-space by construction, so the test could not have failed. Separating 0.489 from 0.500 needs σ_γ ≈ 0.004, and the SPARC-side σ(γ) has never been derived. The agreement is inherited from Λ+MOND, not evidence for Synchronism. Left un-guarded in a proposal, a "0.1σ cross-sector concordance" is exactly the kind of number that gets harvested into the site as a success — so the deflationary reading is now the *ledger* statement.

## Disposition

- **PREDICTIONS.md** — appended to the DE block (keyed to the HEAD fit): my in-repo re-run discharges the "(ii) not re-run" caveat; the substituted branch measured ≡ ΛCDM; both covariant completions fail the fit; and the cross-sector γ concordance registered **deflationarily (no power to fail)**. Count unchanged (6); Bucket 3 unchanged; no novel DE prediction.
- **Gates on dp (routed, not inscribed):** the TEST-26 registration statistic change (quadrant-membership → Δχ²(substituted vs w₀wₐCDM) on DR3, as a *consistency check* not a discriminator) — a governance act, same class as the DR2 three-branch pre-commitment. The count recount also gates on dp (held at 6).
- **Flagged, not mine:** the w_eff sign error (Sessions 100/101, archive lane, unchanged); the SPARC-side σ(γ) derivation (owed before any concordance could ever be priced — but the point is it *can't* discriminate even then).
- **Bucket 0 unchanged (0); count 6; arc AT REST.**

## So what

The DE sector got its first real likelihood contact with data, and it behaved exactly as the reparametrization diagnosis predicted: the branch that *lives* (the algebraic substitution) is ΛCDM to Δχ²=−0.3 and can't touch the crossing DESI prefers; the branches that would *differ* (the covariant completions) are excluded by the fit — one of them is literally pre-1998 EdS. The tempting headline — "γ_cosmo = γ_galaxy at 0.1σ, cross-sector concordance!" — is a no-power test, because both sectors reduce to standard models (Λ, MOND) that already sit at those γ values. Same lesson as every door-#1 session, now on the cosmology axis: Synchronism succeeds precisely where it reproduces a known model, and is excluded precisely where it would say something new. The value I added was not another framing but the deflation: computing that the concordance had no power to fail, and re-running the covariant algebra myself so the caveat is discharged honestly rather than on trust. Bucket 0 stays 0.
