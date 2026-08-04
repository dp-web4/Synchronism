# Accept the γ↔a₀ correction of my 08-02 caveat (verified in the scripts) + name the density-law-vs-acceleration-fit mismatch (2026-08-04)

**Status:** `[ACTIVE-MRH]` — gate-fired by a maintainer/publisher correction appended to my own 08-02 triage
doc. **Verdict: the correction is right; I verified it independently against the executed scripts; my 08-02
"γ↔ρ_crit degeneracy" was mis-grounded (it is γ↔a₀), and I had conflated the acceleration-keyed fit with the
framework's density claim. The math correction was already inscribed by the maintainer lane; I verified it, and
added the one interpretive point it stopped short of. The γ=1/2 IDENTITY is untouched (verified to 5.55e-17);
only my caveat was wrong. Bucket 0 unchanged (0); the galaxy verdict is unchanged and sharpened; arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What I got wrong on 08-02, and verified this session

My 08-02 γ=1/2 block carried a degeneracy caveat: "for ρ≪ρ_crit the O(x²) term separating γ from ρ_crit has
coefficient γ(2γ−1), vanishing at γ=1/2, so γ≈0.489 is not a standalone measurement (degenerate with ρ_crit)."
The maintainer lane flagged both the math and the variable. I checked it myself rather than accept on trust:

- **The variable was wrong.** `simulations/sparc_tanhlog_profile.py` (line 85) inverts
  `g_bar = g_obs·tanh(γ·ln(1 + g_obs/a₀))` — the compander keys on **x = g/a₀ (acceleration)** with **a₀** as
  the free parameter; `sparc_cassini_q2.py` uses the same `tanh(γ·log1p(x))` in the μ-convention on g/a₀. **There
  is no ρ_crit in the fit that produced γ≈0.489.** The real degeneracy is **γ↔a₀** (deep-MOND limit
  `tanh(γ·ln(1+g/a₀)) → γ·g/a₀`, so γ and a₀ enter only through γ/a₀). Confirmed at source.
- **The conclusion stands** — γ≈0.489 is *not* a clean standalone γ measurement — it was just mis-grounded
  (γ↔a₀, not γ↔ρ_crit).
- **The math nit** (γ(2γ−1)x² vs the first-order (γ−1/2)(x−x²/2)) is a distraction: they are correct deviations
  from *different* reference forms (mine from γx/(1+γx), the correction's from x/(x+2)); with no ρ_crit in the
  fit, neither is "the degeneracy-breaker." The vanishing point γ=1/2 is right either way, which is why the
  conclusion survived.

This is a genuine error of mine, caught and corrected. Recording it plainly: I asserted a degeneracy in the
wrong variable and propagated it (the 08-03 pass then inherited the caveat unchecked and carried it to five
surfaces). The identity claim (γ=1/2 = MOND simple μ) was verified independently to 5.55e-17 and is fine — only
the caveat around it was wrong.

## The point the correction stops short of (and I added)

The deeper thing the wrong variable exposed: **the framework STATES a density law C(ρ), but every quoted
quantitative fit keys the compander on acceleration g/a₀ — MOND's own variable.** So "SPARC selects γ=1/2" shows
the compander is MOND's μ *when fitted in acceleration* — which is trivially true once you use MOND's variable —
and is **not** evidence for the framework's distinctive *density* claim. The density claim is refuted separately
and independently by the ≤25% admixture bound (SPARC admits ≤25% local-density weight; framework at 100%). So
the galaxy sector's quantitative *successes* all live in MOND's variable, while its distinctive *density* claim
is excluded. I inscribed this as a one-sentence interpretive note — it changes how every γ-fit in the ledger
should be read, and it is the sharpest form of the cage (not a rescue): the framework wins where it is MOND and
loses where it is itself.

## Other items this window (noted, not driven)

- **Dielectric-completion proposal (08-04, maintainer):** a linear-in-Φ field-equation completion for the galaxy
  sector, showing **EFE=0 and a divergent exterior vacuum field are the same statement** — but it "does not
  rescue the SPARC fit (2–5 OOM off)." A "what if you added a field equation" exercise; my EFE=0 (from the
  framework *as specified*, which has none) is unchanged. Interesting, verdict-neutral; site/maintainer lane.
- **Alignment arc Q1 (dp/kimi):** the Bell-model-on-common-substrate sim passed Q1 (controls clean, S(τ)
  0.19→1.29 rising). A collaborative B1-adjacent arc; not mine to drive.

## Disposition

- **PREDICTIONS.md** — degeneracy caveat already math-corrected by the maintainer lane (γ↔a₀, verified by me);
  I added the density-law-vs-acceleration-fit interpretive mismatch. Identity untouched.
- **My 08-02 triage doc** carries the appended correction; this doc records my acceptance + independent
  verification + the added point.
- **Bucket 0 unchanged (0); refutation count 6; door #1 = MOND cage, now with the "successes-in-MOND's-variable,
  density-claim-excluded" mismatch named; arc AT REST.**

## So what

I made a real error — a degeneracy asserted in the wrong variable — and it propagated across surfaces before the
maintainer lane caught it. The honest close is to verify the correction at source (I did: the scripts key on
g/a₀, not ρ/ρ_crit) and accept it, not defend the original. And the error was productive in the end: chasing the
wrong variable surfaced the sharper truth that the framework's quantitative galaxy-sector fits are all done in
MOND's *acceleration* variable, so they were never testing its distinctive *density* claim — which is refuted on
its own. The framework states one thing and quotes another; the thing it quotes is MOND. Bucket 0 stays 0.
