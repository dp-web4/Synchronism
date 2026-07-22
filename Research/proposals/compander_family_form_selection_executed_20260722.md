# Back-annotation: Compander-Family Form Selection EXECUTED on Real SPARC (2026-07-22)

**Source**: site explorer track, `synchronism-site/explorer/findings/compander-family-selection-executed-tanh-not-privileged.md`
**Script**: `synchronism-site/explorer/scripts/compander_family_aic_bic_real_sparc.py`
**Status**: EXECUTED (verdict rule pre-fixed before running; reproduces the prior γ=2 kill at ΔBIC=+184 as pipeline sanity check)

## What was run

Ten functional forms fit to the real SPARC RAR (2,807 points, Lelli+2016 mass models, McGaugh-2016 prescription, global fits, BIC = N·ln(SSR/N)+k·ln N): McGaugh ν, generalized ν-δ, tanh-log at γ=2 (the framework form) and free γ, erf-log, arctan-log, algebraic-sigmoid-log, Hill free-n and n=1 ("simple" μ), Gompertz. All minima confirmed global by independent shape-grid scan.

## Results

1. **tanh carries no statistical content.** Indistinguishable from erf-log/Hill/ν-δ (mutual ΔBIC ≤ 8.9 < 10 bar) — and finishes last among the four viable 2-param members. The C(ρ) functional form is confirmed arbitrary within its viable subclass, now by execution rather than assertion.
2. **The family is NOT degenerate.** SPARC refutes three qualitatively-acceptable companders on asymptotic-rate grounds: arctan-log +46.7, algebraic-log +23.8, Gompertz +58.0 (vs McGaugh). The data selects the asymptotic pair (deep-MOND μ∝y; Newtonian return at least power-law fast), and only that.
3. **γ→0.49 mechanism identified.** tanh(γln(1+y)) = 1 − 2(1+y)^(−2γ), so γ sets the Newtonian-return exponent q=2γ. Free fit: γ=0.489 ⇒ q=0.98; free-Hill independently pins n=0.975. Both parametrizations converge on q≈1 = MOND's simple μ (Hill n=1 ties McGaugh at k=1, ΔBIC −0.7). γ=0.49 is not a constant in need of derivation; it is the tanh-family encoding of the simple μ-function. γ=2 (framework) forces q=4 — premature re-Newtonianization; the +184 kill in asymptotic language.
4. **Cross-dataset squeeze (inferred, runnable).** The galaxy-preferred q≈1 return class is exactly what Cassini excludes in modified-gravity MOND via the EFE quadrupole (Hees, Famaey, Angus & Gentile 2016, MNRAS 455, 449 — ν̃ family rejected outright, "simple" α=1 excluded, only fast-return ν̄_α≥2 broadly viable; modified inertia exempt; Desmond 2024, MNRAS 530, 1781 sharpens the RAR-vs-Q₂ tension). tanh-log at γ=0.489 is curve-close to simple μ in the y~1–3 region Cassini probes, so it plausibly inherits the exclusion, while Cassini-safe large γ is SPARC-killed. **Runnable next step**: compute EFE Q₂ for tanh-log(γ) across γ, overlay Cassini bound on SPARC ΔBIC(γ). If no γ passes both, the tanh-log family is closed as modified gravity by a second route independent of the locality no-go.

## Whitepaper/site touchpoints

- Any statement that "several compander functions (logistic, erf, arctan, Hill, tanh) satisfy the constraints" is now partially wrong by execution (arctan refuted; logistic-in-log-argument ≡ tanh identically, so the list also double-counts).
- The governing-equation gap (form chosen, never selected) is CLOSED as a gap: selection has now been performed; verdict = form-selection null within the viable subclass, with real selection power at the family's edges.
