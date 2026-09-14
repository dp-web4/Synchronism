# The Dark-Energy Sector Is a Cardassian Model — a 24-Year Prior-Art Gap, and What It Pre-Answers

**Filed**: 2026-09-14 (site maintainer track, WAKE phase — before any site fix)
**Status**: proposal. Prior-art correction plus a literature-imported constraint. **No bucket moves. Count stays 6. Bucket 0 stays 0.**
**Raised by**: visitor researcher persona, `synchronism-site/visitor/logs/2026-09-14.md` (Pass 4, /dark-energy row)
**Bears on**: PREDICTIONS.md Bucket 3 (DE-sector note, 2026-08-10 → 08-18 block), TEST-26 registration (gated on dp),
`Session100_Modified_Friedmann.md`, `Session107_DESI_Forecasts.md`, the A2ACW H1/H2 question
**Script**: `synchronism-site/maintainer/scripts/de_sector_is_cardassian.py` (+ `_output.txt`)

---

## The gap (existence claim, verified at the primary layer)

`grep -ri "cardassian\|freese"` over the **whole** Synchronism repo (Research/, manuscripts/, explorations/,
simulations/, compilation docs), and over the whole site repo: **0 hits.** Per the standing 08-10 rule, this is a
negative-existence claim checked at the primary layer, not the compilation layer.

Session 100 (2025-12-08) writes `H² = 8πGρ_m/(3C)`. That is a modified Friedmann equation `H² = (8πG/3)·g(ρ_m)`
with no new degree of freedom — the defining form of **Cardassian expansion** (Freese & Lewis 2002, Phys. Lett. B
540, 1). From August 2026 on, the archive and the site have executed covariant completions, a DESI likelihood fit,
a locality fork and perturbation channels for this sector, and cited none of the literature on this model class.

## What the identification is, exactly (script, all checks at machine precision)

With `C = tanh(γ ln(1+x))`, `x = ρ_m/ρ_crit`, `1/C = 1 + 2/((1+x)^{2γ} − 1)`:

1. **Identity.** `ρ_m/(ρ_m+ρ_DE) ≡ C` (max |diff| 1.1×10⁻¹⁶ over γ ∈ {0.3, 0.487, 0.5, 2}). Cosmological C *is*
   the model's Ω_m(a). (Same content as the explorer's 09-11 floor identity; restated because the visitor asked
   for it on /dark-energy, where it was absent.)
2. **γ = ½ is exactly a Cardassian member:** `1/C = 1 + 2ρ_crit/ρ_m` (|diff| 2×10⁻¹³), i.e. ΛCDM (n = 0).
3. **High density (x ≫ 1):** `1/C → 1 + 2x^{−2γ}` — the modified-polytropic Cardassian (Gondolo & Freese
   2002–03) with **q = 1, n = 1 − 2γ, ρ_car = 2^{1/2γ} ρ_crit**; w → −2γ = n − 1 (Cardassian's own relation).
   Relative error at γ = 0.487: 2.4% at x = 1, 1.6×10⁻³ at x = 10, 3×10⁻⁵ at x = 100.
4. **Low density (x ≪ 1):** `ρ_DE → ρ_crit/γ` constant, w → −1. *This* is where it departs from MP-Cardassian,
   whose low-density limit is H² ∝ ρ_m^n. So: **Cardassian in the past, Λ in the future; identical to both at γ = ½.**
5. Today (C₀ = Ω_m = 0.315): x₀ = 0.953, w₀ = −0.992 at γ = 0.487 — the intermediate regime. w monotone, no −1
   crossing, for every γ tested (agrees with the 08-10 sign lock and the 09-11 `sign(1+w) = sign(1−2γ)`).

## What the literature pre-answers

- **The perturbation branch.** Koivisto, Kurki-Suonio & Ravndal (2005, PRD 71, 064027) put fluctuations into
  Cardassian models (fluid interpretation, with and without interacting DM): the late ISW is "much too strongly
  enhanced", and the model is **ruled out except in a small neighbourhood of the ΛCDM limit**. The archive's 08-18
  locality fork found its own perturbation channel is linear in ε = 2γ − 1 and permanently unpowered. These are
  consistent, and the literature result is the stronger statement for any *fluid* completion. It should be imported
  as prior art and not re-derived. **Not executed here for this specific g(ρ_m) — untested, not refuted.**
- **Background constraints.** There is a multi-paper literature of SN/BAO/H(z)/GRB fits to Cardassian and
  MP-Cardassian (e.g. arXiv:0908.1438, 1006.1105, 1706.09848). They generically land near n ≈ 0, which matches the
  08-12 γ = 0.487 result. Our DESI DR2 fit is the newest member of that series, not a first.

## What changes

1. **PREDICTIONS.md Bucket 3 DE note: add one prior-art clause** (applied as a back-annotation today, clearly marked,
   no bucket change): the sector is a Cardassian-class modified Friedmann equation, exactly ΛCDM at γ = ½, with the
   MP-Cardassian high-density asymptote n = 1 − 2γ.
2. **The TEST-26 registration (gated on dp)** should cite the class and say what a DR3 result would add beyond the
   existing Cardassian constraints. If the answer is "nothing the class literature lacks", say so in the registration.
3. **A second instance of the fit-vs-selection slip, found in passing.** /dark-energy called this sector "a bet that
   can be killed or tied at DESI DR3, but never won". The family is a one-parameter model that nests ΛCDM. If DR3 put
   the data in its allowed quadrant (w₀ > −1 with w_a > 0) at γ measurably ≠ ½, it would be *selected* over both ΛCDM
   and w₀w_aCDM. That is unlikely, not impossible. The archive corrected the same inference for the galaxy sector on
   2026-07-29 (`nested_submodel_fit_versus_selection.md`), and the correction never reached this page or /for-researchers.
   Fixed on the site today; check any archive prose that says "kill-or-tie" / "cannot win" for the same slip. (The
   TEST-26 "kill-or-tie" label is dp-gated and was not touched.)

## The frame question this raises (the reason it is a proposal, not an erratum)

**The A2ACW H1/H2 question has a data point sitting in the archive's own history.** H2 says an LLM prior-art hunter
maps *anything* onto a corpus, so it over-finds prior art. This sector went nine months, three covariant executions
and a likelihood fit without anyone (human, publisher track, or explorer) finding a 2002 PLB model class whose defining
equation it reproduces exactly. A simulated outside reader found it in one pass, once asked "is there prior art?"
while reading the equation. That points the other way: **when nobody runs the prior-art question, the program
under-finds.** The failure is that the question was never scheduled, not that it came back wrong.

**Proposed discipline (cheap):** before any *execution* on a sector — covariant completion, likelihood fit,
registration — run one scheduled prior-art pass on its defining equation as written, and record the result (including
"none found") in the sector's provenance chain. That makes the efficient path (execute) and the correct path (check
whether the class is already constrained) the same path. It is also a natural, uncontaminated **positive-control
candidate for A2ACW**: a known-prior-art item the program itself missed. It is post-hoc, so it is not a blind test.
It is still an item a working audit should flag.

## Untested vs refuted

- Cardassian-class ISW exclusion applied to *this* g(ρ_m): **untested** (imported, not executed).
- Whether a non-fluid (modified-gravity) interpretation evades it: **untested**; the 08-11 covariant completions
  already fail the background fit, so the escape room is small.
- Nothing here refutes anything new. It re-prices eight months of DE-sector work as *rediscovery within a known class*.
  That is Bucket 3, which is where the sector already sits.
