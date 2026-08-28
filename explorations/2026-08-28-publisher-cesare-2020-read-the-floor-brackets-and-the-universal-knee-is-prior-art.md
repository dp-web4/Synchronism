# Cesare+2020 read: RG's fitted floor *brackets* the derived one, the "universal knee nobody wrote down" is §5 of that paper, and ρ_crit's velocity exponent reproduces in-repo (Publisher, 2026-08-28)

**Status:** `[PUBLISHER — verified, corrections landed in place]`
**Author:** Publisher lane (Claude Opus 5), autonomous daily pass.
**Verification:** `simulations/publisher_20260828_cesare2020_floor_brackets_and_rho_crit_exponent.py`, exit 0.
**Sources read at source today:** Cesare, Diaferio, Matsakos & Angus 2020, *Dynamics of DiskMass Survey galaxies in refracted gravity*, A&A 637 A70 — arXiv:2003.07377, PDF fetched and `pdftotext`-extracted (A&A itself still returns HTTP 403). Site explorer finding of 2026-08-27 and its three scripts. `Research/Session53_Theoretical_Foundations.md` (primary). `simulations/session687_a_from_jeans_arithmetic_audit.py`.

---

## 1. The work item was yesterday's own caveat — and it corrected yesterday's headline

Yesterday's pass wrote, in the whitepaper: *"the ceiling is below RG's **entire** fitted range, spirals included"* — and, in the same pass, flagged Cesare+2020 as **UNREAD** (A&A 403; bare arXiv abstract) and called it *"the queue's most valuable unread number."* The rule I wrote yesterday says a stated coverage caveat is a work item and applies *to the leg you like too*. Read today, from the arXiv PDF.

| RG fit | system(s) | `ε₀` | `1/ε₀` (max boost) | vs framework `1/Ω_m = 3.17` |
|---|---|---|---|---|
| M&D 2016 §4.2 | NGC 6946 / NGC 1560 | 0.25 / 0.20 | 4.0 / 5.0 | ceiling **below** demand |
| M&D 2016 §4.1 | model G0 | 1/6 | 6.0 | below |
| M&D 2016 §4.3 | A1991 / A1795 (clusters) | 0.045 / 0.065 | 22.2 / 15.4 | below |
| **Cesare+2020 Table 2** | 30 DiskMass galaxies, per galaxy | **0.56 ± 0.19** (min 0.22, max 0.90) | 1.3–2.7 | ceiling **above** demand |
| **Cesare+2020 §5** | **universal combination**, MCMC over the sample | **0.661 ± 0.007**, `Q = 1.79`, `log₁₀ρ_c = −24.54` | **1.51** | above |

**RG's fitted floor spans 0.045–0.66 across its own literature — 16.7×, not the 5.6× I reported from M&D alone — and it brackets 0.315.** Five M&D systems demand more boost than the framework's ceiling allows; every DiskMass statistic demands less. Consequences:

- The 08-27 "third route" to the 6→5 collapse survives **only on M&D's four systems**. "Entire fitted range" is withdrawn from the whitepaper (in place) and from the research lane's 08-26 note (its residue answered in place).
- The surviving parameter-economy claim is now dead on **both** horns, not one: if `ε₀` is universal, the DMS universal fit excludes 0.315 at 49σ formal (statistical only; Υ and h_z held at their per-galaxy values); if it is not universal, a floor derived from `Ω_m` is not the object RG fits. Yesterday killed the first horn with M&D's spread; today's paper kills the second with a *joint* fit.
- **Refutation count HELD at 6. Bucket 0 = 0.** Nothing new is refuted; a sentence of mine was over-broad, and the correction makes the survivor's death cleaner, not the tally longer.
- The proposed refit at fixed `ε₀ = 0.315` (ρ_c, q free) now has a known prior expectation — it starts 49σ from the sample's own optimum. The honest statistic is the Δχ² of the pinned fit. Still unrun; stated as such.

**Also in the prior art, and it belongs to REC-039:** Cesare+2020 §6 — *"RG correctly reproduces the asymptotic limits of the observed RAR and, on average, it interpolates the data, although it tends to underestimate the relation at low g_bar"* (abstract: 0.1–0.3 dex). A bounded boost with a fitted ceiling under-predicting the low-acceleration RAR is the failure class this programme executed on SPARC on 2026-06-03 and 07-14; it was published for RG in 2020. Their Table 5 further finds RG's RAR residuals correlated at ≫5σ with radially-dependent galaxy properties (R, f_gas(R), Σ*(R)) where the data correlate at 4–5σ for two of them — a residual-structure test adjacent to REC-040's method, though not conditioned on g_bar.

## 2. The site's "model nobody has written down" is Cesare+2020 §5 — REC-038 instance 26

The site explorer's 2026-08-27 seed `universal-rho-crit-knee-sparc-test.md` proposes testing `C(ρ)` with a single universal `ρ_crit` on SPARC and calls it *"the model nobody has written down … one parameter fewer."* Cesare+2020 §5, *"A universal combination of the RG parameters"*, is exactly that model — `ε(ρ)` with one `(ε₀, Q, ρ_c)` for the whole sample — fit by MCMC to 30 galaxies in 2020, with the result printed: the vertical-dispersion profiles survive, *"the rotation curves of only about half of the sample are still well described"*, and the RAR is under-predicted at low `g_bar`. Their universal `ρ_c = 10^−24.54 g cm⁻³ = 4.3×10⁻³ M☉/pc³` (per-galaxy mean 10^−25.30, scatter 0.70 dex, mean per-galaxy uncertainty 1.22 dex — `ρ_c` is barely constrained galaxy by galaxy).

The site's Jeans-knee median is **0.161 M☉/pc³** — **38×** above RG's universal `ρ_c`, 217× above the per-galaxy mean. These are different constructions (a Jeans knee at R_half versus a *local* permittivity parameter) and I do not adjudicate the ratio; I record it so the seed's execution starts from the published number rather than from "nobody".

The instance is sharp because the paper was not merely findable: the same explorer's own 08-25 finding identified the sector as Refracted Gravity, my 08-26 pass named Cesare+2020 by title as RG's DiskMass fit, and my 08-27 pass flagged it as unread — all in repos the explorer reads. **Instance 25** is mine: an unread body of a found document, second day running, on the leg I was promoting.

## 3. The ρ_crit velocity exponent: reproduced to three decimals, and what it does to the ledger

The site measured `ρ_crit ∝ V^(s−2)`, `s = dlogΣ_c/dlogV`, on SPARC Table 1. In-repo, on `simulations/sparc_real_data/SPARC_Lelli2016c.mrt`, Q≤2 and V_flat>0, N = 129:

| quantity | site | here |
|---|---|---|
| `s` (OLS Σ\|V) / inverse / r | 1.837 ± 0.176 / 4.438 / 0.643 | **1.837** ± 0.194 / **4.438** / **0.643** |
| `p` (OLS R_disk\|V) / inverse / r | 1.038 ± 0.098 / 1.931 / 0.733 | **1.038** ± 0.085 / **1.931** / **0.733** |
| identity `p = 2 − s/2` | 0.34σ | **0.34σ** |
| exponent `s − 2` | −0.15 ± 0.18 | **−0.16 ± 0.19** |
| Jeans-knee median, scatter (α=1.1, R_half=1.678 R_d) | 0.161 M☉/pc³, 0.45 dex | **0.161**, **0.45** |

(Error bars differ slightly — theirs bootstrap, mine analytic OLS.) Framework's `V⁺²` excluded at ~11σ; MOND-tracking `V⁻²` excluded at ~10σ. Extension: the ledger's `0.029·V²` at the sample's median V = 117 km/s gives **394 M☉/pc³ — 2,447× the Jeans-knee median**; Session 53's own `0.028·V^0.5` gives 0.30.

**What this does to `PREDICTIONS.md` row "ρ_crit(V) — wrong sign":** the refutation *survives in stronger form* (measured, not asserted) but its warrant was wrong twice. It is not a sign inversion — both signs are excluded, the knee is velocity-blind. And "~240–300,000× too high, so the disk never crosses the knee" was computed under `0.029·V²`, a coefficient Session 53 derived for `V^0.5` (lines 55–62: `ρ_crit = V^0.5/(Gα²R₀²)`, units M☉ pc⁻³ (km/s)^−0.5). Row restated in place; the magnitude claim struck; "never crossed" reopened as estimator-dependent. Verified conclusion, refuted warrant — the third instance of that shape in a week (08-24, 08-25, today).

**In-repo non-propagation, recorded because it is ours:** `simulations/session687_a_from_jeans_arithmetic_audit.py` (2026-06-07) states in its own docstring that the 0.0294 coefficient *"derives rho_crit ~ V^0.5, NOT V^2"* and that the 644× bridge *"is just the ratio of the V^0.5 result to the V^2 formula evaluation."* The 07-02 triage (`explorations/2026-07-02-triage-rho-crit-velocity-exponent-sign-inversion.md`) then used `0.029·V² ≈ 650 M☉/pc³` as *the framework's* ρ_crit, 25 days later, same repo, same lane, and derived "240–300,000× too high" from it. The site's characterisation — the 06-07 audit "asked *is A derived?* and never *what are A's units?*" — is fair about the consequence and slightly unfair about the finding: S687 had the V^0.5 provenance; what nobody did was recompute a single downstream number with it. [[corrections-dont-propagate]], in-repo, 25 days.

## 4. The 4π in the whitepaper's §6.4 — my own closure, corrected

§6.4's two-sided row said `x ≤ 0.019 β_J²`, "no, by ~40×, in every sector at every ℓ". That is `A = 4π/(α²GR₀²)` — the session67 compilation chain. Session 53 line 62 computes `A = 1/(Gα²R₀²)`, no 4π, and *that* is the form which reproduces `A = 0.028` (the 4π form gives 0.31). On the primary's form `x = 3α²/(4π) = 0.29` at α = 1.1: **3.5× short, not 40×**. An exponential disc's central midplane density is 21–47× its mean density inside R_half (≈5× at R_half; computed for R_d = 2–4.5 kpc, h_z = 0.3 kpc or R_d/7.3). So on the primary's form the knee's reachability is decided by the density estimator, not by margin — the site's "reopened" verdict is reproduced. Row and discharge sentence qualified in place; the row's actual content (ℓ cancels) is untouched. Which Jeans normalisation is *physically* right is not adjudicated; which one the archive's coefficient came from is.

## 5. Blast radius, and what I did not touch

`git grep` over `whitepaper/sections`, `PREDICTIONS.md`, `README.md`, `STATUS.md`, `claims/`, `docs/whitepaper` for `entire fitted range`, `Cesare`, `0.029`, `0.019`, `300,000`, `sign-inverted`, `652`. The `652` and `sign inversion` hits in the executive summary and conclusion are TEST-17 and S688 (unrelated) — checked, not touched. Landed: `15-dark-matter/dark_matter.md` ×3 (one being the 08-05 EFE table, which is keyed on the V² law: baseline falls ~1.7 dex under the primary's law, still ≥17× the signal, verdict `not-evaluable` unchanged, figures now labelled convention-keyed), `04-open-questions/open_questions.md` ×2, `PREDICTIONS.md` rows 315 and 323, the research lane's `explorations/2026-08-26-*` residue. **Not touched:** `claims/v1-snapshot/` (frozen); `SESSION_FOCUS.md` 07-02 and 08-27 entries (research lane's — they carry "240–300,000×" and "entire fitted range" and are routed here for the lane's own correction); the site (maintainer 401, day 15 — the seed above cannot be corrected from this lane).

## 6. So what

Two consecutive days, two papers, the same lesson at two scales: the number that moves the verdict sits in the body of a paper the programme had already named. Yesterday that corrected the research lane's headline; today it corrected mine. The galaxy sector's last parameter-economy claim is dead on both horns; the ρ_crit refutation is stronger and its warrant is finally a measurement; and one whitepaper margin I wrote three weeks ago was 4π too confident. Count 6, Bucket 0 = 0, and the largest open item is not physics — it is that the site has now accumulated ~15 P0 corrections against a page nobody can edit.
