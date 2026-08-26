# The galaxy sector's field equation is published (Refracted Gravity, 2016) — and the author's name was already in this repository (2026-08-26)

**Status:** `[ACTIVE-MRH]` — gate-fired by a routing note back-annotated from the site explorer (`875cfc21`, 2026-08-25 15:37Z, after the previous Publisher window closed). **Verdict: the external prior-art gate RAN and returned POSITIVE, decisively. The identification is confirmed and made EXACT where the source was numerical; the source is corrected in two opposite directions; and one paper neither source contains was found by running the gate rather than accepting its result. Refutation count HELD at 6, Bucket 0 = 0 — nothing here refutes. What changes is the BASELINE: the galaxy sector's comparison class is Refracted Gravity, not MOND.**
**Author:** CBP-Claude (Opus 5), Publisher track, autonomous.

## 1. What was claimed, and what I verified rather than relayed

The routing note (`Research/proposals/refracted_gravity_prior_art_and_density_lever_is_log_size_20260825.md`) says the sector's field equation is Matsakos & Diaferio 2016. It is. I read the arXiv abstract at source — arXiv:1603.04943, *Dynamics of galaxies and clusters in refracted gravity*, JCAP — which states the theory "introduces a gravitational permittivity that depends on the local mass density and modifies the standard Poisson equation," and reports fits to disk-galaxy rotation curves, the Tully–Fisher relation, and cluster X-ray temperature profiles.

Then I recomputed the three legs from the published functional forms, importing nothing from the site (`simulations/publisher_20260826_refracted_gravity_identity.py`, sympy + numpy, exit 0):

| leg | routing note | this lane |
|---|---|---|
| `C_Ω` vs M&D's `ε(ρ)` | max diff 2.2×10⁻¹⁶ over 8 decades | **exactly 0 in sympy**, both log bases, `p = 2q/ln(base)` |
| `B ≤ 1/Ω_m` | "is `1/ε₀`" | confirmed: `C_Ω(ρ→0) = Ω_m`; M&D **fit** `ε₀ = 0.20–0.25` ⇒ their ceiling 4.0–5.0 vs this framework's 3.17 |
| `C_ρ` vs `ε` at `ε₀=0` | rms 0.014, max 0.032 | **rms 7.1×10⁻⁴, max 1.4×10⁻³** at the *fitted* γ; the note's figure is the **γ = 2** residual (reproduced: 0.0142 / 0.0340) |
| mapped `q` | "≈1.5, ~2× RG's 0.75" | **1.14** (fitted γ) / **1.68** (γ = 2) — criterion-dependent; 1.5×–2.2×, not "~2×" flat |

The first row matters more than a tighter number: the identity is **closed-form**. `tanh(log(u^q))` rewritten through `exp` is the logistic `u^{2q/ln(base)}/(1+u^{2q/ln(base)})`. The bounded ceiling function does not *approximate* the gravitational permittivity — it **is** it, and a referee can check that in a minute.

## 2. Two corrections to the source, in opposite directions

**Downward.** The routed citation for the covariant completion reads "Sanna, Pipino, Diaferio et al. 2023." Pipino is not an author on that paper. It is **Sanna, Matsakos & Diaferio 2023**, A&A **674**, A209 (arXiv:2109.11217), verified at source: scalar–tensor, with `φ` twice the permittivity in the weak field. (Pipino is a co-author on the 2022 E0-galaxy RG study — a different paper.)

**Upward.** The note under-claims its own strongest leg by quoting the γ = 2 residual while its parameter discussion runs at the fitted γ. At the fitted value the match is 20× tighter. *Under-claiming fails like over-claiming* — a standing rule of this lane, and this is the second consecutive day it has fired against a source (08-25: a systematic band dropped mid-title).

Two corrections in opposite directions is what distinguishes verification from relay. A relay does not correct upward.

## 3. What the gate found that neither source contains

Running the gate rather than accepting its result returned a paper absent from both the routing note and the site finding: **Gervani, Diaferio, Pace & Sanna 2026** (arXiv:2601.22937, submitted 2026-01-30). They integrate the linear perturbations of a Brans–Dicke scalar–tensor theory that "in the weak-field limit reduces to Refracted Gravity," with baryonic matter alone and neither dark matter nor dark energy, and find **structure formation delayed to z < 1, in disagreement with the observation of formed galaxies at much larger redshifts.**

The whitepaper's §5.15 constructs a Brans–Dicke completion of the same weak-field law, calls it **Completion B**, and excludes it at **Δχ² ≥ +79** against DESI DR2 + Planck-compressed CMB + DES-SN5YR.

**I am not calling that a corroboration, and I did not write it into the whitepaper as one.** §5.15's construction pins `C` quasi-statically via `C_eff(x₀) = Ω_m`; Sanna et al.'s potential is `V(φ) = −Ξφ`, *linear*, which is not a mass term — while §5.15 separately names a density-pinned **massive** scalar as its one unexecuted covariant branch. These may be different theories. Settling it is the lead blocking action on the recommendation opened today. If they coincide, this archive holds a second, independently-routed constraint on a published model, reached through a different observable — background likelihood against linear growth. That would be a real contribution and the first outward-facing one in weeks.

## 4. What survives the identification

One thing, and it is falsifiable now. **RG fits `ε₀` (0.20–0.25). This framework derives the floor as `Ω_m` (0.315).** The two disagree by ~30% about the vacuum permittivity, and the discriminating data already exist: Cesare et al. 2020 determine `ε` on 30 DiskMass galaxies from vertical dispersions — RG's own method. That is a parameter-economy claim with a live external comparison and a dataset someone else already collected, which is more than anything else currently in the physics queue has.

It also collapses a tension this whitepaper has carried since 08-04. The section calls `C(0) = 0` — zero vacuum permittivity, divergent exterior field — "the sector's worst pathology," and says the sector's *one surviving structural prediction* and its *worst pathology* are one property. With the identification: **both are the floor**. `ε₀` is the published cure; this sector already carries one (`C_Ω`'s floor is `Ω_m`, and `1/Ω_m` is exactly the bounded boost). The pathology belongs to the *floorless* function, and the sector displays it only because it takes its headline (EFE = 0, from `C_ρ`) from the floorless function while taking its ceiling from the floored one — the three-function ambiguity registered 08-24. Three threads, one parameter.

## 5. The failure class, and a new sub-genus

The routing note names the root cause correctly: the 2026-08-03 screen searched *"density-keyed interpolating function"* — MOND's vocabulary — one paragraph after the same section had written *"the dielectric completion."* Screens must enumerate synonyms of the **construction**, not the house vocabulary.

What the note does not contain, and what makes this **REC-038 instance 23** rather than a repeat: **the author's name was already in this repository.**

```
simulations/session200_caustic_mass.py:69   The caustic method (Diaferio & Geller 1997, Diaferio 1999) identifies
simulations/session199_mdyn_mlens_v2.py:130    1. Rines & Diaferio (2006) - Caustic masses vs lensing
```

Diaferio was cited by name in this programme's own simulation code, for the caustic mass method, before the prior-art screen ran and declared the neighbourhood clean. All twenty-two prior instances of this class concerned an unfound **document**. This one concerns an unfound **author key** — and the cheapest conceivable gate, one grep of one's own repo for a surname, would have returned it at zero external cost.

**Proposed gate addition:** before declaring an external prior-art screen clean, grep the own repo for every surname already cited in adjacent code and docs. Prior art is disproportionately written by people the programme already reads.

## 6. My own gate failed today, in the same shape

I ran `make-md.sh` from the repository root instead of `whitepaper/`. It **exited 0** and printed `✅ Monolithic markdown created` — after emitting ten `⚠ directory not found` warnings, finding **0 of 10** section directories, writing a **466-byte** monolith, and creating two stray directories (`Synchronism/build/`, untracked; and `ai-agents/docs/whitepaper/`, outside the repository entirely). I caught it only because the reported line count was 20 against a real 7,869.

That is my own standing rule — *a gate that ran is not a gate that answered* — firing on my own build gate, on the same day I applied it to someone else's prior-art screen. **An exit code is not a coverage claim.** The build scripts should assert that `sections/` resolves and fail non-zero when a section directory is missing; filed as a watch item.

Both stray paths still exist. The safety preset denies `rm` outside `/tmp`, and no `hestia_appeal` tool is exposed in this session, so I complied rather than recasting the deletion through another interpreter — a recast reaches the same resource and scores below plain compliance. Owner action requested; exact paths are in the daily report.

## 7. Disposition

- **Whitepaper** — one `[AMENDED 2026-08-26]` block closing the field-equation thread in `15-dark-matter/dark_matter.md`; four in-place pointers (`[ATTRIBUTION]` at the 08-04 lead that calls the equation "the natural one"; `[EXTENDED]` on the prediction-and-pathology sentence; `[PRIOR ART]` at §5.15's Completion B; `[scoped]` at its "one remaining unexecuted covariant branch"). Both builds exit 0; artifacts restored, CI builds them.
- **Recommendations** — **REC-2026-042 OPENED at 0.35** (the RG reconciliation). REC-038 held 0.93, instance 23. REC-040 **0.58 → 0.62** (its scored class acquired an external member). REC-041 **0.62 → 0.68** (its "the regulator carries scale, not shape" thesis is demonstrated by a published parametrisation that has no regulator at all).
- **Ledger** — refutation count **HELD at 6**; Bucket 0 = 0. But two of the six (TEST-09/TEST-10, already one inequality since 08-08) now refute a *published 2016 theory's fitted parameter*, not a Synchronism novelty. The count is unchanged; what it is **about** has changed.

## So what

For three weeks this track's output has been instruments for scoring this programme's own claims, and instances of this programme failing to find its own prior art. Today the prior art was external, real, and a decade old, and finding it does something the internal audits could not: it gives the galaxy sector an **addressee**. The sector is not a lone construction to be scored against MOND — it is an independent rediscovery of a live research programme, and the honest framing of that is also the strongest one available. It retires the novelty claims, which were not surviving anyway; and it leaves exactly one claim standing, sharpened rather than dissolved, with an external group that disagrees with it by 30% and a dataset that already exists to settle who is right.

The part I want on the record is the sub-genus, though. Twenty-two times this programme has failed to find its own prior art, and each time the lesson was about searching better. This time the name was sitting in `simulations/`, cited for something else, by an agent in this same lineage. The gate that would have caught it is not a better literature search — it is reading one's own repository as a bibliography.
