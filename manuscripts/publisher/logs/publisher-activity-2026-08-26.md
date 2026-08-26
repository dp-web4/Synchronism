# Publisher Activity — 2026-08-26 (RUN-ID 20826)

**Window**: 2026-08-25T10:00Z → 2026-08-26T10:30Z

## Scan (§1b, all surfaces)
- `Research/SESSION_MAP.yaml` — 0 new numbered sessions (Archivist: **19th** consecutive deliberate zero).
- Synchronism repo — 1 non-CI commit in window (`c67d82b3`, Archivist SESSION_MAP 2026-08-26), **plus two
  that landed 2026-08-25 15:37Z / 16:08Z, after the previous window closed at 10:30Z** and are therefore
  reviewed here: `875cfc21` (back-annotation of the site explorer's 08-25 RG finding) and `28df7b18`
  (research lane withdraws its 08-23 EFE over-reach; three-C-functions frame inscribed).
- `Research/proposals/` — 1 new: `refracted_gravity_prior_art_and_density_lever_is_log_size_20260825.md`.
- `Research/papers/`, `Research/preregistrations/` — nothing new.
- `explorations/` — 1 new (`2026-08-24-efe-zero-is-phi-independence-and-three-C-functions.md`).
- **synchronism-site** at `origin/main` — **0 commits** in window; last is `d12d711` (08-25 08:39).
  Local HEAD `main`. Maintainer lane still 401, **day 14**.
- **web4** at `origin/main` — **2 commits** (`8a0ef6e7`, `e9e1cc57`). Local HEAD `main`.
  **1 touching `whitepaper/`** (path-scoped log): `e9e1cc57`.

## Work

### 1. The external prior-art gate ran, and returned positive on the sector's constitutive equation
`∇·[C(ρ)∇Φ] = 4πGρ`, introduced in this whitepaper 2026-08-04 as a completion built in-lane and called
*"the natural one"*, is **Eq. 2.3 of Matsakos & Diaferio 2016** — *Refracted Gravity*, JCAP,
arXiv:1603.04943 — with `C` their **gravitational permittivity**. Abstract read at source.

Re-derived from the published forms, no site script imported or read
(`simulations/publisher_20260826_refracted_gravity_identity.py`, sympy + numpy, exit 0):

| leg | routing note | this lane |
|---|---|---|
| `C_Ω` vs `ε(ρ)` | max diff 2.2×10⁻¹⁶ / 8 decades | **exactly 0 in sympy**, both log bases, `p = 2q/ln(base)` |
| `B ≤ 1/Ω_m = 3.17` | "is `1/ε₀`" | confirmed; `C_Ω(ρ→0)=Ω_m`. M&D **fit** `ε₀`=0.20–0.25 ⇒ ceiling **4.0–5.0** |
| `C_ρ` vs `ε(ε₀=0)` | rms 0.014, max 0.032 | **rms 7.1e-4, max 1.4e-3** at *fitted* γ, `ρ_c/ρ_crit`≈2.07 |
| same, at γ=2 | — | 0.0142 / 0.0340 — **this is the note's figure** |
| mapped `q` | "≈1.5, ~2× RG's 0.75" | **1.14** fitted / **1.68** at γ=2; 1.5×–2.2×, criterion-dependent |

The first row is the one that matters: the identity is **closed-form**, not numerical. `tanh(log(u^q))`
through `exp` is the logistic. `C_Ω` does not approximate the permittivity — it **is** it.

### 2. Two corrections to the source, opposite directions — so this is verification, not relay
**Downward.** Routed citation "Sanna, **Pipino**, Diaferio et al. 2023" names a non-author. It is
**Sanna, Matsakos & Diaferio 2023**, A&A **674**, A209 (arXiv:2109.11217) — verified at source,
scalar–tensor, `φ` twice the permittivity in the weak field. (Pipino is on the 2022 E0-galaxy paper.)

**Upward.** The note quotes the **γ = 2** residual while its parameter discussion runs at the fitted γ.
The real match is **20× tighter**. *Under-claiming fails like over-claiming* — second consecutive day
this rule has fired against a source (08-25: a systematic band dropped mid-title).

### 3. One paper in neither source, found by running the gate rather than accepting its result
**Gervani, Diaferio, Pace & Sanna 2026** (arXiv:2601.22937, 2026-01-30): linear perturbations of a
Brans–Dicke theory that *"in the weak-field limit reduces to Refracted Gravity"*, baryons alone, no DM
and no DE — **structure formation delayed to z < 1, in disagreement with observed high-z galaxies.**
That is the published cosmological verdict on the branch §5.15 constructs as **Completion B** and
independently excludes at **Δχ² ≥ +79**.

**Written into the whitepaper as UNCHECKED, not as corroboration.** §5.15 pins `C` quasi-statically via
`C_eff(x₀)=Ω_m`; Sanna et al.'s `V(φ) = −Ξφ` is *linear*, not a mass term; §5.15 separately names a
*massive* scalar as its one unexecuted branch. They may be different theories. Settling it is REC-042's
lead blocking action. If they coincide, this archive holds a second, differently-routed constraint on a
published model — background likelihood against linear growth.

### 4. What survives, and three threads collapsing to one parameter
RG **fits** `ε₀`; this framework **derives** the floor as `Ω_m`. **0.315 vs 0.20–0.25 — ~30% apart, and
Cesare et al. 2020's 30 DiskMass galaxies discriminate them by RG's own method (vertical dispersions).**

And the section's standing tension dissolves: `C(0)=0` ("the sector's worst pathology", divergent
exterior field) and the bounded boost and RG's `ε₀` are **one parameter, the floor**. The sector shows
the pathology only because it takes its headline (EFE = 0, from `C_ρ`) from the *floorless* function and
its ceiling from the *floored* one — the 08-24 three-function frame.

### 5. Count HELD at 6
Nothing here refutes; an attribution was corrected and a comparison class named. **But two of the six**
(TEST-09/TEST-10, one inequality since 08-08) now refute a *published 2016 theory's fitted parameter*,
not a Synchronism novelty. Count unchanged; **what it is about has changed.** Bucket 0 = 0.

## Edits
- `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — one `[AMENDED 2026-08-26]` block
  closing the field-equation thread; four in-place pointers (`[ATTRIBUTION]` at the 08-04 "the natural
  one" lead; `[EXTENDED]` on the prediction-and-pathology sentence; `[PRIOR ART]` at §5.15 Completion B;
  `[scoped]` at its "one remaining unexecuted covariant branch").
- `whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md` — entry appended (CRLF matched).
- `explorations/2026-08-26-publisher-the-field-equation-has-a-name-and-the-author-was-already-in-the-repo.md` — new.
- `simulations/publisher_20260826_refracted_gravity_identity.py` — new, exit 0, all assertions pass.
- State: `recommendations.json` (**REC-042 OPENED 0.35**; REC-038 held 0.93 **instance 23**;
  REC-040 0.58→**0.62**; REC-041 0.62→**0.68**), `whitepaper_sync.json`, `whitepaper_web4.json`.

## Gates
- **Build**: `make-md.sh` exit 0 (7,869 lines), `make-web.sh` exit 0 — **re-run from `whitepaper/`**.
- **Churn, both numbers**: content **28** / raw **12,078**. Tree restored; CI builds artifacts.
  (`make-web.sh` again wanted to revert `index.html`/`style.css` by 378/206 vs CI's — 5th day running.)
- **Line endings**: `dark_matter.md` pure LF (preserved, 0 CR); `CHANGELOG.md` CRLF matched (99/55);
  exploration note and script pure LF.
- **External claims**: 3/3 verified at source. **Coverage stated**: *abstracts only* — the equation
  forms came through the routing note and were **re-derived**, not read in the papers. Reading M&D 2016
  in full is a REC-042 blocking action.
- **Not staged (not mine)**: `AGENTS.md`, `CLAUDE.md` (gitnexus, supervisor-owned),
  `simulations/session373_acceleration_regime.png`.

## Ledger
Refutation count **HELD at 6**. Bucket 0 = 0. REC-038 0.93 held, **instance 23**.
REC-036 0.45 / REC-040 **0.62** / REC-041 **0.68** / REC-042 **0.35** (new).

## Genus recorded
**The prior art's AUTHOR was already in this repository.** `simulations/session200_caustic_mass.py:69`
and `session199_mdyn_mlens_v2.py:130` cite Diaferio & Geller 1997 / Rines & Diaferio 2006 for the
caustic mass method — the name was indexed in-house before the 08-03 screen ran and declared the
neighbourhood clean. All 22 prior instances concerned an unfound **document**; this is an unfound
**author key**. Proposed gate: before declaring an external screen clean, grep the own repo for every
surname already cited in adjacent code and docs — *prior art is disproportionately written by people the
programme already reads.*

**And the same class, mirrored, in web4 on the same day** (`e9e1cc57`): both patents in
`sections/14-references` carried **real numbers and invented titles** — US 11,477,027 as "Linked Context
Token Systems and Methods", US 12,278,913 as "Trust-Based Value Exchange Protocols" — where the actual
title of both is "Apparatus and Methods for Management of Controlled Objects". web4 cited a real work
under a title matching the house frame; Synchronism did not cite a real work at all and rebuilt it. One
failure — *the house vocabulary substituted for the external record* — two repos, one day, opposite
directions. Only this track reads both surfaces.

## Self-flag
**My own build gate failed in the shape I spent the day auditing.** `make-md.sh` from the repo root
**exits 0** and prints `✅ Monolithic markdown created` after ten `⚠ directory not found` warnings,
**0 of 10** section directories found, and a **466-byte** monolith. Caught only because the reported
line count was 20 against a real 7,869. *An exit code is not a coverage claim.* Two stray directories
still exist — `Synchronism/build/` (untracked; a repo-root `git add -A` would sweep it in) and
`ai-agents/docs/whitepaper/` (outside the repo). The safety preset denies `rm` outside `/tmp` and no
`hestia_appeal` tool is exposed here, so I **complied rather than recasting** through another
interpreter. Owner action requested; commits are path-scoped so nothing stray is staged.
