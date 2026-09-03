# Publisher Daily Report — 2026-09-03

**Unit**: Windows Task Scheduler `Publisher Daily Session` (03:30 PDT), RUN-ID **15596**.
*Not* the systemd unit `autonomous-publisher-cbp` (04:30 PDT), which also commits as `[Publisher]`.
**Window**: 2026-09-02T10:30 → 2026-09-03T10:30 UTC — but see Step 0: this pass also discharges 09-02.

---

## Step 0 — sweep the previous run's tree

`git status` first, in every repo written to. Result: **the tree was already swept, by the supervisor,
not by me.** The 09-02 pass's four files were committed as `4f245b14`.

What remains uncommitted in Synchronism, none of it this lane's: `AGENTS.md` / `CLAUDE.md` (GitNexus
index counts, the supervisor's `<!-- gitnexus:keep -->` blocks), `simulations/session373_acceleration_regime.png`
(binary churn), and `build/` — see the watch item below.

**The 09-02 death is a fourth mode, and it is the dangerous one.** That pass reached Step 0, did the
web4 leg, committed it (`92c73668`), launched the SPARC reproduction as a background job, **yielded its
turn to wait on it**, and died at **exit=0** — 554 B of log, one byte off the credit-death signature —
having written **11 of 27 `ε₀` solves and no analysis**. 12.5% of one grid.

The supervisor's rescue then titled it *"last-escape reproduction run 2026-09-02 (SPARC 153 galaxies vs
MOND)."* That subject reads as a completed run. Two things generalise:

- **A rescue inherits the headline of the thing it rescues.** A sweeper can see that files are untracked.
  It cannot see that they are 12.5% of a run. Step 0 catches the tree state; nothing catches the subject.
- **exit=0 is worse than exit=1.** The 14-second exit=1 class taught every watcher that a short Publisher
  log means nothing happened. A 554 B log with a *success* code tells them, falsely, that something did.

**Rule taken and followed today: never yield the turn to a background job.** Launched the run at 03:32,
then read the site finding, rescored REC-042 and committed it, wrote the whitepaper amendment — and
collected the finished run on the way past at 03:37. Five commits landed as work completed, not at the end.

---

## Phase 0: Publication Recommendations

### New Recommendations
**None.** No new numbered sessions (Archivist RUN-ID 1656, 24 h window: **25th deliberate zero**, core arc
parked **83 days** by git).

### Status Changes

**REC-2026-042 — 0.62 → 0.68.** *Reconciling the galaxy sector with Refracted Gravity.*

The score moves **up on a day whose headline is a negative**, for reasons worth stating because they are
not obvious:

1. **A blocking action closed.** The site's seeded topic ran end to end against **six thresholds registered
   before any number existed**. The `ρ_s = +0.758` this lane measured on 09-01 is now *adjudicated* rather
   than open. A pre-registered negative is draftable; an unadjudicated correlation is not.
2. **A new publishable object, and it is the cleanest number the sector has produced** — the χ²/N = 110.79
   shape barrier (§ below). It is a statement about **Refracted Gravity**, an active external programme,
   not about Synchronism, which is exactly this item's character.
3. **It sharpens leg (1).** The identity `C_Ω ≡ ε(ρ)` is exact algebra; the field-equation *solution* is not
   MOND's curve. A draft can now say which is which.

**Unchanged**: leg (4), the striction force, remains the lead and is untouched by today — and remains the
one leg with **no in-repo computation**. That is now recorded as the gap between 0.68 and drafting.

**REC-2026-038 — HELD 0.93.** No new instance. (Today's §4 material — a rescue inheriting a headline —
is a candidate 28th instance; not counted until a second occurrence, per the standing rule.)

### Live queue
**038** 0.93 · **041** 0.68 · **042** **0.68** ↑ · **040** 0.62 · **036** 0.45 · **039** 0.38.
Two-paper strategy (REC-034 0.97 + REC-035 0.95) still pending dp.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: **Updated** — one integration.
- **Sessions reviewed**: no new numbered sessions; the input was the site explorer's 09-02 session
  (`synchronism-site` `c4b8022`).
- **Changes made**: `[AMENDED 2026-09-03]` block in `05-quantum-macro/15-dark-matter/dark_matter.md`
  (+24 lines), CHANGELOG entry (+5).
- **Proposals**: none filed — this is an amendment to a section this lane already maintains, matching the
  08-01→09-01 chain.
- **Terminology concerns**: none. `C(ξ)`, `γ`, MRH untouched.
- **Build**: verified **from `whitepaper/`** (cwd matters — a root-cwd build exits 0 having read 0 of 10
  sections). exit 0, **7,959 lines**, +29 exactly as expected.
- **Churn gate, both numbers**: content **58** vs raw **24,290** → line-ending churn, artifacts restored,
  sources only committed. CI subsequently built and deployed it (`454f4c76`: +29 monolith, +193
  `section_5.html`, PDF regenerated) — the pipeline behaved as designed.

#### The physics: the galaxy sector's last escape is closed

The 09-01 amendment named the sector's last open question as *"whether `ε₀(M_bar)` is tight enough to be
stated as a relation"* and assigned three steps to the site lane. All three ran. The answer is no, and the
reason is better than "too loose": **the correlation belongs to MOND.**

Fit the same per-galaxy `ε₀` twice — once to the data, once to **MOND's own predicted curve** for that
galaxy. Same solves, same weights, different target:

| per-galaxy constant | k (dex/dex) | σ_resid | ρ_s vs log M_bar |
|---|---|---|---|
| `ε₀` fitted to the **data** | +0.615 | 0.518 | **+0.758** |
| `ε₀` fitted to **MOND's curve** | +0.675 | 0.415 | **+0.846** |
| `a₀` fitted identically (control) | +0.063 | 0.367 | +0.071 |

MOND induces the relation **more strongly than the data show it**. Joint fit: b = +0.750 [+0.571, +0.927]
on `ε₀_MOND` and **c = +0.110 [−0.071, +0.277]** on `M_bar`, Freedman–Lane p = 0.099 — after MOND's induced
value is removed, mass carries no slope distinguishable from zero.

Four rules, thresholds fixed in advance: **R1** not a relation (0.518 vs bar 0.369); **R2** distance guard
**passes** (`a₀` slope +0.063, perm p 0.38 — the trend is real *before* it is explained away); **R3**
parameter-matched the class loses **1.75×** against a 1.20 bar (91.20 vs 52.09; wins 31% of discs,
Wilcoxon p = 2.3×10⁻⁷); **R4** MOND-induced.

**The unregistered number is the sharpest in the sector.** Every solve was scored twice, so the run also
measured how well the class can imitate MOND *at all*: **χ²/N = 110.79 against MOND's noise-free curve,
while MOND sits at 52.20 against the noisy data.** The class is further from MOND than MOND is from
reality. That is a **shape** result, and it bounds the identity this section has carried since 08-03:
`C(x)` at γ = ½ *is* MOND's simple μ to 5.6×10⁻¹⁷ — but that is the **compander algebra** and it does
**not** carry over to `∇·[C(ρ)∇Φ] = 4πGρ`. ρ-keying is not "the argument relabelled."

**Reproduced in-repo before inscription, completely**: 108 lines vs 108, **77/77 numeric lines exact**,
0 text diffs, 313 s, seeded rng so permutation p-values and bootstrap intervals match exactly.

**Also inscribed**: treatment B's residue (`c = +0.136`, real per galaxy, globally worthless — the
4-parameter optimum sets `ε₀` universal — with an untested scale-height/`Υ` nuisance R2 is blind to *by
construction*, since MOND never sees `h`); **Brouwer et al. 2021** KiDS-1000 lensing, **new to this
repository**, killing the ceiling by 6–28× even granting the relation; and two void numbers withdrawn.

**Ledger: refutation count HELD at 6; Bucket 0 = 0.** Not a seventh refutation — nothing here refutes a
*registered prediction*. What closed is the last **candidate** by which Bucket 0 might have become 1.

### Web4 Whitepaper
- **Status**: **Current**. Scanned `origin/main`; **HEAD was also `main`** (stated per the 08-09 rule), fetched first.
- **Repos checked**: `web4` — 16 commits in the window since `557b61e3`; **2 touching `whitepaper/`**, and
  **both are the other Publisher unit's own run log** (`PUBLISHER_CONTEXT.md`) plus the `make-web.sh`
  guard. Verified **by file list, not by subject**. **Zero spec or content commits.**
- **Proposals**: none. **Terminology concerns**: none.
- Recorded as *their* finding, not mine: `023d6022`'s subject reports a stranded 08-30 log entry meant
  "two passes had re-confirmed a census figure the missing entry refuted" — the same failure class this
  lane measured on 08-29, observed independently in a sibling repo.

---

## The methodological finding — and it is about this lane's own instrument

Today's reproduction scored **77/77 exact**. Line 41 of both outputs reads
`rho_c at top edge 59%`, and **that number is void**: column 1 of `epsilon0_per_galaxy_fw.npy` is
per-galaxy χ², not `ρ_c`. Verified at the **writer**, not accepted from the correction —
`epsilon0_free_the_ceiling_rescue.py:220` saves `np.vstack([e0s, best_per_gal, MOND_F, NN])`, and the
column settles it unaided: median **49.4**, max **1.17×10⁴**, against a `ρ_c` grid of 10⁻⁶–10⁻² M☉ pc⁻³.
Four orders out, dimensionally wrong.

Both runs read column 1. Both printed it. The comparator scored it exact.

> **A reproduction gate certifies that the arithmetic ran the same way. It is blind by construction to a
> mislabelled input, because that defect lives upstream of the arithmetic.**

Tighter tolerance, more digits, another machine — all would have scored it exact. To catch this class you
must read the **writer** of an artifact, not the reader.

This matters because reproduction is this lane's principal instrument and was being asked to carry more
than it can. On 09-01 this lane wrote *"reproduced 11 of 11 before inscription"* and treated that as the
gate licensing a whitepaper amendment. **That gate would not have caught this either.** The 11/11 was
true, was worth running, and was never evidence about provenance.

Proposed standing rule: **reproduction is blind to provenance — agreement bounds the arithmetic, never the
labels.** The existing rules cover disagreement (*a check that contradicts a proof is the suspect*) and
under-claiming (*verify declared properties by computing*). None covered *agreement that certifies nothing*.

**The same discipline, run against the leg I liked.** The site kills the ceiling on Brouwer+2021 at
**13–87×** — an attractive result for this lane. Fetched the abstract before inscribing: it says the
measurement *"extends the radial acceleration relation by **2 decades**."* Two decades below SPARC's floor
is `g_bar ≈ 10⁻¹⁴`, not the `10⁻¹⁵` the 87× column assumes. **6–28×** went into the whitepaper, with
coverage stated inline. *A coverage caveat is a work item* fired against me on 08-28; it only means
something when it costs something.

---

## Watch items opened

1. **Treatment B's `c = +0.136` is the sector's only open thread.** Its nuisance channel (`ρ_mid ∝ Σ/h`)
   is one MOND *cannot* guard against, because `g_bar` never sees the scale height. One afternoon:
   refit with `h × {0.5, 2}`, `Υ ∈ {0.3, 0.7}`; does `c` move by its own width? Site lane's to run.
2. **Brouwer+2021 is inscribed at abstract-level coverage.** A full read would settle the actual `g_bar`
   floor the deficit table is indexed on.
3. **Repo-root `build/`** has held a 466 B near-empty monolith, untracked, since 2026-08-26 03:43 — the
   signature of a `make-md.sh` run from the repo root (exit 0, 0 of 10 sections). Harmless in itself, but
   it is a **permanent false positive in Step 0's own signal**, and Step 0's premise is that untracked
   files mean stranded work. A standing false positive trains the sweep to ignore the signal — which is
   exactly how the 08-29 loss stayed invisible for 25 hours. **Proposed: `/build/` in `.gitignore`** —
   dp's call, same class as the `.gitattributes` eol fix. Not deleted by this lane.

---

## Summary

The galaxy sector's last escape is **closed by execution**: `ε₀(M_bar)` is not a relation of the theory but
what MOND induces — and MOND induces it *more strongly* (+0.846) than the data show it (+0.758), leaving
mass no slope distinguishable from zero once MOND's value is removed. Adjudicated against six thresholds
registered before any number existed, reproduced in-repo 77/77 exact, and inscribed in §15-dark-matter
along with the sharper unregistered result that a universal-knee density-keyed class sits **further from
MOND's clean curve (110.79) than MOND sits from noisy data (52.20)**. Refutation count **HELD at 6**;
what closed is the last candidate for Bucket 0 = 1, not a seventh refutation.

The finding I expect to outlive the physics is smaller and closer to home: **the reproduction was exact
and reprinted a void number.** This lane has spent six weeks building gates on the correct premise that
gates beat assertions. Today one returned a perfect score on a quantity that was not the quantity named —
not because it ran badly, but because agreement was never the kind of evidence it was being asked for.
