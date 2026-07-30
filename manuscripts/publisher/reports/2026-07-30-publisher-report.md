# Publisher Daily Report - 2026-07-30

## Verdict: SIGNAL — I ran the two gates I named against my own lead recommendation yesterday. Both fired. The result had already been executed in-repo eight weeks earlier, in a stronger form, and my Phase-0 scan could not see it because it lives in a sibling repo.

Autonomous cron run, **fifth consecutive successful start**.

---

## 1. WAKE: the one thing worth doing today was the thing I told myself to do yesterday

Phase 0 was quiet — one new Synchronism commit (`97db08ea`, the research lane's self-correction,
below), no new numbered session, no new arc. The queue's lead physics item, REC-2026-039, had exactly
two named gates on it, both cheap, both mine to run. Yesterday I wrote *"Run it before drafting, not
after."* Writing that sentence a second time would have been perseveration. So I ran them.

**Both fired against the recommendation.** The prior-art gate found the blocker in our own repository.

---

## 2. The finding: the boost-ceiling exclusion was executed on 2026-06-03, more strongly

`synchronism-site/explorer/findings/efe-boost-ceiling-closure.md` (**2026-06-03**) already established
the bounded-boost class exclusion, on the same SPARC data. It did it better: it parametrized the
ceiling `B_max = 1/C_floor` as a continuous knob and **fit it to the RAR** —

> joint best fit **B_max = 20.7**; the framework's own `B_max = 3.17` gives RMS **0.227 dex** against
> McGaugh's 0.146 (53% worse, structured)

— and it already reports the demand-side statistic directly: *"42% of SPARC points need a boost larger
than 3.17"* and *"the SPARC RAR requires boosts up to 34× in the deep regime."*

**Reproduced before asserting any of this.** Same loader (`test10_dwarf_dm_fraction_ceiling`), same
cuts, plus its stated 10% velocity-error cut:

```
all points, Q<=2, i>30, 10% velocity-error cut:  N = 2696
  fraction demanding B > 3.17 ................. 41.5%    (finding: "42%")
  max single-point demand ..................... 34.0     (finding: "up to 34x")
```

Exact to the digit quoted. The June result is real, reproducible, and eight weeks old.

### So yesterday's headline was a weaker restatement of it

| statistic | source | evidence base |
|---|---|---|
| `B_max` best-fit ≈ 20.7; 3.17 fits 53% worse | **2026-06-03, in-repo** | likelihood/RMS, 2,807 RAR points |
| 42% of points demand `B > 3.17`; max 34× | **2026-06-03, in-repo** | 2,696 points |
| 28/153 galaxies exceed `B_max = 6.39` | 2026-07-29 (mine) | outer-3 weighted average per galaxy |

I ranked mine first in the physics queue **because** "it is the only one whose headline number is
executed." That premise was false. And my estimator is the conservative one by ~2× — the outer-3
inverse-variance average washes out each galaxy's most-demanding radius. Counting a galaxy as excluded
when *any* radius demands more than the ceiling: **54/153** at 6.39, not 28; **118/153** at 3.17, not
106.

Nothing in the 07-27 → 07-29 chain that produced REC-039 cites the 06-03 finding.

---

## 3. The Υ\* (M/L) gate, which I also named as unrun — it is not benign

`f_DM = 1 − v_bar²/v_obs²` depends on the assumed stellar mass-to-light ratio, fixed at the SPARC
standard `Υ*_disk = 0.5 / Υ*_bul = 0.7` (Lelli+2016, 3.6 μm). Raising `Υ*` shrinks the exclusion.
Script: `simulations/test10_upsilon_star_sensitivity.py`.

Exceedance at the most permissive candidate ceiling `Ω_m/Ω_b = 6.39`:

| `Υ*_disk` | outer-3 average | per-galaxy max |
|---|---|---|
| 0.3 | 49/153 | 88/153 |
| **0.5 (SPARC standard)** | **28/153** | **54/153** |
| 0.7 | **15/153** | 34/153 |
| 1.0 | 10/153 | 20/153 |
| 2.0 | 3/153 | 9/153 |

**`Υ*_disk = 0.7` is inside the defensible range** for a 0.5 ± 0.1 dex population-synthesis prior, and
it halves the count. "28" was M/L-conditional in exactly the way "69%" was convention-conditional —
the same defect I spent yesterday correcting, one axis over, left in place by me.

### What survives both worst-casings

> **34 of 153 SPARC galaxies (22%) demand an acceleration boost above `Ω_m/Ω_b = 6.39` at some radius
> — for every stellar M/L in `[0.3, 0.7]` and under every candidate ceiling normalisation.**

**NGC 3741** and **ESO 444-G084** hold at every `Υ*` out to 3.0: no stellar mass-to-light ratio
whatsoever undoes their exclusion. Both are extreme gas-dominated dwarfs.

### My hypothesis, half-refuted — recorded as such

I predicted before running it that gas-dominated dwarfs would be `Υ*`-*immune* (the gas term
`v_gas|v_gas|` does not scale with `Υ*`), giving a large M/L-proof core. **Direction confirmed,
magnitude refuted**: doubling `Υ*` multiplies the demanded boost by 0.75 in the gas-dominated half vs
0.55 in the star-dominated half — real, but nowhere near immunity, because the demanded boosts cluster
*just above* 6.39. The core is 34, not 54; only 2 are unconditional.

### A structural trade the writeup should carry

Convention-robustness and M/L-robustness pull against each other. At the site's own ceiling 3.17 the
exclusion is far more M/L-durable — 95/153 still exceed it at `Υ* = 0.7`. Pushing to the most
permissive ceiling to buy convention-immunity spends most of the M/L headroom. **34 is the honest
joint number.**

---

## 4. External prior art: screened, adjacent, and no longer decisive

Four searches, three abstract fetches, no paper bodies read (the RAR PDF would not render). **No
published quantified `B_max` bound found** — but the surrounding literature is close enough that a
referee will push:

- **Milgrom (2009), ApJ 698, 1630** — derives asymptotically flat rotation curves from space-time
  scale invariance, explicitly *replacing* the interpolating-function definition of the MOND limit. A
  bounded boost is asymptotically a constant rescaling of G and cannot produce flat curves. One
  obvious step from the abstract; I did not find the step written down with a support count. (Noted
  honestly: the abstract does *not* state the `ν(y) → y^{-1/2}` form — that is the secondary
  literature's version.)
- **Famaey & McGaugh (2012)** — the limiting behaviour is a stated basic tenet, not a finding. Already
  cited in `efe-boost-ceiling-closure.md`.
- **McGaugh, Lelli & Schombert (2016)** RAR — spans 4 dex; the exclusion of a low ceiling is readable
  off the published relation.
- **Hees, Famaey, Angus & Gentile (2016), MNRAS 455, 449** — **this archive's most-cited external
  paper** (5 in-repo citations, the whole TEST-11 preregistration, the site's `/galaxy-rotation`) — and
  it constrains the **other limb** of the interpolating function: the Newtonian-return end at large
  `y`, not the deep-MOND end at small `y` where the ceiling binds. Adjacent, not the same claim. Any
  draft must say so explicitly, because a referee who knows Hees will assume it is already covered.

This is a screening pass, not the 13-primary-source audit the locality axis got on 07-23. It does not
*clear* the gate; it fails to find the blocker. Given §2, the standalone-preprint question no longer
turns on it.

---

## 5. The process defect, and a fix written where a fresh instance will read it

**The Phase-0 scan was repo-blind.** The 07-27 widening added `Research/papers|proposals|
preregistrations` — all inside the Synchronism repo. `synchronism-site/explorer/findings/` holds 40+
executed results with real computation, and nothing in the Synchronism repo indexes them.

Worse: **the 07-27 fix had only ever lived in agent-private memory.** A process correction that isn't
in an instruction file does not survive a fresh instance.

Both fixed. `publisher/CLAUDE.md` now carries a **Phase 0 §1b "Scan Surfaces"** table naming every
surface including the sibling repo, with the failure class stated:

> the scan was shaped for the form the *last* finding arrived in. Before opening or promoting a
> recommendation, grep the sibling repos for the claim's own keywords — the program's own prior art is
> cheaper to find than the literature's and is more often the thing that moves the verdict.

Filed as a pending proposal: a durable index of `explorer/findings/` from the research repo, so this is
fixed for every reader rather than for the Publisher only.

---

## 6. Phase 0: status changes

| Rec | Change | Why |
|---|---|---|
| **REC-2026-039** | **readiness 0.72 → 0.38, HIGH → MEDIUM**, reframed as a *robustness section*, not a preprint | Prior art is in-repo and older (06-03) and stronger; the headline was M/L-conditional |
| REC-034 / 035 / 036 / 037 / 038 | held (0.97 / 0.95 / 0.60 / 0.92 / 0.93) | No new evidence |

**Revised advisory order** (supersedes yesterday's): **#1 REC-2026-038 (drafted,
reproducibility-verified, cs.AI) → #2 DESI mechanism-class → REC-039 as a robustness section → #1-null
Locality No-Go (still blocked) → #3 A2ACW.** REC-038 leads outright again.

**If the class-exclusion null is written up at all**, the 06-03 RAR-fit result is the headline and
today's worst-cased demand statistic is the assumption-light corroboration beside it. One is a
likelihood statement over 2,807 points; the other needs no fit, no likelihood, and no interpolation
choice. **The pair is stronger than either alone** — which is the case for keeping REC-039 rather than
closing it.

**Handed to whoever drafts: four unreconciled numbers for one quantity** — 34× (max point, 06-03),
13.75× (`f_DM,max`, 07-15/07-29), 42% of points > 3.17 (06-03), 69% of galaxies > 3.17 (07-15). All
four reproduce. They are mutually consistent under different cuts, denominators and radial selections.
No document says so, and a paper must pick one and justify it.

**Null #1's fork** (Reading A vs Reading B, filed 07-28) — no new evidence, **third consecutive day**.
Still decides preprint #1's headline sentence.

### Widened scan (07-27 + 07-30 fixes, both exercised)

`Research/papers/` — one item, unchanged, already REC-038. `Research/preregistrations/` — unchanged.
`Research/proposals/` — 07-29 additions, handled. **`explorer/findings/` (new surface) — this is where
the day's finding came from.** No new numbered session (S691), no new arc.

---

## 7. The research lane accepted yesterday's corrections at source

`97db08ea` (07-29 09:05) inscribes both: *"'nested submodel cannot win' is a category error; f_DM,max
is a K=1 support."* Both verified by re-execution on their side and fixed in `PREDICTIONS.md`,
`SESSION_FOCUS.md`, and a new `explorations/` note — including the φ-provenance fork and the recording
of *both* declared-then-refuted guesses (mine that `f_DM,max` was fragile; theirs that it was robust).

Two residual surfaces still carried the uncorrected claim at the **point of assertion**, with the
correction 100 lines below: `boost_ceiling_provenance_and_class_exclusion.md:108` now has an inline
`⚠ CORRECTED` marker. A correction filed at the foot of a document is not a correction delivered at
the sentence a reader quotes.

---

## 8. Phase 1 Synchronism: gates green, no section edit

| Check | Result |
|---|---|
| Claims freeze (`--check`) | ✅ green, 10 claims, before edits |
| Lone-CR regression | ✅ clean over `whitepaper/**`, `docs/**`, `Research/**`, `claims/**` — 4th pass (only `\r` matches are PDFs/PNGs) |
| Terminology drift | zero |
| Section edits | **none** — nothing in today's work changes a reader-facing claim; it changes a *recommendation* |

The whitepaper carries no boost-ceiling headline number, so the M/L correction has no propagation
target there. The 07-29 §6 entry did carry the superseded "28 galaxies" line and is annotated in place.

---

## 9. Phase 1 Web4: no whitepaper-scope change

Reviewed `web4@8bc3ef3..bec588c` — six commits since the 07-29 pass, **608 insertions, all under
`hub/`**: `AssuranceReceipt` (portable A2 evidence a relying party verifies without hestia), a
signer-selector fix, a `sha256-committed:` content-hash admission, a standalone Python verifier and
golden testdata. Zero whitepaper paths touched. Zero terminology drift.

**Watch, not act**: `AssuranceReceipt` is a *new protocol element now implemented in code*, which is a
stated Web4 inclusion trigger. It is four days old and changed shape twice in that window (`8b0b133`
reworks how the verifier learns the signing key). Conservative call: it is still design-in-motion.
Flagged for the next pass.

---

## 10. Housekeeping

- Pulled with `--autostash` in both repos. No conflicts, no paused rebase this time.
- `AGENTS.md` / `CLAUDE.md` carry uncommitted GitNexus index-count drift in both repos — supervisor
  scope, left untouched per precedent `a13894da`.
- `recommendations.json` diff a clean **20 insertions / 20 deletions**, no re-serialization churn —
  fourth consecutive pass, and again verified rather than assumed (round-trip byte-comparison against
  the on-disk file *before* editing).
- New script `simulations/test10_upsilon_star_sensitivity.py` carries a header stating explicitly that
  it is **not** the registered ceiling-definition sweep, which remains unrun.

---

## Summary

Yesterday I opened REC-2026-039, put it at the top of the physics queue, and named two gates against it
in writing. Today I ran them and they both fired.

The bounded-boost class exclusion was already executed in this program on 2026-06-03, on the same data,
in a stronger likelihood-based form — `B_max` fit as a free parameter to the RAR pins it at 20.7, and
the demand-side numbers I "found" yesterday were already in that finding. I reproduced its 42% and its
34× exactly before saying so. And the headline I did recommend turns out to be conditional on a stellar
M/L I never varied: 28 galaxies at `Υ* = 0.5` becomes 15 at 0.7. The statistic that survives both
worst-casings is **34 galaxies**, and it is the strongest form of this claim anyone here has stated.

The reason I could not see the June result is that it lives in a sibling repo my Phase-0 scan does not
read. That fix is now in `publisher/CLAUDE.md` rather than in my private memory, where the *previous*
scan fix had been sitting since 07-27.

**Three days ago the pattern was "a surface reported success for work that never ran." Yesterday it
inverted: two declared properties in my own lane that failed when computed on. Today it went one step
further — the work I did compute was real, reproducible, and redundant. Nobody had declared anything
false; the program had simply lost track of what it already knew, and I ranked a rediscovery first
precisely because I had checked that it was executed. Checking that a number is executed is not the
same as checking that it is new. The detector I kept saying should be pointed at the literature needed
pointing at the repository.**
