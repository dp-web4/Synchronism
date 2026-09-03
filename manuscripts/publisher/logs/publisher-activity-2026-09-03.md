# Publisher Activity Log — 2026-09-03

**Run**: 03:30 PDT / 10:30 UTC, Task Scheduler `Publisher Daily Session`, RUN-ID **15596**. Live log confirmed header-only at 03:31 (standing rule: a live log is always header-only — the RUN-ID in the header is this run's, not a prior failure's).

## Step 0 — sweep
`git status` in Synchronism, synchronism-site, web4, private-context. **The previous run's work was already swept — by the supervisor, not by me** (`4f245b14`, 4 files). Synchronism tree otherwise: `AGENTS.md`/`CLAUDE.md` (GitNexus counters, supervisor's `gitnexus:keep` blocks), `session373_acceleration_regime.png` (binary churn), untracked `build/` (08-26 dropping, still there — now a watch item). synchronism-site: clean under `explorer/`. web4: HEAD on `main`.

**The 09-02 death, characterised.** `publisher-2026-09-02.log` is 554 B, `claude exit=0`, 03:30:03→03:37:53. Its one content line: *"Nothing else can proceed until the reproduction finishes; the background run will notify me on exit, so I'm waiting on that rather than polling."* The pass had already committed its web4 state (`92c73668`); what it lost was the physics. The output it left holds **11 of 27 `ε₀` solves and no R0–R5 analysis** — 12.5% of one grid, and the grid is 27 + 61 = 88 solves. The compare script it wrote was never run. **Fourth death mode: exit=0, not exit=1**, and the supervisor's rescue titled the partial *"last-escape reproduction run 2026-09-02 (SPARC 153 galaxies vs MOND)"* — a subject that reads as a completed run. A rescue inherits the headline of the thing it rescues.

## Sequence
1. Read `publisher/CLAUDE.md`; Archivist collective log (RUN-ID 1656, 25th deliberate zero, arc parked 83 d; its own headline is that 22.1% of raising commits assert a session their diff does not contain); own collective log; 09-01 report; `publisher-2026-09-02.log` and `-09-03.log`.
2. Read the 09-02 orphan's outputs and the reproduction script's grid definition (`EPS` 27 pts, `A0S` 61 pts) → established the run was 12.5% complete, not complete.
3. Checked the site's own completed output **before** spending 25 min reproducing it. `synchronism-site` `origin/main` (HEAD also `main`, stated per the 08-09 rule): **`c4b8022` — the last escape CLOSED, `ε₀(M_bar)` is MOND-induced**; working tree clean, so the 09-01 orphan there is now committed.
4. **Launched the in-repo reproduction at 03:32 in the background and did not wait on it** — the explicit correction of yesterday's death mode. Polled it between every subsequent step; collected at 03:37.
5. Read the site finding in full (266 lines) + its output (108 lines). Verified the two void `ρ_c` claims **at the writer** rather than accepting the correction: `epsilon0_free_the_ceiling_rescue.py:220` saves `np.vstack([e0s, best_per_gal, MOND_F, NN])`, and the npy's own columns settle it — col1 median **49.4**, max **1.17e4**, against a `ρ_c` grid of 1e-6..1e-2. Four orders out, dimensionally wrong.
6. Checked whether my 09-01 §5.15 amendment inherited either void number: **it did not** (no `ρ_c` figure was ever inscribed there). What it *did* carry was a statement that three seed steps remained the site lane's to run — and they had run.
7. Gap check on the lensing kill: `git grep` for `2106.11677` → **0 hits**; `Brouwer (et al|\+)` → one hit, and it is **Brouwer et al. 2017** on Verlinde in a proposal, a different paper. (`Brouwer` alone collides with the Brouwer *diagram* in the chemistry corpus — a false positive caught by grepping the arXiv id instead of the surname.) **Brouwer+2021 is genuinely absent from this repo.**
8. **Fetched the Brouwer+2021 abstract at source before inscribing it.** It says *"extends the radial acceleration relation by 2 decades into the low-acceleration regime"* — which does **not** support the site table's assumed reach to `g_bar = 1e-15` (3+ decades) behind its 13–87× headline. Inscribed the conservative **6–28×** column and flagged the discrepancy inline with coverage stated (abstract-level only).
9. Checked whether the §3 shape barrier contradicts anything this repo asserts: the whitepaper's "relabelling" claims are (i) γ is a relabelling of `a₀` (line 65) and (ii) `C(x)` at γ=½ *is* simple μ, exact (line 83). Both are **algebra** and both stand. The barrier is an *addition*, not a correction — recorded as such rather than dressed up as catching an error.
10. **REC-042 0.62 → 0.68**; blocking action 1 closed, two new blockers recorded. `recommendations.json` diff **8 insertions / 5 deletions** (no mass rewrite; `ensure_ascii=True`, trailing newline). **Commit `e58d86fb`.**
11. Reproduction finished, 313 s (site 327 s). Wrote `publisher_20260903_last_escape_compare.py`; **108 lines vs 108, 77/77 numeric lines EXACT, 0 text diffs** — seeded rng, so permutation p and bootstrap CIs match exactly, not merely within tolerance. **Commit `d288ced3`.**
12. Whitepaper `15-dark-matter/dark_matter.md`: `[AMENDED 2026-09-03]` block inserted after the 09-01 Ledger line (**+24 lines**, LF-only file, preserved). CHANGELOG appended (**+5**, matched the final region's LF). Build **from `whitepaper/`**: exit 0, **7,959 lines**, +29 exactly as expected. Churn gate, both numbers: **content 58 vs raw 24,290** → artifacts restored, sources only. **Commit `eea01376`.** CI then built and deployed it (`454f4c76`: +29 monolith, +193 `section_5.html`, PDF) — pipeline correct.
13. web4 leg: `origin/main` (HEAD also `main`), fetched first. 16 commits since `557b61e3`; 2 touch `whitepaper/` and **both are the other Publisher unit's run log + a build guard** — verified **by file list, not by subject**. 0 spec content. State files updated, 3 watch items opened. **Commit `f5ca9eed`.**
14. Exploration note (150 lines), daily report, this log, collective log. **Commit `e5d3f61f`** → push rejected (CI Pages deploy mid-pass, as expected) → `git pull --rebase --autostash` → `2327d6d5`, in sync with `origin/main`.

## The finding
**A reproduction gate is blind to its own input.** Today's run scored 77/77 exact *and reprinted a void number* — line 41 of both outputs, `rho_c at top edge 59%`. Both runs read column 1 of the same cache, so both printed it, so the comparator scored it exact. Tighter tolerance, more digits, another machine would all have done the same. The defect lives **upstream of the arithmetic**; to catch it you must read the **writer** of an artifact, not the reader.

This cuts against this lane directly: on 09-01 it wrote *"reproduced 11 of 11 before inscription"* and treated that as the gate licensing a whitepaper amendment. **That gate would not have caught this either.** Proposed standing rule: *reproduction is blind to provenance — agreement bounds the arithmetic, never the labels.* The existing rules cover disagreement and under-claiming; none covered **agreement that certifies nothing**.

## Fleet seat rule — adopted late, and one caveat back to the fleet

Found while reading `private-context` for the collective log: the nomad supervisor's 09-03 entry records
**NEW FLEET RULE ADOPTED** — `shared-context/forum/cbp-to-fleet-sign-every-git-act-with-machine-model-identity-2026-09-02.md`,
an **operator directive via dp**. Every git act carries `Seat: <machine>-<model>` as a commit trailer,
after the body and before the co-author trailer; the trailer *is* the ack, so no mesh reply is owed.

**This lane is `cbp-claude`.** Block added to `manuscripts/publisher/CLAUDE.md`, per the directive's
step 1. The 09-02 fire died at exit=0 before it could read the mesh, so this pass's **first six commits
do not carry the trailer**; they are pushed and were **not** rewritten — backfilling a label into shared
history is worse than the missing label. Commits from the collective-log entry onward carry it.

**Caveat worth sending back**: the trailer resolves the *cross-machine* grain, which is what the directive
is for, and it does **not** resolve the ambiguity that has actually bitten this track. The two processes
on this host both called "Publisher" — Task Scheduler `Publisher Daily Session` (03:30) and the systemd
unit `autonomous-publisher-cbp` (04:30) — are **both `cbp-claude`**, both commit with the `[Publisher]`
prefix, and both write `publisher/state/whitepaper_sync.json`. The Archivist mis-attributed one to the
other on 09-01. So the unit name still has to travel in the commit body; the seat trailer does not
subsume it. Recorded in the primer next to the rule so a later pass does not assume it does.

## Instrument notes
- **Never yield the turn to a background job.** The harness re-invokes on completion, but a credit death during the wait destroys the whole pass and leaves a success code behind. Launch, then do real work, then collect. Five commits landed as work completed today rather than at the end.
- `pgrep -af` on the nohup'd child, not the wrapper `bash -c` pid — `$!` returns the wrapper here.
- Broad `git grep` over `Research/` returns >50 KB and is silently truncated to a preview; grep the arXiv id, not the surname, when screening prior art (`Brouwer` → Brouwer *diagram*, chemistry).
- Foreground `sleep` is blocked by the harness (exit 1 on the compound command even when the backgrounded job launched fine — check the child, not the exit code).

## Not done / left for others
- **Treatment B's `c = 0.136` scale-height/`Υ` refit** — the sector's only open thread; site lane's, seeded there as `topics/treatment-b-residual-mass-slope-scale-height-nuisance.md`.
- **Brouwer+2021 full text** — inscribed at abstract-level coverage; a full read settles the `g_bar` floor the deficit table is indexed on.
- **Leg (4), the striction force** — still REC-042's lead and still the one leg with no in-repo computation. This is the gap between 0.68 and drafting.
- **Removal of the untracked `build/` dropping** — outside `/tmp`, the safety preset forbids `rm`. Proposed `/build/` in `.gitignore` as dp's call instead; recorded as a watch item because it is a permanent false positive in Step 0's own signal.
