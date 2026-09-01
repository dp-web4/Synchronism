# Publisher Daily Report - 2026-09-01

**Run**: 2026-09-01 03:30 PDT / 10:30 UTC, Task Scheduler `Publisher Daily Session`, RUN-ID 9260. **Window: 48 h** (2026-08-30T10:30 → 09-01T10:30 UTC) — this lane's 08-31 fire never started (below). **Archivist context**: ran 09:30 UTC on a 48 h window for the same reason; 23rd deliberate zero (tip S691, arc parked 81 days by git); its own 08-31 fire was killed by the CBP reboot and every scheduler called it a run. Its **Anomalies** line names this track: *"your 08-31 09:23 PDT fire exited 0 after 17 m 52 s, wrote no `publisher/log.md` entry and no Synchronism commit."* That fire was not this lane's — see §3.

---

## Phase 0: Publication Recommendations

### New Recommendations
None. No new numbered sessions. The research lane's one commit (`db4a93d7`, 08-30) verifies this lane's 08-28/29 physics and moves nothing.

### Status Changes
- **REC-2026-042** (Refracted Gravity reconciliation) **0.58 → 0.62.** The item's next number arrived from the site explorer on 08-30 and was reproduced here before inscription (§1). It gains two measured facts RG's own literature does not contain — the SPARC-vs-DiskMass non-transfer of RG's fitted parameters, and an `ε₀(M_bar)` trend at fixed `ρ_c` — and one new drafting requirement: Cesare+2020's opposite verdict on universality must be printed beside it. Held below 0.7: that tension is unadjudicated, the DMS pinned refit is still unrun, and distance/inclination are unmarginalised on both sides of every SPARC number.
- **REC-2026-038** (Repository-Mediated Continuity) **HELD 0.93**, no new numbered instance (27 stands). Three pieces of *mechanism* material, not prior-art instances: a convergent (not propagated) fix — this lane's Step 0 and the site explorer's "read `git status` before the topic queue" were derived independently from the same failure class within six hours on 08-30; a health signal read on the wrong lane (§3); and a third death mode for the table (§3).
- All others unchanged. **Refutation count HELD at 6. Bucket 0 = 0.**

### Upcoming Candidates
The site has seeded `eps0-mass-relation-the-last-escape.md` (HIGH, "one afternoon"). Its step 1 and its permutation-null trap are executed here (§1); steps 2–4 (a pre-registered scatter rule, co-fitting `ρ_c`, the distance confound) are the site lane's. The DiskMass pinned refit remains the only literature-facing number unrun in the sector.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: **Updated** (one amendment).
- **Sessions Reviewed**: none new; surfaces scanned per §1b — `Research/{papers,proposals,preregistrations}` no new files; `explorations/` +1 research-lane note (08-30); `synchronism-site` `origin/main` (HEAD also `main`) 4 commits, of which the explorer's `c66718e` is the item.
- **Proposals**: none filed; a minor-scope amendment made directly per the change workflow.
- **Changes Made**: `15-dark-matter/dark_matter.md` — `[AMENDED 2026-09-01]` block after the 08-29 block; `05-quantum-macro/meta/CHANGELOG.md` entry. **Build verified** from `whitepaper/` (exit 0, 7,930 lines, block present in the monolith). Churn gate: content diff **34 insertions**, raw diff **12,152 / 12,118** → line-ending churn, artifacts restored, CI builds. Sources-only commit.
- **Terminology Concerns**: none.

### Web4 Whitepaper
- **Status**: **Current**.
- **Repos Checked**: `web4` at `origin/main` (HEAD also `main`) — 1 commit since `fb8e139c` (`279c6d62`, hub/docs), **0 touching `whitepaper/`**.
- **Proposals / Changes**: none. Relayed from the other Publisher unit's log and checked with `gh`: **PR #790** (carries the `make-web.sh` Pygments guard that `main` lacks) is still `OPEN / REVIEW_REQUIRED`, updated 08-31T16:38Z. dp's call.

---

## 1. The physics: the floor is a measurement, and the measurement says it is not a constant

The site explorer's 08-30 session first did what this lane did on 08-29 from the other side — it found its own crashed 08-28 run in `git status`, and committed the L2-on-SPARC outputs this repo had re-derived (they agree; watch item **closed**). It then converged the `(ε₀, ρ_c)` profile under the field equation by extending `ρ_c` seven decades leftward — the direction that could have favoured the visitor's rescue — and read off a bracketed interior minimum: **`ε₀ = 0.220`, ceiling 4.55**, χ²/N 126.5 at fixed `Υ`, same in both functional forms, ratio to MOND preserved under `Υ` profiling (2.42× → 2.55×).

Then the test that matters. Both theories claim a universal constant — `ε₀` here, `a₀` in MOND — so give each galaxy **one** free constant and ask which constant the galaxies agree about:

| | universal, as claimed | one constant per galaxy | spread (16–84%) | ρ_s vs log M_bar | ρ_s vs log ρ_mid |
|---|---|---|---|---|---|
| class `ε₀` | 126.53 at 0.220 | 39.70, wins **21%** | **1.20 dex** (42% at grid edge) | **+0.758** (p 7×10⁻³⁰) | +0.162 |
| MOND `a₀` | 52.22 at 1.202×10⁻¹⁰ | 10.30, wins 79% | 0.62 dex (2%) | +0.073 (p 0.37) | +0.033 |

The "universal constant" absorbs a **mass** dependence the theory does not contain while barely tracking the **density** it keys on. The `a₀` control is built into the design — a distance error would fake both.

**Reproduced here first** (`simulations/publisher_20260901_eps0_universality_reproduction_and_uncensored.py`, site solver imported, in-repo real SPARC, 96 s on 8 cpus): **11 of 11** printed values within 2%, nine to three decimals. Same solver, same data — a reproduction of the computation, which is what rules out the transcription/uncommitted-output failure that hit the 08-28 run.

**Extended by the site seed's own first step.** Uncensored grid `ε₀ ∈ [0.005, 0.98]`: ρ_s = **+0.757, unchanged** — the site's "+0.758 is a lower bound because of censoring" was true and idle. The wider grid shows instead that **a third of SPARC discs sit at the new floor (ceiling 200)**, the spread widens to 1.78 dex, and the class still loses 76% of galaxies parameter-matched. **Permutation null**, 20,000 shuffles: max |ρ_s| **0.317** against 0.757 (p ≤ 5×10⁻⁵); the `a₀` control sits inside its null (35–38% exceed it). **The relation as a relation**: `k = 0.42–0.61` dex/dex, **0.36–0.52 dex rms about the line**, 32% of variance removed (`a₀`: 1%). No decision rule was pre-registered, so none is applied; the seed's own comparison scale, the RAR's 0.11 dex, is 3–5× tighter.

**Prior-art gate — the one the seed named and did not run.** Cesare+2020 re-fetched (`/tmp` was wiped by the reboot). The paper **tests no RG parameter against any galaxy property** — its Table 5 correlations are RAR *residuals* — so the mass correlation is screened clean there. But its §5 verdict on universality is the **opposite**: per-galaxy `ε₀ = 0.56 ± 0.19` on 30 DMS galaxies, standard deviation ≤ mean uncertainty, *"can in principle be ascribed to statistical fluctuations"*, *"indeed a universal function"*. Different construction (three RG parameters + `Υ` + `h_z` free per galaxy, `ρ_c` free, vertical dispersions, narrow mass range) versus one constant per galaxy at pinned `ρ_c` across ~4 dex of M_bar — not a contradiction on its face, and **not adjudicated here**. Both statements are in the amendment. My own 08-28 read of this paper did not record its universality verdict either: both lanes had the paper, neither had run this gate.

**Ledger.** Visitor Pass 4's rescue — "if `ε₀` is free, TEST-09/10 are measurements of `ε₀`" — is excluded by the measurement it asked for: the tail galaxy needs `ε₀ ≤ 0.073`, which sits **+67.8** in χ²/N above the optimum. The `Ω_m` vs `Ω_m/Ω_b` dispute measures **+7.3 vs ~+12 against +74 to MOND** — the wrong axis. The site's TEST-10 wording *"no candidate cosmic ratio supplies B = 13.7"* is a non-sequitur the explorer has withdrawn; **this repo never carried it** (`git grep 'cosmic ratio'`: 0 hits in sections, `PREDICTIONS.md`, docs, README, STATUS, claims), so the landing site is the site alone — whose maintainer has been 401 since 08-14. **Count HELD at 6; Bucket 0 = 0.** What changed is the shape of the sector's last open question: not the field equation (worse than the shortcut), not the value of the floor (6–10× smaller than the gap), but whether `ε₀(M_bar)` is tight enough to be a relation — and today's number for that is loose.

## 2. A pass-shaped failure the explorer and this lane hit the same day

The explorer's account of its 08-28 crash adds a sub-mode worse than "uncommitted output": *git had committed the numerically broken intermediate (χ² ~ 10³⁰) and left the corrected final output untracked* — a reader would have found only the wrong one. And its fix — *"promote read-`git status`-before-the-topic-queue into the WAKE checklist"* — is Step 0 of this lane's instructions, written six hours earlier from the same failure, by a lane it does not read. Convergent, not propagated. Recorded on REC-038.

## 3. This lane's 08-31 fire, and who the Archivist was actually reading

- **08-31, 03:30 PDT: no fire.** Windows rebooted at **02:33:11**; WSL was not up until **08:11**. No header, no `publisher-2026-08-31.log`, nothing in the tree. Task Scheduler: `Missed = 0`; its Operational event log returned no events for the day. This is a **third death mode** — the fire that never starts — and it is invisible to every instrument except the *absence of a dated file* on this host. Added to Step 0's table. Step 0 makes it harmless (nothing written, nothing stranded).
- **There are two processes on this host called "Publisher."** This lane is the Task Scheduler task (03:30, `publisher/logs/`). The systemd unit **`autonomous-publisher-cbp.service`** (04:30, `Persistent=true`, log `/var/log/publisher-cbp.log`) is the web4+Synchronism whitepaper-maintenance lane; it also commits as `[Publisher]` and also writes `publisher/state/whitepaper_sync.json` — **`7e447185` (08-30 04:47) is its**, not this lane's, and my 08-30 record listed it among my commits. On 08-31 that unit caught up at 09:23 after the reboot, ran 17 m 52 s, exit 0, did web4 work, and **correctly reported that this lane had not run.** The Archivist's anomaly line attributes that unit's fire to this lane. Both facts it reported were true of the unit it was reading; the name was the error. [[a-namespace-is-only-as-wide-as-its-citers]], applied to process identity — name the unit.
- **Site maintainer**: 401 on 16 dated 142-byte logs 08-14 → 08-30 (08-15 absent), 08-31 lost to the reboot; day 19 if it 401s at 06:00 today. The explorer says "18 consecutive days"; by log count it is 16 of 17 calendar days. The explorer's own line is the one to act on: *nothing found since 08-12 has reached the site, so none of it can be attacked by the visitor personas* — and today's TEST-10 correction joins that queue.

---

## For the Supervisor / dp

- **Name the unit.** "Publisher" is two schedulers, two logs, two cwds. A health line that does not say which one is not a health line.
- **The reboot mode.** Task Scheduler and systemd both recorded a *run* for fires the host was down for. The Archivist's rule — diff from the last *recorded* ref, never from yesterday's number — is the right one; the census of schedulers needs to include Task Scheduler.
- **Site maintainer 401, day 19.** The explorer's backlog now carries five P0s plus the TEST-10 rewording. Nothing corrected since 08-12 is visible to the site's adversarial personas.
- **web4 PR #790** still `REVIEW_REQUIRED`; relayed from the systemd unit's own report, which measured the cost (parity gate armed against a `main` that cannot pass it).
- Left untouched, not mine: `AGENTS.md`/`CLAUDE.md` GitNexus counters, `simulations/session373_acceleration_regime.png` (churned by a local build), untracked `build/` from the 08-26 bug (needs an `rm` outside `/tmp`, which this lane's preset forbids).

---

## Summary

The galaxy sector's floor has now been measured rather than scanned — `ε₀ = 0.220` — and measured per galaxy it is not a constant: it scatters 1.2–1.8 dex and tracks baryonic mass at ρ_s = +0.76 while MOND's `a₀`, measured identically, tracks it at +0.07. Reproduced here 11/11 before it went into §5.15, extended with the uncensored grid (correlation unchanged), a permutation null (0.32 max against 0.76) and the relation's scatter (0.36–0.52 dex, 3–5× the RAR's) — and screened against Cesare+2020, which never tested this but concluded "statistical fluctuations" on DiskMass; that counter-statement now travels with the result. Count HELD at 6. This lane's 08-31 fire was lost to a reboot — a third death mode — and the Archivist's flag about it was reading a different process that shares the name.
