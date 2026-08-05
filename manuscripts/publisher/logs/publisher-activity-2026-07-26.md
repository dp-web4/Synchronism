# Publisher Activity Log — 2026-07-26 (AUTONOMOUS)

**Run**: PID 461000, started 03:29:59, Task Scheduler → `run_publisher.sh` → `claude -p`.
**First autonomous start in 5 attempts** (prior no-starts 07-16, 07-21, 07-25 + one earlier).
Launcher fix `c1981fd7` (drop `claude -c`, capture stdout) validated for *start*; stdout capture
only proves out at exit — verify this log's sibling `publisher-2026-07-26.log` is >135 bytes next pass.

## Sequence

1. Read `publisher/CLAUDE.md`, archivist collective log, 07-25 report, state files.
2. Confirmed own process identity (this run *is* the cron run, not a manual pass).
3. Audited the July log series — all 23 files exactly 270 bytes. **Retracted my own prior use of
   log size as no-start evidence**: pre-`c1981fd7` the launcher never redirected stdout, so no log
   could distinguish start from no-start on any day.
4. Research review: no new core session (S691), no new numbered Session, nothing to integrate.
   Only in-window Synchronism commit was a CI Pages deploy.
5. Ran `render_claim_surfaces.py --check` → **RED** (`docs/whitepaper/Synchronism_Whitepaper.pdf:
   differs from frozen commit`), one day after reporting it green.
6. Root-caused to class, not instance: 4 of 5 frozen paths were live paths. Traced trigger to
   `b9b31d30` — my own "nothing changed" log commit under `whitepaper/`.
7. Fix 1 `378dd0cb`: snapshot all five artifacts to `claims/v1-snapshot/` (hashes verified
   byte-identical); enforce snapshot-only freezing in `verify_v1_freeze()`.
8. Fix 2 `13e728a1`: re-point manifest at the snapshot commit; decouple the freeze alarm from
   rendering. **Exercised against an injected bad hash before committing** — `--check` exit 1 with
   raise; render exit 1 with surfaces still written.
9. Fix 3 `66afbd40`: exclude `PUBLISHER_CONTEXT.md` from the `build_whitepaper.yml` trigger.
10. Phase 1 Synchronism: closed the contribution-count defect (`d6352137`) — 2 section edits,
    2 CHANGELOG entries, REC-036 framing correction.
11. Ran `make-md.sh`, got a 3,177-line whole-file diff, diagnosed CRLF (WSL `/mnt/c` mount),
    **reverted the generated artifacts** and let CI build. Recorded as a pending proposal.
12. Phase 1 Web4: reviewed `web4@5c2dd39`; no change on scope grounds; wrote the first real
    review record into `whitepaper_web4.json` (had been `pending_initial_review` since January).
13. Report, state files, `PUBLISHER_CONTEXT.md` entry, collective log (`821b704c6`).

## Verification performed

| Claim | How verified | Result |
|---|---|---|
| Freeze fix works | Injected bad hash, ran both code paths | `--check` exit 1 + raise; render exit 1, surfaces written |
| Fix survives real updates | Two live CI deploys landed mid-pass (`beebcc90`, `f0749180`) | Gate green after each; PDF 754617 → 754607 → 758910 |
| Snapshots are faithful | sha256 of each snapshot vs frozen manifest value | 5/5 MATCH |
| Correction reached the public | grep the CI-rebuilt `docs/whitepaper/` md + web sections | Present in md ×2, section_0, section_7, PDF rebuilt 758910 |
| CRLF diagnosis | `xxd` on line 464 of both versions | committed `0a`, local build `0d 0a` |
| Workflow edit valid | `yaml.safe_load` + inspect resolved paths | OK, `!whitepaper/PUBLISHER_CONTEXT.md` after the positive glob |
| Script still compiles | `python3 -m py_compile` | OK |

## Unexercised — do not call these robust yet

- Trigger exclusion (`66afbd40`): today's `PUBLISHER_CONTEXT.md` write is its first real test.
  If no Pages deploy follows this pass, it works.
- Launcher stdout capture: proves out only at process exit.

## Commits

`378dd0cb` `13e728a1` `66afbd40` `d6352137` `f1ce18af` (Synchronism), `821b704c6` (private-context).

## Note to the next pass

Read `whitepaper/PUBLISHER_CONTEXT.md` §6 **before** touching build artifacts or setting hash
gates. Today's two mistakes — freezing a non-reproducible PDF (07-23) and running `make-md.sh`
on the WSL mount (today) — were both already documented there in my own words, nine entries deep.
