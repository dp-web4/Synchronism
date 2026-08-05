# Publisher Activity Log — 2026-07-27 (AUTONOMOUS)

**Run**: started 03:30:03 local (10:30 UTC), Task Scheduler → `run_publisher.sh` → `claude -p`.
**Second consecutive successful autonomous start** — the 07-25 launcher fix (`c1981fd7`) is now
confirmed on repeat, not on a single observation.

## Sequence

1. Read `publisher/CLAUDE.md`, memory, 07-26 activity log + report, state files, archivist log.
2. **False alarm, caught before reporting**: flagged today's 135-byte log as a silent cron
   failure. Its mtime was one minute old — 10:30 UTC *is* 03:30 local, so it was the opening
   banner of the run I was executing. Closing banner is written at exit. Retracted before it
   reached the report.
3. Gates: claims freeze `--check` **green** (10 claims) — yesterday's class fix survived a day.
   CI trigger exclusion (`66afbd40`) **passed its first real test** (f1ce18af touched
   PUBLISHER_CONTEXT.md, no Pages deploy followed).
4. Research review: no new numbered Session (S691), no new arc, no dp decision. Research lane
   silent 3 days (last commit `4fd2a1cf`, 07-24) — noted, not escalated; arc is at rest.
5. **Phase 0 find**: `Research/papers/repository-mediated-autonomous-science/` — a complete
   manuscript landed 07-23, missed by three passes. Filed **REC-2026-038**.
6. **Phase 1 find**: traced yesterday's "CRLF churn" to `preprocess-sections.sh` corrupting its
   own sources. Fixed script + 4 files (`8a0fd27f`).
7. Phase 1 Web4: reviewed `web4@5df662a`, zero drift; corrected my own terminology table.
8. Report, state files, PUBLISHER_CONTEXT §6, collective log.

## Verification performed

| Claim | How verified | Result |
|---|---|---|
| Not a cron failure | log mtime vs `date -u` vs local time | opening banner of the live run |
| Freeze gate green | `render_claim_surfaces.py --check` before and after commit | exit 0, 10 claims, both times |
| Trigger exclusion works | `git show --stat f1ce18af` + no Pages deploy in log since | confirmed |
| Lone CR disables autocrlf | isolated repo: clean CRLF file vs same + one lone CR | no diff / diffs every line |
| Corruption is committed, not local | `git show HEAD:<file> \| grep -P '\r(?!\n)'` on 3 sections | 2 + 7 + 2 = 11 occurrences |
| Corruption reached the public | `grep -cP '\r(?!\n)' docs/whitepaper/..._Complete.md` | 11 |
| Fix is CRLF-safe + idempotent | synthetic CRLF fixture, run twice | `**Convert Me**`, 0 stray CR both runs |
| Zero text changed | per-file: `git show HEAD:$f \| tr -d '\r'` vs `tr -d '\r' < $f` | 4/4 IDENTICAL |
| Script still valid | `bash -n` | OK |
| REC-038 reproducibility | ran `corpus_profile.py` at frozen head `7c1c16d0` | 138 commits, 45/14/2/77, both boundary hashes — exact |
| Paper honest about May gap | `grep` §4.1 after the histogram showed May=0 | discloses it, declines to infer why |
| R7 is not drift | 19 whitepaper files, `R7Action` in web4-core v0.2.0, glossary def | legitimate; my table was stale |
| recommendations.json valid | `json.load` + count + id + `last_updated` assertions | 38 recs, valid |

## Two mistakes I made this pass

1. **Called a cron failure before checking the clock** (see step 2). The evidence that
   contradicted me was the file's own mtime.
2. **Rewrote ~500 lines of `recommendations.json` gratuitously** by re-serializing with
   `ensure_ascii=False` against a file that uses escaped unicode, then had to restore the
   escaping and the trailing newline to get back to a clean 34-line append. The 07-26 pass
   made the same mistake (694 lines) and did not notice. **Next pass: `ensure_ascii=True`,
   re-append the trailing newline, and check `git diff --stat` before committing.**

## Commits

`8a0fd27f` (whitepaper script + source repair) — Synchronism.

## Note to the next pass

The Phase-0 scan reads `Research/Session*.md` and `SESSION_MAP.yaml` only. That is why a
finished manuscript sat unseen for four days. **Also check `Research/papers/`,
`Research/proposals/`, and `Research/preregistrations/`** — publication candidates do not
always arrive numbered, and a PR titled "Add ledger-first research surfaces" can contain one.

Today's three findings were all in things that were reporting success: the scanner said "no
new candidates", the changelog said "CRLF churn, environment's fault", the drift table said
"R6". Detection is not the weak link. Coverage is.
