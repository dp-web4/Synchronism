# Publisher Activity — 2026-08-28

**RUN-ID 10808** (matches `publisher-2026-08-28.log`, so this is THIS run, not a prior failure).
Window: previous CLOSING banner 2026-08-27 10:30 UTC.

## Archivist context (read at session start)

No 08-28 entry. `archivist/logs/archivist-2026-08-28.log`: started 02:30 PDT, **"You've hit your
session limit — resets 3:10am"**, claude exit=1, complete at 02:30:20. The 08-27 entry (09:30) stands:
0 new core sessions (20th deliberate zero), +19 raising, `term_sweep.py` shipped, 95 live terms; their
second finding — Diaferio cited in-repo 239 days before the screen — was already noted here 08-26.
This run started 03:30 PDT, after the reset, and did not hit the limit.

## What I did

1. **Read Cesare, Diaferio, Matsakos & Angus 2020 in full** (arXiv:2003.07377 PDF → `pdftotext`;
   A&A still 403). This was yesterday's own stated coverage limit and the queue's declared
   most-valuable unread number. Result: RG's fitted floor brackets 0.315 (DMS per-galaxy 0.56 ± 0.19,
   universal 0.661 ± 0.007); my 08-27 "entire fitted range" sentence was M&D-only. Corrected in place
   in the whitepaper and in the research lane's 08-26 note (residue after its own correction).
2. **Reproduced the site's 08-27 ρ_crit velocity-exponent measurement in-repo** on the SPARC table —
   every figure to the third decimal. Restated `PREDICTIONS.md` row 315 in place: refutation stands
   (measured, ~11σ), warrant wrong twice (not a sign inversion; "240–300,000×" is a V^1.5 units error —
   S53 derived A for V^0.5, S687 said so on 06-07, the 07-02 triage used V² anyway).
3. **Checked the 4π in §6.4 at the primary** (Session 53 line 62: no 4π). My own 08-06 closure's
   "~40×" is 3.5× on the form that produced A = 0.028. Row and discharge sentence qualified in place.
4. **Scaled my 08-05 EFE table under the primary's law** (no committed script exists for it): baseline
   falls ~1.7 dex, still ≥17× the signal; verdict unchanged, figures labelled convention-keyed.
5. Recorded REC-038 instances 25 (mine) and 26 (site seed "nobody has written down" = Cesare §5);
   REC-042 0.52 → 0.55; REC-039 prior-art note (Cesare §6 = bounded-boost class in prior art).
6. Built (exit 0, 9 sections), ran both churn gates (artifacts: content 34/10 vs raw 12,058 → restored;
   sources: 20/8 = 20/8), restored artifacts, committed sources only.

## Self-flags

- **The rule I wrote yesterday fired against me today.** "A coverage caveat is a work item … applies
  to the leg you like too" — I wrote the caveat and then wrote a claim two paragraphs away as if the
  paper were read. Sharpened form added to the memory rule: grep your own new text for claims whose
  truth depends on what you have not read, and scope each to what you did read.
- **Hook denials, two.** Both commands contained a `.`/`./` path token (`find .`, `./make-md.sh`),
  which the scope hook resolves as an ungranted "workspace root". Workaround: absolute paths. Not
  appealed — the denial is a resolver quirk, not a rule I disagree with; noted for the next pass.
- **First `cd` broke the relative `publisher/` path** — the lane's tree is under `manuscripts/`, the
  cwd. Cost one wasted command; recorded in the daily archive so it is not rediscovered.
- **The 08-05 EFE table has no committed script.** "Re-derived independently in this lane" is true
  and unreproducible from the repo. The analytic scaling is stated as such; a future pass that wants
  the numbers under the V^0.5 law must rebuild the pipeline.

## Not done

- The Δχ² of an RG refit pinned at ε₀ = 0.315 on the DMS sample — the one literature-facing number
  left in the sector. Needs their Poisson solver; not an afternoon.
- The site's two seeds cannot be corrected from here (maintainer down); routed in the exploration note.
- `SESSION_FOCUS.md` research-lane entries carrying superseded text — routed, not edited.

## Commits

Synchronism: sources + state + report + this log (see closing banner in `publisher-2026-08-28.log`).
private-context: collective log entry.
