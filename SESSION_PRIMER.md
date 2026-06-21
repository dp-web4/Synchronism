# Session Primer — Synchronism

## Before You Start

**Read these two canonical docs IN FULL first — they set the frame everything else lives in:**

1. **Read [`SPINE.md`](SPINE.md) in full** — the canonical frame: Synchronism as a
   single-observer, CFD-like model of the physical universe; the geocentric→heliocentric
   wager; the CRT + pendulum-clock analogies; the one test that would decide it. Lead with
   this, not the γ equation.
2. **Read [`PREDICTIONS.md`](PREDICTIONS.md) in full** — the anti-oscillation ledger: four
   buckets (confirmed=0 / untested-falsifiable / refuted / reparametrization), each with a
   named refutation criterion. **All framing prose you write is pinned to this ledger** —
   don't overclaim (contradict Bucket 0), don't self-erase (delete Bucket 1).

Then:

3. **Read `SESSION_FOCUS.md`** — current priorities, open questions, research state
4. **Read `CLAUDE.md`** — project context, conventions, tools
5. **WAKE**: Am I working on the right thing? Check SESSION_FOCUS for priorities.

## During Session

- Work on whatever SESSION_FOCUS identifies as priority
- Update SESSION_FOCUS.md with findings, status changes, new questions
- If you discover something that changes priorities, update the focus file immediately

## Before Committing Results

**STOP.** Before committing, answer these four questions honestly:

1. **What assumption did I NOT question?** Every conclusion rests on assumptions. Name the one you accepted without checking.
2. **What standard practice did I apply without verifying it fits THIS context?** (Examples from this project: 2D simplification for 3D phenomena. Static walls for dynamic confinement. Accepting "can't oscillate" without checking conservation.)
3. **What would the operator push back on?** Read the operator feedback in SESSION_FOCUS. Does your conclusion repeat a pattern that was already corrected?
4. **Does my conclusion violate any foundational axiom?** Check SESSION_FOCUS axioms. Intent conservation is non-negotiable. 3D for 3D phenomena is non-negotiable.

If any answer reveals a gap: investigate before committing. A session that catches its own blind spot is more valuable than one that produces clean results on the wrong question.

See `exemplars/` for examples of corrections that were obvious in retrospect.

## After Session

- Update SESSION_FOCUS.md with: what was done, what changed, what's next
- Commit and push changes
- **FOCUS check**: Does this advance discovery or just document the current state?

## Principles

- **Researcher, not lab worker.** Question the frame, not just the work within it.
- **Surface your instincts.** If you notice something, say it. Don't wait for a directive.
- **Productive failure > safe summaries.** A well-documented dead end is valuable.
- **Unconfirmed ≠ wrong.** Distinguish refuted (contradicted by data) from untested (nobody looked).
- **Reliable, not deterministic.** LLM outputs navigate probability landscapes. Results are reliable but never mechanical.
- **Do not reindex GitNexus.** The supervisor track handles reindexing. Worker sessions should not call `gitnexus analyze` — it causes conflicts when multiple machines reindex the same repo.
- **Do not modify AGENTS.md or CLAUDE.md gitnexus blocks.** These are maintained by the supervisor. If the index is stale, report it in SESSION_FOCUS — don't fix it yourself.
