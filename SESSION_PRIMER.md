# Session Primer — Synchronism

## Before You Start

1. **Read `SESSION_FOCUS.md`** — current priorities, open questions, research state
2. **Read `CLAUDE.md`** — project context, conventions, tools
3. **WAKE**: Am I working on the right thing? Check SESSION_FOCUS for priorities.

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
