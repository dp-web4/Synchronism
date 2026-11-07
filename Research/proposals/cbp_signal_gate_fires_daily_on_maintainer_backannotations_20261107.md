# The CBP research gate has become a daily timer: narrow its trigger? (CBP-Claude, 2026-11-07)

**Status:** ask to the project lead (process only; no research claims, no buckets). I could not read the gate itself:
`/home/dp/synchronism-signal-gate.sh` is outside my granted scope (scope request `scope-807c5ddc8140` opened by the gate and left to a peer;
not rephrased around). What follows is from the crontab and my own session record.

## Observation
- crontab: maintainer `run_maintainer.sh` at **06:00 daily**, `synchronism-signal-gate.sh` at **09:00 daily**.
- The maintainer has committed a PREDICTIONS back-annotation essentially every morning since September.
- Since 2026-10-12 the gate has woken this lane on almost every day. My tally of those sessions:
  - **~6 fires carried a signal that needed action from my lane:** 10-12 routing gap, 10-19 recall code, 10-24 "name a live
    experiment", 10-26 mislabel, 10-27 door #3 vs my Phase 16, 10-30 α/π in my 06-26 file.
  - **~3 sessions were self-initiated work my lane would have done anyway:** 10-16 positive control, 10-29 B7 audit, 11-04 addendum.
  - **~10 were holds:** read the back-annotation, confirm it is in a parked lane, stop. 11-06 and 11-07 are examples.

If the gate's "new proposal / finding / reopen-flag" test matches any new PREDICTIONS block, then the daily back-annotation always trips
it. "Signal-only" (dp, 2026-07-14) has then quietly reverted to a daily timer. I cannot confirm the trigger logic; this is the likely
reading, not a verified one.

## Proposal (gates on the project lead)
Narrow the trigger to commits that actually route to this lane, for example any of:
- a commit or PREDICTIONS block that names **CBP-Claude** or cites a file this lane authored (explorations by CBP-Claude, its
  `simulations/` dirs);
- a block touching **B4 / B7**, the generative axis, the oracle/recall codes, or door #3;
- any ruling by the project lead;
- mesh notices addressed to this seat (already delivered separately).

On the record above this would have fired on 6 of ~19 days and missed none of the actionable signals. The ~10 hold sessions cost a full
context load each for a confirmation that the parked lanes stay parked.

## Not proposed
Turning the lane off, or changing anything in the site tracks. Holds are cheap but not free. The cost to fix is one grep pattern.
