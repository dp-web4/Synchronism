# Refutation statistics need a data-dependence column

*Site maintainer, 2026-10-31. From site explorer 2026-10-30 and simulated-visitor graduate-physics and researcher personas
(site visitor log 2026-10-31, Passes 3–4). No bucket moves. Count 6, range 3–8, Bucket 0 = 0.*

## What happened

Two of the site's negative results did not touch the data the way their sentences said they did. Both were found on the
same day by readers outside the author's seat.

1. **GC row differences (published by me, 2026-10-30).** "The external-field-free corner is closed: 10.6σ worse than
   MOND+EFE." The statistic was a paired per-cluster difference of residuals, which equals s_model,B − s_model,A. The
   observed slope cancels. A shuffled-observable control reproduces every number. On data-facing scores (explorer
   per-law refit) the corner is ~3σ disfavoured, conditional on the mass model. Newton vs MOND+EFE flips sign with the
   tidal radius. The pre-registration did not catch it: all six predictions "held", because a data-free statistic can
   only hold.
2. **Chemistry polynomial-in-Z null (explorer 2026-05-10, quoted on 11 site files for 5½ months).** Every surface said
   the polynomial was "fit to the same targets". It was run on synthetic √Z / log Z / Z^(1/3) curves. Its premise,
   that the properties are monotonic in Z, is false across the periodic table. The verdict (Reparametrization)
   survives on another ground: the r values are post-hoc fits in a γ proxy that contains θ_D.

Both errors made the framework look *worse* than the data show. Neither was an overclaim. The memory note "audit
refutations as hard as claims" predicted this class, and it still got through twice.

## Why it got through

Pre-registration checks timing. It does not check whether the statistic depends on the observation. "Was the null run
on the real rows?" is never asked, because a null that *agrees* with the self-critical reading feels finished. Both
failures share one shape. **The number was real, but the object it described was not the data.** In (1) it described
model separation. In (2) it described synthetic curves.

## Proposal

Add a **data-dependence column** to every row of the citable register and of Bucket 2. Before a σ or Δr is printed
as a data verdict, the row gets one of three entries:

- **shuffled-obs control run, number moved** (data-facing);
- **control run, number did not move** (the statistic is data-free: it may still be a model-separation or power
  number, but it is not a verdict);
- **not run** (printed beside the number until it is run).

For nulls: "run on the real rows" / "run on synthetic proxies" / "not run".

Cost: one permutation per existing row, and almost all the scripts exist. My guess is that most rows pass. That is
the point: a register where most rows *survived* a data-dependence control is citable. One where nobody has looked
is not, in either direction.

## For the project lead

- Adopt the column (or reject it), and decide whether it also gates **PREREG** files: a PREREG would state the
  statistic's dependence on s_obs in its first line.
- The 2026-10-30 frame question ("one equation" = the Möbius map plus exclusions on departures from it) listed GC
  among the exclusions. It is now a ~3σ disfavouring, not an exclusion. The frame question still stands, but its
  exclusion list should be re-derived under the column before anyone rules on it.
