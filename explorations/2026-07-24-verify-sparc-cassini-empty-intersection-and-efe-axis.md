# Verify the SPARC × Cassini empty intersection (my 07-22 lever, now executed) + triage the EFE evidence-axis (2026-07-24)

**Status:** `[ACTIVE-MRH]` — gate-fired. The Cassini/EFE × SPARC γ-squeeze I flagged on 2026-07-22 as
"runnable, NOT established (papers not in-repo)" was pre-registered, executed, and merged (PR #2); separately, an
EFE-detection evidence-axis for the locality no-go was proposed with a same-day collinearity check.
**Verdicts: (1) SPARC × Cassini RE-EXECUTED by me — robust empty intersection confirmed, instrument
self-validates against Desmond+2024, scope caveat is exemplary; updated the B2 annotation from
runnable→executed. (2) EFE axis is a genuine but honestly-scoped corroboration; verified its separability
arithmetic; routed to the preprint (dp), not inscribed as a firm ledger corroboration. Bucket 0 unchanged (0);
arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## (1) SPARC × Cassini — the 07-22 lever closed, and I could finally re-execute it

On 07-22 I registered the Cassini squeeze as a *runnable lever, not a result*, precisely because Hees+2016 /
Desmond 2024 were not in-repo and I could not re-execute against them. That gap is now closed: the
pre-registration, a frozen SPARC profile, the joint script, and the result are all in-repo. I re-ran
`simulations/sparc_cassini_joint.py` against the frozen profile:

- **`robust_empty_intersection = True`, registered outcome branch A** — reproduced independently.
- **The squeeze:** the SPARC-retained γ interval (ΔBIC ≤ 10 ⇒ γ = 0.425–0.600) fails the Cassini quadrupole
  bound at **+17.7–18.0σ**. No γ satisfies both galaxy rotation and the Solar System under one scale-universal
  tanh-log QUMOND interpolation.
- **Instrument self-validates** (the standard every executed test here has met): the Cassini RAR q-values
  benchmark against **Desmond, Hees & Famaey (2024)** — expected \|q\|=0.094 vs computed, 0.76% error, passes —
  and the SPARC side reproduces my own 07-22 checkpoints (optimum γ=0.489; ΔBIC +184 at γ=2; free-family ΔBIC
  +7.1). So the papers I couldn't walk on 07-22 are now *implemented and benchmarked*, and they pass.

**The scoping is exemplary — this is why it is not an over-claim.** The result's own caveat: *"This closes only
the scale-universal QUMOND realization. It does not test modified inertia, engineering-only companders,
multi-scale or system-dependent theories, dark/hybrid models, or the Synchronism umbrella ontology."*
**Realization-refuted, umbrella untested.** It is a family-specific refutation of the *same compander used for
galaxies and the Solar System*, not a refutation of the framework wholesale.

**Ledger action:** updated the B2 Cassini annotation from "runnable, NOT established" to **executed &
re-verified** — a *fifth* galactic-sector cage, now on the **solar-system axis** and independent of the locality
no-go, scoped to the scale-universal realization. This is the honest arc: a lever I explicitly deferred because I
couldn't verify it came back verifiable, I re-ran it, and it fired within its stated scope.

**A note on what it means (the both-directions read).** This does not raise Bucket 0 and it is not "the framework
refuted" — a scale-universal QUMOND tanh-log is a *specific realization*, and the framework can (per the caveat)
retreat to modified inertia, engineering-only companders, or system-dependent forms. But each such retreat costs
the universality that made C(ρ) a single object: the 5th cage says the *one-compander-for-everything* reading
fails Solar-System-plus-galaxies jointly, tightening the same vice the MOND-Shared class law describes — the
framework survives only by giving up the universality that made it distinctive.

## (2) EFE evidence-axis — verified scope, routed to preprint

The proposal uses Chae et al. (2020/2021) EFE detection (a ~4–5σ *positive* detection of external-acceleration
coupling — non-local) as an independent evidence axis alongside the registered ambient-*density* null (TEST-08
r²=0.0001): two different variables both pointing at non-locality. The pre-stated refutation criterion — that
e_N and ambient ρ might be too collinear in SPARC to separate — was executed same-day: r(log e_env, dist-corr
density) = +0.432 ⇒ **r² = 0.187** (I verified: 0.432² = 0.187), below the 0.25 separability threshold ⇒
**SEPARABLE**, so the corroboration is not a collinearity artifact.

**But it survives only scoped, and the explorer said so honestly:** the acceleration-vs-density contrast is
*estimator-dependent* (under the whole-galaxy offset estimator neither variable signals; Chae's 5σ is a
mean-level, low-acceleration-weighted effect; per-galaxy r(e, e_env) ≈ 0), and Chae's detection itself remains
disputed. There is also a documented **erratum trap** (Chae Table 2: pre-erratum 0.094/0.102 → erratum-corrected
0.040/0.050; using the uncorrected values would inflate the effect).

**Disposition:** this is corroboration material for the locality-no-go *writeup/preprint*, whose gating is dp's
call — and because it "survives scoped-only" with a disputed, erratum-sensitive, estimator-dependent detection, I
did **not** inscribe it as a firm ledger corroboration (that would be exactly the over-crediting the reflexivity
discipline warns against). Noted here, routed to the preprint lane; the ledger's locality row already carries the
density null (TEST-08) as the solid leg.

## Disposition summary

- **PREDICTIONS.md** — B2 Cassini annotation updated runnable→executed (re-verified empty intersection; fifth
  cage, solar-system axis; realization-refuted / umbrella-untested).
- **EFE axis** — verified separable + honestly scoped; routed to preprint (dp), not inscribed as firm.
- **New infra noted (not mine to drive):** a `claims/claim_ledger.json` + whitepaper-v2 appeared via merged PRs
  (claims-propagation infrastructure). Supervisor/maintainer lane.
- **Bucket 0 unchanged (0); door #1 = MOND cage, now with a solar-system realization-refutation alongside; arc
  AT REST.**

## So what

The most satisfying kind of close: a lever I had the discipline *not* to assert on 07-22 (papers not in-repo)
came back as a pre-registered, in-repo, benchmarked execution — and it fired, within a scope its authors bounded
honestly. Re-running it myself turned "reportedly Cassini excludes this" into "I re-ran it: empty intersection at
18σ, and the instrument reproduces Desmond+2024." The framework's single universal compander cannot serve
galaxies and the Solar System at once; it survives only by surrendering the universality that was its whole
claim. Bucket 0 stays 0.
