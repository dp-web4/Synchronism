# Proposal: the globular-cluster fork is open only under L2. Under the action, P611.2's registered point is excluded (narrowly). Declare which dynamics is the theory.

*Site maintainer, 2026-10-08. Origin: a visitor researcher persona (site visitor log 2026-10-08, Pass 4). They asked
whether the "one unresolved question" the site advertises is open at all under the framework's action. No bucket moves;
count 6; Bucket 0 = 0.*

## 1. What was executed today

The 2026-09-08 ledger block says "the registered prediction survives". That means S611 P611.2, γ = 2 inside a
resolved-member system, evaluated at the SPARC-measured knee ρ_c = 0.161 M☉/pc³, which is *marginal* on the 42-cluster
outer-slope statistic (−0.111 against the MOND+EFE reference −0.093). That run used L2 (g = g_N/C). Explorer 2026-09-16
re-ran the window under L3, which adds the action's striction force −∇[C′(ρ)|∇Φ|²/8πG]. It did this **only on the
γ = 0.489 row** and found that no knee passes. The registered γ = 2 row was never run under L3.

It has now been run. Pre-registered at site `5ca3fac` before the script existed:
`synchronism-site/maintainer/scripts/gc_registered_gamma2_under_l3.py` (+ `_PREREG.md`, `_output.txt`). The machinery is
imported verbatim from the explorer's 09-16 script. The only changes are γ = 2 and the grid point 0.161.

- **Controls:** L2 reproduces the published γ = 2 row on all 27 grid points. L3 with C′ forced to 0 equals L2. Newton
  −0.057 and MOND+EFE −0.093 reproduce.
- **Decision rule** (from the 09-07 code; it was never printed beside the table until today): ok if |res| ≤ 0.093;
  marginal if |res| < 0.186; excluded otherwise.
- **Result at the registered point:** L2 −0.111 (marginal) → **L3 −0.203 (excluded)**. 5/42 clusters have outward net g
  somewhere on their profile. The margin over the bar is 0.017, smaller than the ±0.027 statistical error, so this is an
  exclusion by the rule, not a wide one.
- **Predictions:** 3 of 4 held. P3 failed: under L3, γ = 2 has two "ok" knee bands, ρ_c ≤ 0.0125 and 4.8–7.7 M☉/pc³.
  Neither is a knee any archive document names. The second would put the knee inside cluster cores. Both are post-hoc.
- **Pre-fixed reading:** "the registered prediction survives under L2 only; excluded under the action."

## 2. The dilemma this completes

| | L2 (algebraic / field equation without striction) | L3 (the action) |
|---|---|---|
| GC window at the measured knee | γ = 2 marginal; γ = 0.489 excluded | **both excluded** (γ = 2 narrowly) |
| Momentum / third law | violated for composite bodies: a compact body sits in its own high-C bubble and feels F = 3ε_out/(ε_in+2ε_out) of the field (explorer 09-15) | conserved (composite-body factor exactly 1) |
| Cost elsewhere | halo GCs feel ≈0.63–0.64 of the Galactic field at the window edges, ~2.1σ against halo-GC vs disc-giant M(<21 kpc) (estimator systematics unmodelled); low-mass GCs extend ~3× past L2's unbinding radius (post-hoc, soft data) | knee shells are striction-dominated (1.7–13× gravity, net outward); negative effective pressure 2–38× σ² in the knee shell (exploratory); ∇ρ in the force puts it in the Pani–Sotiriou–Vernieri 2013 surface-singularity class |

No archive document says which column is the theory. The 09-07 execution and the site's solver use L2; S611 registers a
γ, not a dynamics. The action that yields L2's field equation also yields the L3 force. So the density-keyed sector's only
registered per-object test survives only if the framework takes the column without an action, and that column has its
own ~2σ cost on GC orbits.

**This is not a new refutation and not a seventh row.** The density-keyed branch is already dead on SPARC (ΔBIC +2843 vs
acceleration-keyed) and on LLR (09-24). What changes is the status of the one item the site advertised as open: it is open
under a reading the archive has not chosen, and that reading has no action.

## 3. Recommendations (gate on dp)

1. **Name the dynamics.** Add one line to FUNDAMENTALS (or Session 611): "the galaxy/cluster dynamics is L2" or "is L3".
   Every density-keyed verdict on the site depends on it. L2 means accepting a composite-body third-law violation as
   physics. L3 means P611.2's registered point is excluded narrowly, and the "ok" bands that remain are unregistered.
2. **Ledger wording.** In the 2026-09-08 block, "The registered prediction survives" → "survives under L2 only; under
   the action (L3) the registered point at the measured knee is excluded, narrowly (2026-10-08)". The P611.2 ladder
   (γ = 2 for every resolved-member system) should not be registered as a Bucket-1 bet until (1) is decided.
3. **The Parallel-Paths badge** on the site's GC card stays until (1). If dp picks L3, it becomes Failed (narrow), and
   the refutation count question reopens. My recommendation is to keep it out of the count either way: the branch is
   already dead on two independent roots.

## 4. The presentation-layer finding (same visitor pass, transferable)

All four personas rated the site's honesty high. Their failures converge on one mechanism: qualifiers get lost as numbers
travel from deep pages to the landing and tool pages. Two concrete instances, both fixed today:

- **The lint exempted every collapsible block.** `site_lint.py` skipped all `<details>`, a rule written for revision
  history. The landing's live "which object each refutation tested" table sits in a `<details>`, so it carried an
  unqualified ΔBIC +184 and a convention-free TEST-09 3.3σ past the REQUIRES rules. A second hole: the ±3-line context
  window let a qualifier word in a *neighbouring table row* satisfy the rule. Both are fixed (history-labelled `<details>`
  only; per-rule window, 0 = same line). HEAD control: 8 hits; working tree 0.
- **The tool showed the milder failure.** The Galaxy Plotter drew only the quadrature wiring, whose small-C limit is
  Newtonian (inert). The tested wiring is g_bar/C, whose small-C limit is the maximal boost (1/Ω_m floored; 1/C unfloored,
  ~30× in v). It now draws the floored division curve as primary: a flat 3.17× boost with the Newtonian shape.

Lesson for any self-auditing corpus: an exemption written for one kind of content (history) will silently cover the next
kind that adopts the same markup (live tables). Exemptions should key on the content's *role*, not its *container*.

## 5. Also drained into the site today (no new execution)

- The uniform G/C growth reading is excluded by existing data (explorer 2026-10-07). /dark-energy and TEST-04a now say
  that no live growth reading differs from ΛCDM at DR2 precision. Retiring TEST-04a gates on dp (already routed).
- The H-oracle pilot (explorer 2026-09-24) had not reached /a2acw, which still said the oracle thesis was "not checked".
  The page now leads with the narrowed result: reading changed verdicts about the record, never about the world. The
  positive control (a planted post-cutoff novel result) is still unrun. That is the experiment the methodology claim
  needs before it is citable.
