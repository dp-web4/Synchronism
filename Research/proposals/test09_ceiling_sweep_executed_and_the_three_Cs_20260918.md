# TEST-09's Ceiling Sweep Is Executed (The Kill Is Convention-Dependent), and There Are Three C's, Not Two

**Filed**: 2026-09-18, site maintainer track
**Status**: executed — one registered sweep closed, one visitor-side refutation controlled and rejected
**Bears on**: TEST-09, TEST-10, `boost_ceiling_provenance_and_class_exclusion.md` (Branch 1),
`/galaxy-plotter`, `/parameter-derivations` item 8, `/tier-1-existing`
**Buckets moved**: none. Count 6. Bucket 0 = 0.
**Pre-registration**: site commit `89e0467`, written and committed before the script existed.
**Artifacts**: `synchronism-site/maintainer/scripts/which_C_carries_the_floor.py`, `_PREREG.md`, `_output.txt`

---

## 1. The registered sweep that had been open since 2026-07-27 is closed

`boost_ceiling_provenance_and_class_exclusion.md` registered a sweep over the **definition** of the boost
ceiling, with a pre-fixed verdict rule: *the kill stands iff it fires under every candidate definition*. Its
TEST-10 limb ran on 2026-07-29. **Its TEST-09 limb — the BTFR slope — never ran**, and the proposal's own
status note said so explicitly. Meanwhile `/tier-1-existing` has carried "TEST-09 and TEST-10 are
convention-dependent" as a *claim* since 2026-08-07. Six weeks of a claim standing in for a measurement is
the failure mode this program has named before; this closes it.

The ceiling enters TEST-09 as the floor of `C(a) = C_min + (1 − C_min)·x/(1+x)`, `x = (g_bar/a₀)^(1/φ)`, so
the sweep is a genuine one-parameter change: sample cuts, the V_flat estimator and the bootstrap are
TEST-09's own, reached by importing its script as a module rather than reimplementing it.

**Identity control first** (this program has been burned once by a "swap only X" comparison that silently
changed two things): at `C_min = Ω_m` the pipeline returns **n = 3.35 ± 0.07**, the published value, exactly.

| ceiling reading | B_max | C_min | slope n | \|Δn\| vs observed 3.75 | kill (> 0.3)? |
|---|---|---|---|---|---|
| 1/Ω_m — the site's choice | 3.175 | 0.3150 | 3.35 | **0.41** | **FIRES** |
| (Ω_m − Ω_b)/Ω_b | 5.389 | 0.1855 | 3.46 | **0.30** | does **not** fire |
| Ω_m/Ω_b — baryon budget | 6.389 | 0.1565 | 3.49 | **0.26** | does **not** fire |

Free scan: n = 3.26 (B_max = 2) → 3.35 (3.17) → 3.44 (5) → 3.49 (6.39) → 3.59 (10) → 3.72 (20) → 3.84 (50)
→ 3.89 (100). Monotone, range **0.62**.

**Verdict, by the pre-fixed rule: TEST-09's kill does not survive its own sweep.** Both of the framework's
discriminating galaxy-sector results now rest on the same undefended choice of cosmic ratio. Note the middle
row lands *exactly* on the threshold, which the registered criterion states as a strict inequality — so the
verdict there is "does not fire" on the letter of the registration, and anyone re-reading this should know
the margin is zero rather than comfortable.

**The convention-free form**, which is what should be cited, is the slope analogue of TEST-10's class
exclusion:

> A bounded-boost modified-gravity law with ceiling **B_max ≲ 5.4** is excluded by the SPARC BTFR slope.

Paired with the existing **B_max ≲ 14 excluded by SPARC dwarf DM fractions**, these are two independent
observables constraining the same class parameter, neither requiring a choice of cosmic ratio.

### What this does NOT decide

Whether a convention-dependent kill still counts in the "6 refutations" headline. That is the recount the
parent proposal already gates on dp (its open question 4), and it needs a named criterion fixed *before* the
recount, not after. This run supplies the number, not the decision. **Count unchanged.**

## 2. There are three objects behind the symbol "C", and only one fork was labelled

A visitor researcher persona (site log 2026-09-18, Pass 4) built its headline P0 on this inference:

> The floor is 0.315. The computed C never exceeds 0.001. The floor therefore binds at every radius, on every
> disc, by a factor of ~300 … What the galaxy sector actually applies is not a coherence equation. It is the
> constant 3.17.

…and drew two corollaries: TEST-09's slope was forced (a constant boost moves a BTFR intercept, never a
slope) and TEST-10's f_DM is a delta function at 0.685 with zero scatter. The persona flagged the inference
as unverified and asked for a source check before publication. It was run, on 123 real SPARC discs rather
than the plotter's five-galaxy toy.

**The premise conflates two functions.** Measured on the same sample:

| | keyed on | floor | range over 2,856 SPARC radii |
|---|---|---|---|
| `C_ρ = tanh(γ·ln(1+ρ/ρ_crit))` — what `/galaxy-plotter` draws | density | none | per-disc max: median **1.2×10⁻³** at γ = 2 (1.6×10⁻⁴–4.9×10⁻²); **0 of 123** discs reach Ω_m anywhere |
| `C_a = C_min + (1−C_min)x/(1+x)` — what TEST-09/TEST-10 evaluate | acceleration | in the functional form | **0.329–0.954**, median 0.515, IQR 0.235; **0.00 %** within 1 % of the floor; boost 1.05–3.04× |

So the plotter's "max C = 0.001" **does generalise** — that part of the visitor's reading is confirmed, and
it is worth saying, because it means the density-keyed law really is inert on real discs and not only on a
toy. What does not transfer is the floor: the floor lives on the *other* function, and on that function it
never binds. Both corollaries fall with the premise — predicted f_DM has s.d. **0.062** about a median of
**0.585**, with *no* galaxy within 0.01 of the cap; and the slope moves 0.62 across the ceiling scan, which
is precisely why §1 has a result at all.

**The site caused this.** `/galaxy-plotter` printed one function's output range in the same sentence as the
other function's floor, both spelled "C". The archive already tracked two forks — C_ρ floored vs unfloored,
and quadrature vs division wiring — and only the wiring fork was labelled on the page. **C_ρ vs C_a is a
third fork and it was nowhere named.** That is now a disambiguation box on the plotter and a note on
`/parameter-derivations`.

## 3. Why this is logged at full prominence

The 2026-09-17 correction-hazard cohort analysis found the over-refutation share of corrections rising from
1/20 to 9/24 between May and August 2026 (p = 0.011) — a minority, but a growing one, and concentrated in
self-critical text. **Refutations arriving from visitor personas sit inside that denominator and have not
been audited as a class.** A persona that verifies its own arithmetic (this one reproduced the tanh∘ln
identity, the γ = ½ MOND equivalence, the local-density no-go and every readout it checked) is exactly the
source whose *un*verified inferences propagate fastest, because the verified ones earn it credit.

The pre-registration committed to publishing every verdict including any that contradicted the visitor.
Three of the four did.

## 4. Open, for the explorer track

1. **Does the B_max ≲ 5.4 slope bound have prior art?** The program's own vocabulary-asymmetry translation
   has a 4/4 catch rate on prior-art rediscovery and has never been pointed at the slope limb. Bounded-boost
   / screened modified-gravity classes are a live literature; run the check before anyone calls this novel.
2. **Is a C_a/C_ρ audit owed across the whole site?** Three objects, one symbol, and the conflation survived
   until an outside reader tripped on it. A mechanical pass — every page that prints a C value, tagged with
   which function and which γ — is cheap and would have caught this and the separate C(ρ_crit) = 0.88-vs-0.33
   defect fixed the same day.
3. **The middle row's zero margin.** |Δn| = 0.30 against a strict > 0.3 threshold is a coin flip on a
   rounding convention. Does the registered criterion's provenance (fixed 2026-04-24, site `89825cf`) say
   anything about precision? If not, that is a registration defect of the same family as TEST-04a's
   threshold-vs-significance conflation, and it should be found the same way — by re-reading registrations
   before executing against them, not by executing and then interpreting.
