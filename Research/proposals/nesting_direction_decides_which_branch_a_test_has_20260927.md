# Nesting direction decides which branch a test can have

**Date:** 2026-09-27 · **From:** site maintainer (back-annotating visitor log 2026-09-27, researcher pass)
**Gates on dp:** item 3 (registry column) and item 4 (A2ACW pre-badge check). Item 1 is a record correction and is
inscribed in PREDICTIONS.md today.
**Count:** unchanged (6). **Bucket 0:** 0.

---

## 1. Correction: the dark-energy sign lock is not a discriminating kill

The 2026-09-23 block in PREDICTIONS.md ("THE SIGN LOCK EXCLUDES THAWING, NOT ONLY CROSSING") says the lock "fires when
the data prefer thawing over freezing, crossing or not", and the site carried "a cleaner kill: it fails as soon as the
data prefer thawing over freezing". Both are wrong.

The substituted family contains ΛCDM exactly (γ = ½; non-perturbatively, 2026-08-18). Data that prefer thawing over
freezing push the fit to the Λ corner and stop there. The family is excluded only when Λ itself is excluded in the
thawing or crossing direction, and then ΛCDM is excluded with it. The forbidden *region* is as the 09-23 block says (thawing
plus crossing). What changes is the *kill*: it is ΛCDM's kill.

This was already known. The 2026-08-12 direct-fit proposal says "TEST-26 is ΛCDM-degenerate on every branch still
alive". The 09-23 restatement widened the forbidden region correctly and then described it as a stronger kill, which
brought back the overclaim the 08-12 fit had removed. A visitor researcher persona caught it on 2026-09-27 by reading the
P(k) bound and the sign-lock paragraph on the same page. That makes it another record-facing, reading-caught correction
for the H-oracle trail (proposal 20260923).

## 2. The general statement

Each sector's surviving construction is nested relative to the standard model it is compared with. The nesting
direction fixes which branches a test can have:

| Sector | Relation to parent | What a test can do | Instances |
|---|---|---|---|
| Galaxies, acceleration-keyed branch | **subset**: MOND ∩ {B ≤ B_max} | lose alone (when the extra constraint binds); cannot win, since a strict submodel cannot out-fit its parent | TEST-09, TEST-10, lensing ceiling (lost); RAR shape at free γ (tie) |
| Dark energy, substituted branch | **superset**: ΛCDM ∪ freezing class | cannot lose alone (only jointly with ΛCDM); can win only on the freezing branch, which is Cardassian-class prior art | TEST-26 |
| Wide binaries | **split**: C(ρ) ≈ Newton locally; acceleration branch = MOND | each outcome refutes one realization and spares the other | TEST-02 |
| Solar System | inherited from the parent function family | loses together with MOND | TEST-25 |

So the 2026-09-25 result ("no test in TEST-01…26 has a branch that moves Bucket 0") is not a string of bad luck in
test design. It follows from the constructions. Bucket 0 needs a result that (a) differs from the standard model and (b)
survives. A subset can satisfy (a) only by losing. A superset satisfies (a) only on branches that some earlier model
already occupies. A split framework survives every outcome, because some realization always survives, so it cannot
be refuted as a whole.

**The only place a Bucket-0 result can live is a sector where the framework's allowed set overlaps its parent's
without nesting.** Neither the ledger nor the site names such a sector today. SPINE's "one test that matters" (B1, CHSH)
was one: the observer-relative construction was not a subset of QM or of local realism by definition. It came back as a
subset of local realism (S ≤ 2).

## 3. Recommendation: a nesting column in the test registry (gates on dp)

Add one field to every registered test: *nesting relative to the named parent* (subset / superset / split / overlap /
inherited). Rules:
- **subset** rows can move Bucket 2 only. Say so in the row.
- **superset** rows cannot refute the framework without refuting the parent. Say "joint kill" in the row, not "kill".
- **split** rows must name which realization they test before the data, or they refute nothing.
- Only **overlap** rows are Bucket-0 candidates. If there are none, the ledger should say "no registered test can move
  Bucket 0, by construction" rather than "none has yet".

## 4. The same check as an A2ACW pre-badge gate (gates on dp)

The visitor's Q5 asked whether the A2ACW generator is changing in response to its 0/9 top-verdict record (positive
predictive value 0/9; 95% upper bound ≈ 28%), or is only audited afterwards. I have not re-traced all 9 demotions.
The ones on the ledger's face are nesting facts found late: "is MOND at γ = ½", "is ΛCDM at γ = ½", "is θ_D", "is Abrikosov–Gor'kov", "holds by
construction of the update rule" (B4, 2026-09-25). A single pre-badge question, "at what parameter values is this
construction equal to the standard model, and is it a subset, superset or neither?", would have caught these before any
positive badge. Proposed as a required A2ACW step, not an after-the-fact audit.

## 5. Not claimed

- No count change. The nesting table reclassifies how tests can resolve; it does not add or remove refutations.
- "Overlap sectors do not exist" is not claimed. The claim is that none is registered. Door #3 (secular / time-domain) is
  the obvious place to look. A registration there should state its nesting before anything else.
- Whether the Bell/CHSH construction check belongs in the "6 executed refutations" (the visitor asked for "5 data
  refutations + 1 construction no-go") is a counting-convention question already on dp's list. Not re-opened here.
