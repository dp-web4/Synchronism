# Proposal — TEST-02's amplitude band comes from a knee nobody uses; and the A2ACW null still has no positive control

*Written 2026-09-10 by the site maintainer track, during WAKE, **before** any site edit.*
*Routed to dp. Nothing in the ledger is moved by this document.*

**Source:** visitor pass 2026-09-10 (`visitor/logs/2026-09-10.md`, Pass 4 "leading-edge researcher"),
plus my own arithmetic. Two of the visitor's three P0/P1 items turn out to be real research-level
questions rather than site-copy problems; a third is a **false premise** I want on the record so the
next reader does not re-derive it.

---

## Item 1 (new, verified here) — TEST-02's `0.05–0.4%` band requires a ρ_crit that appears nowhere else in the framework

### What the site says

`/tier-1-existing` TEST-02 states, in one block:

- mechanism: `g_eff = g_N / C(ρ)`, "low ρ → low C → larger boost";
- amplitude: "~0.05–0.4% velocity deviation (Newtonian null level)";
- verdict: "~80× below Gaia DR3 reach → practically untestable, not just difficult".

The visitor read the first two as contradictory ("these are reciprocals") and asked us to "pick a
coupling." **That diagnosis is wrong, and the real one is sharper.** The coupling is not ambiguous —
it is `g_N/C`. What is unstated is **at which ρ_crit the 0.05–0.4% was evaluated**, and the answer is:
at none of the ones the framework uses.

### The arithmetic

Script: `simulations/test02_amplitude_is_knee_conditional.py` (+ `_output.txt`; also mirrored in
the site repo at `maintainer/scripts/`). Local disc density `ρ_local = 0.09 M☉/pc³`;
`C = tanh(γ ln(ρ/ρ_crit + 1))`; acceleration excess `1/C − 1`, velocity excess `1/√C − 1`.

| ρ_crit (M☉/pc³) | provenance | γ | C(ρ_local) | velocity excess |
|---|---|---|---|---|
| 1.52×10³ | `0.029·V_flat²`, MW V_flat = 229 — the **published** calibration (`/galaxy-plotter`, `/key-claims`) | 0.489 | 2.9×10⁻⁵ | **+1.8×10⁴ %** |
| 1.52×10³ | same | 2 | 1.2×10⁻⁴ | +9.1×10³ % |
| 0.161 | measured velocity-blind knee (2026-08-27, GC fork) | 0.489 | 0.214 | +116.3 % |
| 0.161 | same | 2 | 0.710 | +18.64 % |
| 8.3×10⁻³ | Refracted Gravity E0 knee (explorer 2026-09-09) | 0.489 | 0.836 | +9.4 % |
| 8.3×10⁻³ | same | 2 | 0.9999 | +0.0051 % |
| 3.2×10⁻⁴ | bottom edge of the ρ_c grid SPARC refutes | 0.489 | 0.992 | +0.40 % |

**The knee window that actually produces the quoted band** (velocity reading):
`ρ_crit ∈ [3.8×10⁻⁵, 3.2×10⁻⁴]` at γ = 0.489, or `[1.6×10⁻², 3.0×10⁻²]` at γ = 2.
(Acceleration reading: `[1.9×10⁻⁵, 1.6×10⁻⁴]` and `[1.3×10⁻², 2.4×10⁻²]`.)

Both windows are **disjoint from every knee in use**: seven orders below the published
`0.029·V_flat²` calibration, ~3 orders below the measured 0.161, and — at γ = 2 — a factor ~2–4 above
RG's 0.0083, which yields 0.005%, ten times under the band's own floor. The one near-coincidence is
that the γ = 0.489 window's upper edge (3.2×10⁻⁴) is *exactly* the bottom edge of the ρ_c grid that
`PREDICTIONS.md` (2026-09-09) records as killed on SPARC by placement.

### Why this matters for the ledger, not just the page

The "practically untestable" verdict — and therefore the decision **not to execute** TEST-02 — is
carried entirely by the soft amplitude. At the framework's own published calibration the local
prediction is a factor ~3.5×10⁴ boost in `g`, which the Oort limit and Solar-System ephemerides exclude
by orders of magnitude. That is not a new refutation: it is the **same** statement `/key-claims`
already publishes as *"rotation curves under this reading do not fail to flatten — they blow up"* and
as the AQUAL-1984 vacuum-singularity attribution. What is new is that it applies **in the solar
neighbourhood, where ρ is well measured and not small**, so the vacuum-limit framing is not needed to
get there.

**Decision requested (gates on dp):** does TEST-02's *density branch* convert from
"unexecutable / non-discriminating" to **executed and excluded locally**, at the published
calibration? Three readings, all defensible:

- **(a) No ledger change; site copy only.** TEST-02 stays a non-discriminator; the page gains the
  arithmetic and the sentence "the quoted band assumes a knee ~7 orders below the published
  calibration." This is what I have implemented site-side, pending your ruling.
- **(b) Scope it.** TEST-02 splits into TEST-02ρ (density branch — excluded locally at the published
  knee, ~10⁴×) and TEST-02a (acceleration branch — still waiting on the Chae/Banik adjudication).
  The refutation count is unchanged because the exclusion is the already-counted vacuum-singularity
  result evaluated at a new radius.
- **(c) Count it.** A seventh executed refutation. **I recommend against (c)** — over-refutation is
  this program's live failure mode (the 09-09 retraction of the "Oort ∩ GC disjoint" no-go was
  exactly this shape), and the physical content here is not independent of the result already booked.

**My reading: (b), count unchanged.** It is the one that neither hides the number nor double-counts it.

---

## Item 2 (not new; still unrun) — the A2ACW null has no positive control, and the 1.4% yield is uninterpretable without one

`/for-researchers` citable artifact 2 reports: 3,308 sessions → ~47 internally-consistent candidates
(1.4%) → 6 externally audited → **0 survived**; Youden's J = 0 with CI [−0.46, +0.46]; Clopper–Pearson
admits a true survival rate up to ~0.39. The page already says all of this, and badges it
`audited-negative / Registered Null — Pending Cross-Vendor Control`.

What no page says is that the **negative** control (can the protocol catch known demotions? 6/6) has
never been paired with a **positive** control: *feed the protocol a verified discovery published after
the models' training cutoff, citation-stripped, and measure the demotion rate on known-good physics.*

A 100% demotion rate is as suspicious as a 100% confirmation rate. Without the positive control,
nothing in this archive distinguishes:

- **H1** — the framework genuinely produced nothing novel (the reading the site publishes), from
- **H2** — an LLM challenger rewarded for finding prior art maps *almost anything* onto a corpus,
  including real discoveries.

Under H2 the 1.4% yield is a property of the protocol and says nothing about Synchronism. **And H2
would itself be the more publishable result** — a measured prior-art-illusion rate for adversarial
LLM audit is a finding about AI-assisted research methodology, independent of whether any of the
physics here holds. Per SPINE's own "what it's already good for" framing, the applied/methodological
axis is where this program has delivered; this is the one unrun experiment on that axis.

**Decision requested:** either (i) authorise the positive control as an explorer-track run (it needs
only post-cutoff papers and the existing protocol — no instruments), or (ii) rule that it will not be
run, in which case the A2ACW null should be relabelled *permanently uninterpretable* rather than
*pending cross-vendor control*, because the cross-vendor control does not address this.

I have seeded `explorer/topics/a2acw-positive-control-post-cutoff-discoveries.md` with a concrete
protocol. It does not execute until you rule.

---

## Item 3 (a false premise, recorded so it is not re-derived) — DESI DR2 **BAO** is not DESI DR2 **full-shape**

The visitor's top P0 was *"re-execute TEST-04a on DESI DR2/DR3 — DR2 (2025-03) and DR3 (2026-07-30)
have since shipped."* **That is wrong, and the site's presentation invites the error.**

- TEST-04a's registered statistic is **DR2 full-shape fσ₈ at z ≈ 0.51**. Formal DR2 full-shape
  parameter papers are **not published** (~Spring 2027); `/for-researchers` states this correctly,
  including the April 2026 PIRSA:26040071 preliminary presentation and the open prospectivity caveat.
- **DR2 BAO** (arXiv:2503.14738) *did* publish, and this archive has already used it — the
  2026-08-12 `fit_gamma_family_to_desi_dr2.py` likelihood fit is DR2 BAO + Planck priors + Dovekie SN.
- **"DESI DR3 (2026-07-30)"** does not exist in this archive or in the literature it cites; DR3 is
  referenced throughout as ~2027–2028 and only as the *proposed* TEST-26 venue. Treat the date as a
  fetch artifact.

So the registered criterion is **not** stale and **not** overdue: its trigger has not fired. But a
careful reader got the opposite impression from the site in one pass, because the same site
prominently cites "DESI DR2 (arXiv:2503.14738)" on `/honest-assessment` while burying "DR2 full-shape
growth remains unpublished" at the end of a ~2,000-word alert block on `/tier-1-existing`. **A site
that manufactures a false impression of a stale pre-registration is damaging the one methodological
claim the program rests on.** I have fixed this site-side (BAO-vs-full-shape named at every point of
use). No ledger change; recording it here because the same conflation has now been made by an
outside-persona reader and will be made again.

---

## Item 4 (new argument for an old question) — TEST-25 is framework-specific, and the 09-09 compander number is what settles it

`/honest-assessment`'s classification table calls TEST-25 (Cassini/SPARC, +17.95σ) *inherited from
MOND*: the excluded object is the RAR-preferred interpolating-function family, and MOND uses it too.
The 2026-09-10 researcher persona argued this is wrong. **I think they are right, and the explorer's
own 09-09 result is what closes the argument** — the two were produced a day apart and nobody had put
them together.

The asymmetry is that **MOND and Synchronism are not equally free to abandon the excluded function**:

- MOND's μ is a *free function*. Excluded at Cassini, MOND selects a different μ and survives. That is
  what "inherited" is supposed to mean — a failure of a shared component the other party can shed.
- Synchronism cannot shed it. Its free-γ SPARC fit lands at **γ = 0.489**, and at γ = ½ the compander
  is *identically* μ_simple — the excluded function — for every ρ_crit. **The data drive it into the
  exclusion.**
- And there is nowhere to go. The explorer's 2026-09-09 isolation of the compander *form* at a
  non-binding ceiling gives χ²/N = **108.10** vs MOND μ's **51.45** on the same 153 discs — **2.10×**.
  So leaving γ = ½ is not an escape route; it is a 2.10× penalty. **The escape hatch that makes
  TEST-25 "inherited" for MOND is closed for this framework, quantitatively.**

A refutation your own best fit drives you into, and that you cannot leave without paying 2.10×, is
framework-specific.

**Decision requested:** reclassify TEST-25 as framework-specific, making the split
**3 mechanism roots + 1 refuted registration + 1 theorem**. Executed count stays **6**; Bucket 0 stays
**0**. Note this moves *against* the framework, and note it points the opposite way from the landing
page's other correction today (which says the "6" over-advertises breadth). Both are true: the site
over-counts *coverage* and under-attributes *this particular root*. I have published the argument on
`/honest-assessment` as an open question with the classification unchanged, pending your ruling.

---

## Summary of what gates on dp

| # | Item | Ask | My recommendation |
|---|---|---|---|
| 1 | TEST-02 amplitude is knee-conditional; published knee gives ~10⁴× locally | Ledger ruling (a)/(b)/(c) | **(b)** scope into 02ρ / 02a, **count unchanged** |
| 2 | A2ACW positive control never run | Authorise, or relabel the null permanently uninterpretable | **Authorise** — cheapest remaining novel result |
| 3 | DR2 BAO ≠ DR2 full-shape | None (framing note) | Site fixed; no ledger change |
| 4 | TEST-25 classification | Reclassify inherited → framework-specific | **Yes** — the 09-09 compander 2.10× closes MOND's escape hatch; 3 roots, count unchanged |

Count stays **6**. Bucket 0 stays **0**.
