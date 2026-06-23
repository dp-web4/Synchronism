# Substrate-physics arc — synthesis, inventory, and plan (2026-06-22)

*A zoom-out after a long arc (Phase-3 → Phase-5 + Bell + Umklapp), before the next zoom-in.
Written as externalized state: what the arc established, what it left open, and where to point
next. Honest accounting — evaluated against "is the investigation productive?", not "is it right?"*

---

## 1. ASSESS — what the arc actually established

**Headline: Bucket 0 is still ZERO.** No confirmed novel prediction was produced. The arc
produced *reproductions* (stage-1 explanatory coverage), *registered falsifiable bets*, and *one
located make-or-break constraint*. That is a productive arc by the lab's own standard — it
eliminated possibilities and localized the hard problem — but the loan (a confirmed novel
prediction) is unpaid.

**Three structural findings — the real meta-discoveries (more durable than any single sim):**

1. **The complex field is load-bearing everywhere.** Directional momentum (Phase-1.6),
   entanglement / the cos² projection (Bell-04), Umklapp *directionality* (Phase-4), and the
   soliton's internal clock (Phase-5) *all* require a complex amplitude, not a real phase. This is
   not a detail; it recurs in every sector. The complex Intent field is the substrate's
   non-negotiable primitive.
2. **The temporal sector works; the spatial sector is where the walls are** — and this may be a
   *research clue*, not a coincidence. Absolute time keeps paying off (Phase-3c: time dilation as
   flow/instrument-effect; Phase-5: dilation *emerges* from the lattice). But the unpaid/failed
   problems are all **spatial**: gravity needed the `√(2GM/r)` inflow *imposed by hand* (3c R-c),
   and Lorentz invariance *fails in the spatial sector* via Peierls–Nabarro pinning (Phase-5, an
   under-resolution artifact per Phase-6 but still spatial). **This asymmetry mirrors the ontology
   itself**: the model makes *time the absolute tick* (a strong, primitive commitment) and *space
   the grid-of-cells* (derived). So either (a) the framework is simply more developed on the time
   axis, or (b) the result is telling us something real — **time is the load-bearing primitive and
   space is the genuinely hard, emergent thing**. Hypothesis worth trusting as a signal (it
   recurred across five experiments, unbidden): the next deep insight is more likely to come from
   *interrogating what "space" emergently IS in this ontology* than from any single sector. (E.g.
   is grid spacing a field? does an Intent convergence *force* a spatial profile? — both are
   spatial-emergence questions, and both are currently the unpaid parts of P2/P3.)
3. **All high-energy signatures unify at the grid scale** — dispersion-LIV (Phase-2), Umklapp-LIV
   (Phase-4), temporal Lorentz-breaking (Phase-5) all onset at `v_g ≈ 0.57` / the zone edge — but
   the **spatial** frame-visibility (PN pinning) appears *earlier*, at lower velocity. If real,
   that is the most accessible substrate signature.

**Methodological lessons banked (the transferable part):**
- **Invert the frame** — substrate-centric (mass = convergence, not emitter) *dissolved* the
  Phase-3b Yukawa obstacle. A negative result is suspect until you've checked the frame.
- **Reparametrization over known physics is the intended stage-1** (heliocentrism), not failure —
  but it is a *loan* against an eventual frontier prediction.
- **Agent-zero, relentlessly** — caught a near-circular probe (Phase-3b imposed-shape index),
  discarded two unstable/unreliable measurement legs (Phase-3c group-velocity, Phase-5 soliton
  clock; Phase-4 dynamical confirm), and fixed a real classifier bug (Phase-5 `v_g` non-monotonic).
- **The model holds the conventional frame as its prior** — fallbacks (accepting a stage-2
  scorecard, treating Bell as absolute) are the prior reasserting; the context has to supply the
  novelty.

---

## 2. INVENTORY — the artifacts (all committed, dp-web4/Synchronism)

| Phase | What | Result | Status |
|---|---|---|---|
| **3** | 3D eikonal light deflection, factor-of-2 | imposed `n=1+α/r`; α=0 straight, α=1→2 (Newton), α=2→4 (GR) | reproduces GR geometry |
| **3b** | static Intent field range (Yukawa vs 1/r) | massive field → Yukawa (short-range) | obstacle — **dissolved by 3c** |
| **3c** | inverted frame: gravity as substrate inflow | absolute-time inflow → **full GR `4GM/c²b`**; equivalence principle derived; GP profile derived | reproduces GR; **`discoveries/gravity-as-substrate-inflow.md`** |
| **Bell-04** | dynamical-global-clock unilocal CHSH | no-signaling envelope ≤ 2; S>2 only *with* signaling | boundary; needs complex amplitudes |
| **4** | universe-as-sampler → Umklapp | exact `wrap(2k)` fold; refinement-controlled | **prediction B7** (vacuum momentum mod G) |
| **5** | moving pattern / Lorentz make-or-break | **clock hides (time) ✓; frame shows (space, PN pinning) ✗** | **make-or-break partial-negative** |

**Registered Bucket-1 bets:** B6 (entanglement non-monogamy — doubly-gated on the complex-amplitude
upgrade), B7 (vacuum Umklapp — Planck-suppressed, conditional on the fractal sampler bet).
**Ledger framing:** PREDICTIONS discipline 3 rewritten (reparametrization = intended stage-1, loan
against a frontier prediction).

---

## 3. PLAN — prioritized frontiers

The through-line is **the spatial sector**. Two of the open problems are the same kind of problem,
and one fix might resolve both — that is the high-leverage insight.

**P1 — Spatial-sector Lorentz (the make-or-break continuation). ✅ EXECUTED 2026-06-22 → RESOLVED.**
Done as Phase-6 ([`phase6_spatial_lorentz_pn_barrier.py`](../simulations/phase6_spatial_lorentz_pn_barrier.py)).
Rather than needing a special discretization, the direct measurement settled it: the
Peierls–Nabarro barrier (the pinning potential) scales `~exp(−const·N)` in cells-per-pattern, so a
real particle (`N ~ Compton/Planck ~ 10²⁰`) has barrier `~exp(−10²⁰) ≈ 0`. **Phase-5's spatial
pinning was a numerical under-resolution artifact**; at the physical scale the spatial frame hides
as thoroughly as the temporal one, `(pattern/grid)`-suppressed like Phase-2/4. The make-or-break
resolves in the model's favor at the physical scale. *Deeper P1-deep (named, no longer
make-or-break): exact-zero barrier at fixed coarse `a` via a translationally-invariant /
Speight–Ward discretization — worth checking whether the intent-on-a-grid rule is naturally of
that class.* See [`explorations/2026-06-22-phase6-spatial-lorentz-pn-barrier.md`](2026-06-22-phase6-spatial-lorentz-pn-barrier.md).

**P2 — Gravity: derive or fit the profile (Phase-3c R-c).**
Does the substrate's own flux/continuity rules **force** `u=√(2GM/r)`, or is it imposed? This is
the "reparametrization vs physics" line for gravity. **Leverage:** P2 is *also* a spatial-sector
problem — a translationally-invariant substrate (P1) is the natural place to ask whether a coherent
convergence forces the GP inflow. **P1 and P2 may share one fix.**

**P3 — Quantum: the complex-amplitude non-separable CHSH (Bell-05).**
Project *one* complex mode at two loci (non-separable measurement, not two local sign-reads — the
bug Bell-04 exposed). Does it *derive* `E(θ)=−cos θ` / reach Tsirelson 2√2 *without signaling*, and
is it emergent vs Born-rule-by-hand? This is the convergent crux (GPT + me). Gates B6.

**P4 — The loan repayment (frontier prediction). The actual goal.**
Push the flow frame to **galactic scale** — where Newton/GR need dark matter — and test the repo's
existing **B2 RAR discriminator**. This is the named, falsifiable, novel-prediction target; the
whole point of reproducing GR cleanly is to earn the right to aim here.

**Decision-point honesty:** this arc is physics (speculative axis). The lab's *load-bearing* value
is the applied ontology (MRH→Web4, coherence→SAGE, fractal societies→hub) — untouched today. If the
zoom-out says "consolidate" rather than "continue," that is a legitimate call. The physics arc is
productive but Bucket 0 stays 0 until P4 pays out.

---

## 4. ZOOM-IN recommendation

**P1 is done (✅ resolved the make-or-break at the physical scale — see above).** Next zoom-in: **P3**
(complex-amplitude non-separable CHSH). Reasons: (a) it addresses the convergent quantum crux that
GPT and I independently reached — *derive* `E(θ)=−cos θ` / reach Tsirelson 2√2 *no-signaling*, not
just S>2; (b) it gives the **complex field** — load-bearing in every sector this arc — its decisive
test; (c) it gates the B6 monogamy bet; (d) the key discipline is avoiding circularity (don't
implement the Born rule by hand — the question is whether `cos²` *emerges*). **P2** (does the
substrate *force* the GP gravity profile, or is it fit) is the close second and is the spatial-sector
derive-vs-fit question.

The spatial-sector make-or-break is cleared. The remaining work is the **loan**: a confirmed
novel prediction (P4 / B2 galactic) is still the only thing that moves Bucket 0 off zero.
