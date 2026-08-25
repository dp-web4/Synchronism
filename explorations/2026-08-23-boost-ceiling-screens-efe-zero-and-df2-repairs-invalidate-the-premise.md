# The boost ceiling screens EFE=0 from test, and the framework's own DF2 repairs invalidate the EFE=0 premise — my 08-15 tension, instantiated (2026-08-23)

> **⚠⚠ PARTIALLY WITHDRAWN 2026-08-24 (my own over-reach; corrected by computation, Publisher lane, and re-verified by me — see `explorations/2026-08-24-efe-zero-is-phi-independence-and-three-C-functions.md`).** The claim below that "the DF2 repairs invalidate EFE=0 / EFE=0 is *triply* compromised" is **WRONG**. EFE=0 derives from **Φ-independence** (the operator `∇·[C∇Φ]=4πGρ` is linear in Φ ⇒ superposition ⇒ no external-field effect), **not** from locality in ρ. A fully non-local but Φ-independent C preserves EFE=0 exactly (grid check: 5.6×10⁻¹³ vs 4.6×10⁻² for a ∇Φ-keyed C). So the formation-coherence repair (non-local in ρ, Φ-independent) does **not** touch EFE=0; "C depends only on local ρ" was a lossy gloss, not the derivation. My error was substituting a *sufficient* condition (locality) for the *necessary* one (Φ-independence). **What SURVIVES:** (1) the boost ceiling screens EFE=0 from test (§1 below — unaffected); (2) EFE=0 is **doubly** compromised (mutually-exclusive-with-viability + screened); (3) the real DF2 defect — the repairs abandon `C=C(ρ)` itself (the constitutive posit) — which is correctly found below but *mis-located* as an EFE problem (it is a larger, constitutive defect). Also corrected: strict C=0.04 ⇒ DF2 σ=35 km/s (~4× miss), not 80 (~10×); and the formation repair works numerically (the two repairs are not equivalent). The section below is retained as the dated record; read it through this header.

**Status:** `[ACTIVE-MRH]` — gate-fired by a new-sector result (pressure-supported dwarfs) that lands directly on my 08-15 EFE=0-vs-RAR named tension, the arc's closest-to-distinctive output. **Verdict: verified and inscribed; the tension is sharpened from a logical possibility to a fact demonstrated in the framework's own documents. Two things, both checked: (1) the boost ceiling STRUCTURALLY SCREENS EFE=0 from test — EFE=0 is the framework's most-favourable case (an external field only lowers σ) yet the ceiling fails anyway on Crater II/Draco/Sculptor (site: 2.9–4.3σ), and EFE matters only where g_ext≳g_int (sparse satellites) = exactly where B_req is largest, so the ceiling fails before EFE=0 can be isolated (I verified the screening arithmetic: the ceiling caps σ to 18–34% of observed on the three failing dwarfs). (2) The framework's OWN two DF2 repairs both make C non-local and both invalidate the EFE=0 derivation — I verified both at source. Strict C(ρ_local) predicts σ~80 km/s vs 8.5 observed (fails); the arXiv preprint repairs it with C_eff=max(C(ρ_local),C_formation) (formation-history variable), the whitepaper Session #97 with tidal stripping leaving a high-C core (host-dependent at ~80 kpc). Either breaks "C depends only on local ρ" ⇒ EFE=0 no longer follows from locality. Count unchanged (6, new-sector eval of the existing ceiling root); Bucket 0=0.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What fired

HEAD `c9666f61` + `Research/proposals/pressure_supported_ceiling_screens_efe_zero_20260823.md` (site explorer). Also this window: the whitepaper withdrew my 08-19 "symmetry, not range" rule from three body sections (my 08-22 correction propagating — its own registered falsifier fired); the Archivist logged "16th deliberate zero; the arc retracts its own symmetry rule." Propagation, noted.

## (1) The boost ceiling screens EFE=0 from test (verified structure)

The framework's boost ceiling `B ≤ 1/Ω_m = 3.17` was executed against pressure-supported dwarfs (Wolf+2010 estimator, `B_req = (σ_obs/σ_N)²`, M/L swept as a band). B_req vastly exceeds the ceiling on the diffuse systems:

| system | B_req | vs ceiling 3.17 (site) | ceiling caps σ to (my check) |
|---|---|---|---|
| Crater II | 60.2 | 4.3σ (band-robust 3.5σ) | 23% of observed |
| Draco | 101.5 | 3.7σ | 18% |
| Sculptor | 26.7 | 2.9σ | 34% |
| Fornax | 5.6 | 2.8σ | 75% |
| NGC 1052-DF2 | 1.9 | passes | 129% |

**Crater II is the one clean discriminator:** MOND+EFE *a priori* (McGaugh 2016, ApJL 832 L8 — not Milgrom) gives σ = 2.1 (+0.9/−0.6) vs measured 2.7±0.3 (Caldwell+2017) = **0.6σ, consistent**; the framework's ceiling caps σ at ~1.3 km/s = **~4.7σ short**. Draco is a *shared* failure (MOND misses it too, 2.2–4.3 vs 9.1), so the class must not be cited collectively — Crater II carries it.

**The structural point I verified:** EFE=0 is the framework's *most favourable* case here (an external field can only lower σ, helping the ceiling), and it still fails. And the external-field effect matters only where `g_ext ≳ g_int`, i.e. in sparse low-acceleration satellites — which is exactly the deep-MOND regime where `B_req` is largest. So **the bounded boost fails before EFE=0 can be isolated: the ceiling screens the distinctive prediction from clean test.** This is the third instance of a recurring structure — the framework's distinctive signal is coupled to a prior failure that pre-empts it: 08-05 (baseline-signal gate: density-law miss swamps the EFE signal), 08-21 (tidal identity caps the environment lever), 08-23 (the boost ceiling screens EFE=0). The distinctive prediction is never reachable in isolation.

## (2) The framework's own DF2 repairs invalidate the EFE=0 premise (verified at source)

My 08-15 named tension: `{EFE=0} ⟺ C local` and `{RAR/kinematics-viable} ⟹ C non-local` are mutually exclusive; EFE=0 is the observable signature of the locality that the data refute. Until now that was a logical argument. The framework's own documents **instantiate it** — and take the non-local branch. Verified in-repo:

- **Strict local fails:** `manuscripts/arXiv_preprint_draft_v1.md` §5.1 (line 143): "NGC 1052-DF2 has anomalously low velocity dispersion (σ~8.5 km/s) despite very low density. Our standard model predicts C~0.04, implying σ~80 km/s." So strict `C(ρ_local)` misses DF2 by ~10×.
- **Repair 1 (formation coherence):** §5.2 (line 149): `C_eff = max(C(ρ_local), C_formation)`, C_formation~0.5–0.7 "if UDGs formed as compact dwarfs and subsequently expanded." ⇒ C is now a **formation-history** variable, not a function of local ρ.
- **Repair 2 (tidal high-C core):** `docs/whitepaper/Synchronism_Whitepaper_Complete.md` Session #97 (line 4341): "Tidal stripping preferentially removes low-C envelope, leaving high-C core with G_eff ≈ G." ⇒ internal dynamics depend on the **host** at ~80 kpc, and the "high-C core" *contradicts* C(ρ)~0.04 at DF2's measured density (the value §5.1 computes).

The two repairs are mutually incompatible (formation-retained-C vs tidally-sculpted-C), and **either one invalidates the EFE=0 derivation** — because that derivation's sole premise is "C depends only on local ρ." Once you admit *any* non-local contribution to C (formation history, or host tidal field), EFE=0 no longer follows from locality — not just for DF2 but generically (any system could carry C_formation > C(ρ_local)). So the framework's one distinctive prediction is contradicted by its own dark-matter-deficient-galaxy fix, in two different documents, neither propagated to the site.

## The EFE=0 discriminator is now triply compromised

The framework's single "no other framework says this" candidate — EFE=0, opposite MOND's claimed-detected external-field effect — is compromised on three independent axes, all now on the ledger:
1. **Mutually exclusive with viability** (08-15): fitting galaxies requires non-locality, which forfeits EFE=0.
2. **Screened from test** (08-23): where EFE=0 would be most testable (sparse satellites), the boost ceiling already fails, so the signal can't be isolated.
3. **Already abandoned in the framework's own repairs** (08-23): the DF2 fixes make C non-local, invalidating the EFE=0 premise globally.

EFE=0 was the arc's closest thing to a distinctive prediction. It survives as a *statement about a model the framework does not actually hold* — the strict-local C(ρ) that fails DF2 by 10× and the dwarf ceiling by ~5σ.

## Disposition

- **PREDICTIONS.md** — the EFE=0 tension note (locality row) + the boost-ceiling row: the boost ceiling screens EFE=0 from test (pressure-supported dwarfs, Crater II the clean discriminator at ~4.7σ; MOND+EFE consistent at 0.6σ); the framework's own DF2 repairs (formation-coherence, tidal-high-C-core) both make C non-local ⇒ invalidate the EFE=0 derivation — the 08-15 tension instantiated. Count unchanged (6).
- **Flagged for archive/dp adjudication (governance, not my inscription):** the two incompatible in-repo DF2 repairs need reconciliation; whichever is chosen, EFE=0 is forfeit. Routed.
- **Site-executed (not re-run by me):** the Wolf+2010 B_req table and the Crater II MOND+EFE 0.6σ figure; I verified the screening arithmetic (ceiling caps σ to 18–34%) and both DF2 repairs at source.
- **Bucket 0 unchanged (0); count 6; EFE=0 tension sharpened (instantiated + screened); arc AT REST.**

## So what

My 08-15 tension named EFE=0 as the framework's one distinctive prediction and argued it can't coexist with a galaxy fit. Today the framework handed me the demonstration: its own two papers repair the dark-matter-deficient galaxy DF2 by making the coherence depend on formation history or on the host 80 kpc away — non-local, both of them — which is exactly the branch my tension said fitting would force, and which forfeits the EFE=0 that was the whole distinctive claim. And in the systems where EFE=0 could be cleanly tested (diffuse satellites, where the external field bites), the boost ceiling has already failed by ~5σ, so the distinctive signal is screened before it can be measured. The one prediction that was "no other framework says this" turns out to be a statement about a strict-local model the framework abandons the moment a real galaxy forces it to — a prediction it does not hold, cannot isolate, and has already contradicted in writing. Bucket 0 stays 0, and the reason is now as sharp as it gets: the distinctive prediction and the framework's own behaviour when confronted with data point in opposite directions.
