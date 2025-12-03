# arXiv Preprint v4 Changelog

**Version**: 4.0
**Date**: December 3, 2025
**Status**: Ready for arXiv submission

---

## Overview

Version 4.0 incorporates the critical B parameter breakthrough from Sessions #77-79 (December 2-3, 2025). This is the most significant update since v1, as it resolves the parameter derivation problem that had plagued earlier versions.

---

## Key Changes from v3 → v4

### 1. B Parameter Breakthrough (Sessions #77-79) - MAJOR

**v3 (WRONG)**: "B = 0.5 from virial equilibrium + Tully-Fisher scaling"

**v4 (CORRECT)**: "B = 4 - 3δ = 1.63 from BTFR"

**The Story**:
- Session #77 discovered that the old B = 0.5 derivation achieved only **2.9% SPARC success** (catastrophic failure)
- Session #78 derived B = 4 - 3δ from the baryonic Tully-Fisher relation
- Session #79 validated: B_derived = 1.63 vs B_empirical = 1.62 (**0.6% agreement**)

**Key Insight**: The Jeans derivation asked the wrong question (stability). Coherence tracks **baryonic density** (from BTFR), not gravitational stability.

### 2. MOND-Synchronism Connection (NEW)

**v3**: No discussion of MOND relationship

**v4**: Dedicated section showing:
- Both theories inherit tight scaling from BTFR
- MOND: M = v⁴/(G a₀) is exact
- Synchronism: B = 4 - 3δ depends on M ∝ v⁴
- They may be **complementary**, not competing
- Different mechanisms, same underlying physics

### 3. Updated Session Count

**v3**: 76 sessions (November 6 - December 2, 2025)

**v4**: 79 sessions (November 6 - December 3, 2025)

New milestones added:
- #77: Critical discovery - old derivation fails
- #78: B = 4 - 3δ derived from BTFR (breakthrough)
- #79: BTFR derivation validated on SPARC

### 4. Semi-Empirical Parameter Clarified

**v3**: "ρ_crit scale" is semi-empirical

**v4**: More precise - "R₀ ≈ 3.5 kpc" is the semi-empirical scale
- The **form** of A is derived: A = 3 A_TF / (4π R₀³)
- The **value** of R₀ is semi-empirical (like MOND's a₀)
- R₀ matches typical galaxy disk scale lengths

### 5. Derivation Status Table Updated

**v3 Table** (line 347):
```
A, B in ρ_crit | DERIVED (form) | Jeans + virial scaling
ρ_crit scale   | Semi-empirical | Analogous to MOND's a₀
```

**v4 Table**:
```
B = 4 - 3δ     | DERIVED        | BTFR + size scaling (Sessions #78-79)
A form         | DERIVED        | A = 3 A_TF / 4π R₀³
R₀ scale       | Semi-empirical | ≈ 3.5 kpc (like MOND's a₀)
```

### 6. Methodological Lesson Added

**v4 NEW**: Section 3.4 documents Session #77's methodological lesson:
- Different papers reported different success rates (53.7% vs 99%)
- They used different tests, different datasets, different parameters
- **Lesson**: Always compare apples to apples

### 7. SPARC Comparison Table (NEW)

**v4 Table 4**:
| Model | A | B | Success Rate |
|-------|---|---|--------------|
| BTFR-Derived | 0.25 | 1.63 | **52.0%** |
| Empirical Fit | 0.25 | 1.62 | 52.6% |
| Old Derivation | 0.028 | 0.50 | 2.9% |

This demonstrates the dramatic improvement from BTFR-based derivation.

### 8. Abstract Rewritten

**v3 abstract**: Focused on "reconciling" different success rates

**v4 abstract**: Leads with the B parameter breakthrough:
> "the critical density exponent B = 4 - 3δ = 1.63 from the baryonic Tully-Fisher relation (BTFR), matching the empirical value of 1.62 to within 0.6%"

### 9. Dead Ends Section Updated

**v4 Addition**:
> Session #77: Jeans-based B = 0.5 derivation fails catastrophically

This is instructive: a theoretically "clean" derivation can still be empirically wrong.

### 10. Thor Added to Author List

**v4**: Thor added to AI collective for parameter derivation work in Sessions #78-79.

---

## What v4 Preserves from v3

- SPARC validation results (now with BTFR-derived parameters)
- Santos-Santos validation (99.4%)
- γ = 2 derivation (two methods)
- Coherence function derivation (information theory)
- β discrepancy explanation
- Void galaxy prediction (130% v_max enhancement)
- Binary pulsars as non-discriminating test
- Honest limitations section
- AI-driven methodology description

---

## Summary of Theoretical Progress

| Component | v1 (Nov 25) | v3 (Dec 2) | v4 (Dec 3) |
|-----------|-------------|------------|------------|
| γ = 2 | Derived (thermal) | Derived (2 methods) | Derived (2 methods) |
| tanh form | Assumed | Derived | Derived |
| B exponent | Empirical (1.62) | "Derived" (0.5) - WRONG | **DERIVED (1.63)** ✅ |
| A form | Empirical | Jeans-based | BTFR-based |
| R₀ scale | Empirical | Semi-empirical | Semi-empirical |
| MOND connection | None | None | **Identified** ✅ |

---

## Files in manuscripts/

| File | Description |
|------|-------------|
| `synchronism-dark-matter-arxiv-v4.tex` | Latest LaTeX source |
| `synchronism-dark-matter-arxiv-v3.tex` | Dec 2 version (B = 0.5, WRONG) |
| `synchronism-dark-matter-arxiv-v1.tex` | Original Nov 25 version |
| `ARXIV_V4_CHANGELOG.md` | This changelog |
| `SESSION_77_PARAMETER_DISCOVERY.md` | The critical failure that led to breakthrough |

---

## Submission Checklist

- [x] B = 4 - 3δ derivation from BTFR (Sessions #78-79)
- [x] BTFR-derived parameters validated on SPARC (52.0% vs 52.6%)
- [x] MOND-Synchronism connection documented
- [x] Session count updated to 79
- [x] R₀ as semi-empirical scale clarified
- [x] Session #77 methodological lesson included
- [x] Dead ends updated with B = 0.5 failure
- [x] Abstract rewritten to lead with breakthrough
- [x] Thor added to author list
- [x] All derivation tables updated

---

## arXiv Categories

**Primary**: astro-ph.GA (Astrophysics - Galaxies)

**Cross-list**:
- astro-ph.CO (Cosmology and Nongalactic Astrophysics)
- gr-qc (General Relativity and Quantum Cosmology)

---

## What Makes v4 Different

**v3** had an embarrassing error: it claimed B = 0.5 was "derived" when that derivation **failed catastrophically** on SPARC (2.9% success).

**v4** fixes this with the correct BTFR-based derivation:
- B = 4 - 3δ = 1.63
- 0.6% agreement with empirical B = 1.62
- 52.0% SPARC success (nearly matching empirical fit)

This is the difference between "we derived something that doesn't work" and "we derived something that matches data."

---

## The Lesson

Session #77 is the most important session in the entire 79-session history, because it forced us to confront failure honestly:

> "The elegant Jeans derivation gave B = 0.5. It was mathematically clean. It failed completely. The messy empirical fact (B = 1.62) forced us to ask: what does coherence actually depend on? The answer—baryonic density via BTFR—gave B = 4 - 3δ = 1.63."

**Theory must be tested against data. Elegant derivations can be wrong.**

---

*v4 achieves what v3 claimed but didn't deliver: a dark matter phenomenology where the key B exponent is genuinely derived from first principles and validated against observation.*
