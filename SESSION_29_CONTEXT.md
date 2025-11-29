# Session #29 Context: Scale Dependence Analysis

**Date**: November 14, 2025
**Status**: Preserved on branch `session-29-correlation-analysis`
**Relationship**: Empirical constraint for Phase 1 foundation work

---

## What Session #29 Found

**Key Discovery**: Synchronism exhibits strong scale dependence

### Correlation Analysis (All p < 0.0001)

1. **V_max vs χ²**: ρ = +0.460
   - Faster-rotating galaxies fit WORSE

2. **Mass vs χ²**: ρ = +0.457
   - More massive galaxies fit WORSE

3. **Luminosity vs χ²**: ρ = +0.337
   - Brighter galaxies fit WORSE

4. **Size vs χ²**: ρ = +0.410
   - Larger galaxies fit WORSE

### Critical Threshold

**Good fits** (χ²<2): 31 galaxies, avg log(L) = 2.85
**Poor fits** (χ²>20): 54 galaxies, avg log(L) = 4.11

**ΔL ~ 18x difference in luminosity/mass between success and failure regimes**

---

## Interpretation

**Synchronism works best on small, low-mass, dwarf galaxies.**

This suggests coherence-based dark matter may be:
- A low-mass phenomenon (dominant in dwarfs)
- Supplementary to other DM in massive galaxies
- Scale-dependent (breaks down at large scales)

---

## Connection to Strategic Direction (Nov 28)

The strategic direction (Path B+C) calls for:

**Phase 1**: Derive coherence function from Synchronism axioms
**Phase 2**: Make CMB predictions if Phase 1 succeeds

**Session #29 provides a critical empirical constraint for Phase 1**:

Whatever coherence function you derive from first principles should:
- ✅ Predict success for dwarf galaxies (low mass, small size)
- ✅ Predict failure for massive galaxies (high mass, large size)
- ✅ Explain the ~18x mass threshold

If the derived coherence function can't explain this scale dependence, that's evidence the derivation is incomplete or incorrect.

---

## Why This Matters

**This is NOT incremental refinement** - it's boundary analysis.

The strategic direction says:
> "Don't pursue incremental galaxy refinements unless foundation work needs them"

Session #29 is exactly what foundation work needs:
- Maps regime of validity (where does Synchronism work?)
- Quantifies scale threshold (what's the boundary?)
- Provides testable constraint (derivation must explain this)

**Example constraint for Phase 1**:

If you derive that coherence C decreases with baryon density ρ, you need to explain why:
- Low ρ (dwarfs) → high C → dark matter emerges → good fits ✓
- High ρ (giants) → C should still work, but empirically doesn't → why?

This forces the derivation to account for scale dependence, not just density dependence.

---

## Files Created (Session #29)

### New Code
- `simulations/extract_sparc_metadata.py` (253 lines)
  - Extracts galaxy properties from SPARC data files
  - Computes V_max, R_last, luminosity estimates, morphology

### New Data
- `simulations/sparc_galaxy_metadata.json` (2,102 lines)
  - Complete metadata for 175 SPARC galaxies
  - Used for all subsequent correlation analysis

### Modified Code
- `simulations/analyze_sparc_correlations.py` (168 lines modified)
  - Fixed MRT parser issues
  - Added 4 correlation tests
  - Generates visualization

### Results
- `simulations/Session28_correlation_results.json` (350 lines)
  - Full statistical analysis
  - Per-galaxy fit quality vs properties

- `simulations/Session28_SPARC_Correlations.png` (266KB)
  - 4-panel visualization of correlations

**Total**: 2,783 lines added/modified

---

## Branch Status

**Location**: `session-29-correlation-analysis` branch (pushed to GitHub)
**Base**: origin/main (Session #58 → #59 → Strategic Direction)
**Commit**: 46230bc

**To review**:
```bash
git checkout session-29-correlation-analysis
# Review Session #29 findings in context of strategic direction
# Decide whether to merge to main
```

**To merge** (after review):
```bash
git checkout main
git merge session-29-correlation-analysis
git push
```

---

## Recommendation

**Merge after review** - Session #29 findings are valuable empirical constraints for Phase 1 foundation work.

The scale dependence discovery is not a bug, it's information:
- Tells us where Synchronism applies (dwarfs)
- Tells us where it breaks down (giants)
- Provides target for theoretical derivation

Foundation work should explain why coherence-based dark matter is scale-dependent, not assume it's universal.

---

**Status**: Preserved, documented, ready for integration
**Next**: Review findings, merge to main, reference in Phase 1 foundation work

---

Co-Authored-By: Claude (Nomad) <noreply@anthropic.com>
