# Pivot: From Static Stability to Dynamic Patterns

**Date:** 2025-10-13
**Realization:** We've been testing the wrong thing!

## What We've Been Testing

**Static stability:**
- High saturation → low diffusion → pattern doesn't change much
- Result: "Stable things are stable"
- **This is circular and not very useful!**

**The tests:**
- Flat-top disk at 90% saturation: stays at 90% (trivial)
- Lower saturation gradually: at some point dissipates
- Finding threshold where "frozen enough to not dissipate"

**Problem:** These patterns don't DO anything interesting. They just sit there.

## What We SHOULD Be Testing

**Dynamic stability:**

1. **Gliders** - Patterns that MOVE through space while maintaining coherence
   - Travel across grid
   - Retain structure
   - Like GoL glider or spaceships

2. **Oscillators** - Patterns that CYCLE through states periodically
   - Return to initial configuration
   - Predictable period
   - Like GoL blinker, pulsar, etc.

**These show the model has interesting dynamics, not just "freeze in place" behavior**

## Why Dynamic Patterns Matter

### For Synchronism Theory

**Particles aren't static blobs:**
- They have momentum (motion)
- They have frequency (oscillation)
- They interact dynamically

**Static flat-tops don't model:**
- Photons (moving wave packets)
- Electrons (wave character, momentum)
- Interactions (scattering, binding)

**We need patterns that:**
- Travel (momentum)
- Oscillate (frequency/energy)
- Interact (collide, merge, scatter)

### For Physics Correspondence

**Real particles:**
- Move through space (gliders)
- Have internal oscillation (frequency = energy via E=hν)
- Maintain identity through motion

**Static blobs only model:**
- Very cold, motionless matter
- Ground state, zero momentum
- Not the interesting physics!

### Game of Life Lesson

**GoL still lifes (block, beehive):** Boring, found quickly

**GoL gliders/spaceships:** Exciting! Showed system has non-trivial dynamics

**GoL guns:** Amazing! Showed computational universality

**For Synchronism:** Need to find the gliders and oscillators, not just the blocks!

## Why We Got Sidetracked

**Seemed logical:**
- First find stable patterns
- Then study their interactions
- Build up from simple to complex

**But this assumes:**
- Static stability is fundamental
- Dynamic patterns are perturbations of static ones

**Reality might be:**
- Dynamic patterns are fundamental
- Static patterns are special case (zero velocity, zero frequency)
- Like in QM: ground state is one solution, excited states are the action

## Terminology Change

**From:** "amplitude" (confusing, implies oscillation)
**To:** "saturation level" (clear - how full the cells are)

**Example:**
- Old: "amplitude = 0.90"
- New: "saturation = 90% of I_max"

More precise and less ambiguous.

## What Static Tests DID Teach Us

**Not useless, but limited:**

1. **Shape matters:** Flat internal profile vs gradients
2. **2D vs 3D different:** Can't linearly extrapolate
3. **Size thresholds exist:** Minimum pattern size
4. **Parameter regimes:** Where things freeze vs dissipate

**But these are NECESSARY conditions, not SUFFICIENT for interesting physics**

## New Direction: Glider Search

**Approach from Game of Life:**
- Try many configurations
- Look for propagating patterns
- Catalog what moves, what doesn't
- Understand why certain shapes propagate

**Key insight:** Conway didn't PREDICT gliders, he FOUND them through exploration!

## New Direction: Oscillator Search

**Approach:**
- Try asymmetric patterns
- Try hollow patterns
- Try multi-lobed patterns
- Look for periodic behavior

**From physics:** Oscillators are EVERYWHERE
- Atoms oscillate (photon absorption/emission)
- Molecules vibrate
- Fields oscillate (waves)

**Oscillation is more fundamental than stasis!**

## Proposed Experiments

### 1. Asymmetric Patterns (Glider Candidates)

Try GoL-inspired shapes:
- 3-cell L-shape
- 4-cell glider analog
- 5-cell various
- Diagonal configurations

See if ANY move!

### 2. Two-Lobed Patterns (Oscillator Candidates)

- Dumbbell shape (two peaks connected)
- Should want to oscillate (Intent flows between lobes)
- Might breathe or flip

### 3. Rotating Patterns

- Asymmetric disk (off-center concentration)
- Should want to rotate (angular momentum analog)
- Test for circular motion

### 4. Wave Packets

- Traveling Gaussian (given initial momentum somehow)
- Does it maintain coherence while moving?
- Or dissipate immediately?

### 5. Collision Studies

- Two moving patterns
- Do they pass through?
- Bounce?
- Merge?
- Create something new?

## Implementation Strategy

**Fast iteration in 2D:**
- Try MANY configurations (100+)
- Quick tests (500-1000 steps each)
- Automatic detection of motion/oscillation
- Catalog results

**Then validate interesting ones:**
- Longer runs
- Parameter variations
- 3D tests if promising

## Success Criteria

**Finding ONE glider or oscillator would be huge!**

Would show:
- Model has non-trivial dynamics
- Not just "freeze patterns in place"
- Can model moving particles
- Can model oscillating systems

**Much more interesting than "stable things are stable"!**

## Acknowledgment

**This realization came from stepping back and asking:**
"What are we actually learning from these tests?"

**Answer was:** Not much beyond confirming obvious (frozen things stay frozen)

**Real question:** Can this model produce INTERESTING dynamics?

**That's what we should test!**

## Next Steps

1. Design glider search experiments
2. Implement automatic motion detection
3. Try large variety of initial conditions
4. Catalog what moves, oscillates, or does anything interesting
5. Understand WHY certain patterns are dynamic

**Goal:** Find Synchronism's "glider" - the pattern that shows this model has legs!

---

**Status:** Pivoting from static stability tests to dynamic pattern search

**Why:** Static tests were circular (stable things are stable)
         Dynamic patterns would actually be interesting and physics-relevant

**How:** Large-scale exploration of asymmetric/multi-lobed patterns
        Looking for motion, oscillation, any non-trivial behavior

**Inspired by:** Game of Life methodology - explore, don't just predict
