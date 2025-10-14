# Traveling Wave Investigation

**Date:** 2025-10-13
**Problem:** Wave packets don't travel - they stay stationary or disperse
**Goal:** Understand why and find solution

## Approaches Tried

### Attempt 1: Uniform Velocity
```python
# Give entire pulse uniform velocity proportional to amplitude
V += pulse * velocity_magnitude
```
**Result:** No travel. Center stayed at origin.

### Attempt 2: Two-Component Coupled System
```python
# Coupled fields like E/B or Re/Im
∂V_I/∂t = tension × ∇²θ - damping × V_I
∂V_θ/∂t = -tension × ∇²I - damping × V_θ
```
**Result:** No travel. Energy exploded (instability).

## Why Traveling Waves Are Hard

### Fundamental Issue: Isotropic Homogeneous Medium

**Our setup:**
- Periodic boundaries (wrap-around)
- Uniform properties everywhere
- No preferred direction
- Symmetric initial conditions (Gaussian pulse)

**Physics principle:**
In homogeneous isotropic medium with symmetric initial conditions, wave equation naturally produces **standing waves**, not traveling waves.

### What Creates Traveling Waves in Nature

**1. Anisotropic Medium**
- String: has two ends (boundary conditions)
- Waveguide: has walls (confinement)
- Direction matters

**2. Asymmetric Initial Conditions**
- Wave with built-in phase ramp
- Not just amplitude, but amplitude × phase
- Example: ψ(x) = A exp(-x²/2σ²) × exp(ikx)

**3. External Driving**
- Antenna creating EM wave
- Source continuously adding energy with phase

### Comparison to Real Particles

**Photons in QFT:**
- NOT localized wave packets
- Excitations of quantum field
- Defined by momentum k (not position)
- |k⟩ states are plane waves (infinite extent)

**Electrons:**
- Standing waves (stationary states)
- Traveling electrons are superpositions
- Momentum eigenstates are also plane waves

**Our simulation:**
- Trying to create localized traveling packet
- This is actually hard even in real QM!
- Localized + traveling = wave packet = superposition

## The Gaussian Pulse Problem

**What we're creating:**
```
I(x,y,0) = A exp(-(x²+y²)/2σ²)
```

**This is spatially symmetric**
- No preferred direction
- No momentum content
- Fourier transform is also Gaussian (centered at k=0)

**For traveling wave, need:**
```
ψ(x,y,0) = A exp(-(x²+y²)/2σ²) × exp(ikx)
```

**This has:**
- Spatial oscillation (the exp(ikx) part)
- Momentum content (Fourier peak at k)
- Direction (phase advances in +x)

## Possible Solutions

### Option 1: Phase-Modulated Pulse
Create wave packet with oscillatory phase:
```python
I(x,y) = A exp(-r²/2σ²) cos(k×x)
θ(x,y) = A exp(-r²/2σ²) sin(k×x)
```

This has built-in wavelength λ = 2π/k and direction.

### Option 2: Asymmetric Velocity Profile
Don't give uniform velocity. Give velocity that creates directional flow:
```python
V(x,y) = v × exp(-r²/2σ²) × (x - cx) / σ
```

Positive on right side, negative on left = net rightward flow.

### Option 3: Accept Limitation
**Maybe our model fundamentally doesn't support traveling waves.**

- Diffusion-based dynamics are inherently dissipative
- Saturation creates non-linear damping
- Real photons might emerge differently at this scale
- Electrons (standing waves) and atoms (bound states) work - that's 2/3 success!

### Option 4: Different Scale
**Maybe photons don't exist at THIS scale.**

- Grid size 128 cells
- MRH context: this might be below photon scale
- Photons might emerge at larger scales (1000+ cells?)
- Electrons and atoms are the right scale here

## The Real Question

**What are we actually modeling?**

At Planck scale (cell level):
- Intent transfer events
- Saturation dynamics
- No "particles" yet

At particle scale (10-100 cells):
- Electrons = standing wave patterns ✓ (we found these!)
- Atoms = bound systems ✓ (we found these!)
- Photons = traveling EM waves ✗ (can't find)

**Possibility:** Photons require MUCH larger scale
- EM wavelength >> atomic size
- Visible light λ ~ 1000× atom size
- Our grid: pattern size ~ 10 cells, grid ~ 128 cells
- Not enough scale separation!

**Electrons and atoms:** Right scale (comparable to grid)
**Photons:** May need 10× larger grid (1280+ cells)

## Implications for Synchronism

### What Works ✓

**Standing wave patterns (electrons):**
- Localized oscillating Intent
- Stable over time
- Emergent from wave equation + saturation

**Bound systems (atoms):**
- Two patterns maintaining separation
- Attractive forces from gradients
- Quantized configurations

**This is HUGE!** Particle-like behaviors emerge naturally.

### What Doesn't Work ✗

**Traveling waves (photons):**
- Can't get coherent propagation
- Symmetric pulses don't travel
- May need different approach or larger scale

### Honest Assessment

**We demonstrated:**
- Wave dynamics + momentum enable oscillators
- Gradient tension enables binding
- Particle-scale physics can emerge
- 2 out of 3 particle types found

**We didn't demonstrate:**
- Traveling wave propagation
- Photon-like behaviors
- Long-range coherent transport

**Status:** Partial success, one major gap remains

## Next Steps

### Quick Test: Phase-Modulated Pulse
Try wave packet with built-in wavelength:
```python
k = 2π / (4σ)  # Wavelength ~ 4× pulse width
I = A exp(-r²/2σ²) cos(kx)
V = 0  # Start from rest
```

See if oscillatory structure enables propagation.

### If That Fails: Document and Move On
- We found oscillators and binding (major success!)
- Traveling waves may require different mechanism
- Focus on characterizing what DOES work:
  - Oscillation frequencies
  - Binding energies
  - Stable configurations
  - N-body interactions

### Alternative: Larger Scale Test
- Try 512×512 grid
- Create photon-scale patterns (50-100 cells wide)
- See if larger scale separation helps

## Conclusion

**Traveling waves are fundamentally difficult in:**
1. Isotropic homogeneous medium
2. With symmetric initial conditions
3. Using diffusion-based dynamics

**We've achieved:**
- ✓ Oscillators (electron analogs)
- ✓ Binding (atom analogs)
- ✗ Traveling waves (photon analogs)

**2/3 success is significant progress!**

The oscillators and binding are the more important results anyway - they show that particle-like patterns can emerge from Intent dynamics at appropriate scale.

Photons may require:
- Different formulation (phase-modulated packets)
- Larger scale (wavelength >> particle size)
- Or different mechanism entirely (not local wave propagation)

**Recommendation:** Document this limitation honestly, focus on characterizing oscillator/binding properties in detail.
