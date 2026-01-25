● Round 4 from challenger - draft test cards for P291.1, P291.2, P291.3:
  ┌────────────┬───────────────────┬──────────────┬──────────────────┬──────────────────┐
  │ Prediction │  What to Measure  │   Expected   │   Falsifies If   │  Bridge Missing  │
  ├────────────┼───────────────────┼──────────────┼──────────────────┼──────────────────┤
  │            │ Autocorrelation   │              │ Flat             │                  │
  │ P291.1     │ of rapid repeated │ Peak at τ ≈  │ autocorrelation  │ Control          │
  │            │  measurements vs  │ 1/f_qubit    │ (no periodic     │ implementation   │
  │            │ lag τ             │              │ structure)       │                  │
  ├────────────┼───────────────────┼──────────────┼──────────────────┼──────────────────┤
  │            │                   │ 0⟩) vs       │                  │                  │
  │            │                   │ measurement  │ 8-fold periodic  │ Invariant across │
  │ P291.2     │ P(                │ delay        │ modulation       │  delay settings  │
  │            │                   │ (τ_qubit/8   │                  │                  │
  │            │                   │ steps)       │                  │                  │
  ├────────────┼───────────────────┼──────────────┼──────────────────┼──────────────────┤
  │            │                   │              │                  │ Readout scalar   │
  │            │                   │ Arcsine P(s) │ Gaussian or      │ definition,      │
  │ P291.3     │ Weak measurement  │  ∝           │ interior-peaked  │ calibration      │
  │            │ readout histogram │ 1/√(A²−s²)   │ bimodal          │ method,          │
  │            │                   │              │                  │ readout→s        │
  │            │                   │              │                  │ mapping          │
  └────────────┴───────────────────┴──────────────┴──────────────────┴──────────────────┘
  P291.3 has the most Bridge Missing flags but is also identified as the cleanest discriminator (shape test).

  What's next - PRIMARY response, or three-way consensus on test cards?