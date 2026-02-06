# Synchronism: Supplementary Code Material

This directory contains the code supplementary material for:

**"Synchronism: Dark Matter Phenomenology from Quantum Coherence in Galactic Systems"**

## Contents

### Main Code Module

- `synchronism_validation_code.py` - Core validation routines

### Key Functions

| Function | Description |
|----------|-------------|
| `SynchronismModel.predict_dm_fraction()` | Core DM fraction prediction |
| `SynchronismModel.predict_from_galaxy_params()` | Galaxy-level prediction |
| `validate_star_cluster()` | Star cluster validation |
| `validate_galaxy_cluster()` | Galaxy cluster with ICM correction |
| `icm_coherence()` | ICM coherence calculation |
| `parameter_sensitivity()` | Parameter sensitivity analysis |

## Model Parameters

| Parameter | Value | Status | Source |
|-----------|-------|--------|--------|
| A | 0.028 M☉/pc³ | SEMI-DERIVED | Jeans criterion |
| B | 0.5 | SEMI-DERIVED | R ∝ V^0.75 scaling |
| γ | 2.0 | DERIVED | Decoherence theory |

## Installation

```bash
# No dependencies beyond numpy
pip install numpy

# Run demonstration
python synchronism_validation_code.py
```

## Usage Examples

### Basic Prediction

```python
from synchronism_validation_code import SynchronismModel

model = SynchronismModel()
f_DM, C, rho_crit = model.predict_dm_fraction(rho=0.01, V=200)
print(f"f_DM = {f_DM:.4f}")  # ~0.95 for spiral galaxy
```

### Galaxy Validation

```python
f_DM, C, info = model.predict_from_galaxy_params(
    M_star=1e10,  # M_sun
    R_eff=3.0,    # kpc
    V=220         # km/s
)
print(f"Regime: {info['regime']}")  # low_density → DM-dominated
```

### Star Cluster

```python
from synchronism_validation_code import validate_star_cluster

result = validate_star_cluster(M=1e6, R_pc=5, sigma_v=10)
print(f"f_DM = {result['f_DM_pred']:.4f}")  # ~0.0 (no dark matter)
```

### Galaxy Cluster with ICM

```python
from synchronism_validation_code import validate_galaxy_cluster

result = validate_galaxy_cluster(
    M=1e15, R_200=2000, sigma_v=1000,
    f_icm=0.15, T_icm=8.5e7, n_e=3e-3, B_uG=5,
    f_DM_obs=0.87
)
print(f"Error reduction: {(1 - result['error_corrected']/result['error_orig'])*100:.1f}%")
```

## Validation Summary

| System Type | N tested | Success Rate | Mean Error |
|-------------|----------|--------------|------------|
| Rotation curve galaxies | 160 | 99.4% | 3.2% |
| Early-type galaxies | 10 | 70% | 14.1% |
| Star clusters | 19 | 100% | 0% |
| Galaxy clusters (w/ ICM) | 6 | 100% | 3.2% |

## Cross-Scale Validation

The model is validated across 13 orders of magnitude:
- 10² M☉ (open clusters) → 10¹⁵ M☉ (galaxy clusters)

## References

- Session #52: Parameter recalibration
- Session #53: Theoretical derivation of A, B
- Session #54: Star cluster validation
- Session #55: Galaxy cluster testing
- Session #56: ICM coherence, sensitivity analysis

## License

MIT License

## Contact

[TBD for publication]
