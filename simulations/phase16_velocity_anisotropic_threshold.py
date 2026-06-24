#!/usr/bin/env python3
"""
Phase-16 (velocity-anisotropic phase-transition threshold vs isotropy bounds).

dp's phase-transition framing (Phase-15 unifying frame): a pattern holds only while its
internal frequencies ("temperature") sit within a quantized band whose thresholds are
SUBSTRATE properties (global-/CMB-frame-fixed). Translation at velocity v relative to the
substrate dilates ALL internal frequencies by S = 1/gamma = sqrt(1 - v^2/c^2) (Phase-5,
Part-A shared-S result). Measured against FIXED thresholds, the drift is:

    delta f / f = 1 - S = 1 - sqrt(1 - beta^2) ~= (1/2) beta^2     (beta = v/c)

Because v is a DIRECTED velocity relative to the substrate, the *full* directional structure
of the observed internal frequency (what a lab co-moving with the matter would see comparing
to substrate-fixed thresholds, if we resolve direction) is governed by the standard
relativistic Doppler / aberration of the boost. The leading isotropic piece is the
(1/2)beta^2 dilation (a MONOPOLE in 1-cos, i.e. velocity-magnitude only); a genuine
*anisotropic* (direction-dependent) modulation appears at the SAME order ~beta^2 once you
keep the angular dependence of the threshold-crossing condition for a moving pattern.

This script:
  (1) computes the magnitude (1/2)beta^2 for v = 370 km/s (solar system vs CMB);
  (2) shows the dipole (beta, ~1e-3) and quadrupole (beta^2, ~1e-6) angular structure;
  (3) places it against the existing isotropy / LIV bounds (clock-comparison, Hughes-Drever,
      optical-cavity Michelson-Morley) at ~1e-18 .. 1e-22 on fractional frequency anisotropy.

No claim of novelty here -- this is the magnitude/structure calc behind the adjudication doc.
"""
import numpy as np
import json
import os

c = 299792.458  # km/s
v_cmb = 370.0    # km/s, solar system barycenter vs CMB rest frame (Planck dipole)

beta = v_cmb / c

# ---- (1) magnitude of the isotropic (velocity-magnitude) threshold dilation ----
S = np.sqrt(1.0 - beta**2)          # = 1/gamma, the shared-S dilation factor
gamma = 1.0 / S
dilation_exact = 1.0 - S            # 1 - sqrt(1-beta^2)
dilation_approx = 0.5 * beta**2     # leading term

# ---- (2) angular structure of the boost-induced frequency modulation ----
# For a pattern's internal oscillation viewed across the boost, the relativistic Doppler
# factor as a function of angle theta (between v and line of sight) is
#     f_obs/f_0 = 1 / [ gamma (1 - beta cos theta) ]   (relativistic Doppler)
# Expand to O(beta^2):
#     ~ 1 + beta cos theta + beta^2 cos^2 theta - (1/2) beta^2 + ...
# -> DIPOLE term  : beta cos theta            ~ beta      ~ 1.2e-3
# -> QUADRUPOLE   : beta^2 (cos^2 theta - 1/3) (traceless) ~ beta^2 ~ 1.5e-6
# -> MONOPOLE     : (1/2) beta^2 = 1 - S       (the pure time dilation, isotropic) ~ 7.6e-7
theta = np.linspace(0, np.pi, 361)
doppler = 1.0 / (gamma * (1.0 - beta * np.cos(theta)))
mod = doppler - 1.0  # fractional deviation from rest internal frequency

# decompose amplitudes
dipole_amp = beta                       # leading direction-dependent (odd) amplitude
quadrupole_amp = beta**2                # leading even anisotropic amplitude
monopole_amp = 0.5 * beta**2            # pure isotropic dilation (== dilation_approx)

# peak-to-peak anisotropy across the sky at O(beta^2) even part:
even_part = 0.5 * (mod + (1.0/(gamma*(1.0+beta*np.cos(theta))) - 1.0))  # average fwd/back kills dipole
anisotropy_quad_p2p = even_part.max() - even_part.min()

# ---- (3) existing isotropy / LIV bounds on fractional frequency anisotropy ----
bounds = {
    "optical_cavity_MM_modern (Nagel 2015, dim-4 SME)": 1e-18,
    "Hughes-Drever / clock-comparison (nuclear/atomic)": 1e-22,
    "He-Ne / rotating cavity (Herrmann 2009)": 1e-17,
}
tightest_bound = min(bounds.values())
loosest_bound = max(bounds.values())

# orders of magnitude the predicted anisotropy sits ABOVE the bounds
oom_above_tightest = np.log10(quadrupole_amp / tightest_bound)
oom_above_loosest = np.log10(quadrupole_amp / loosest_bound)
oom_dipole_above_tightest = np.log10(dipole_amp / tightest_bound)

out = {
    "v_cmb_km_s": v_cmb,
    "beta": beta,
    "gamma_minus_1": gamma - 1.0,
    "S_dilation_factor": S,
    "isotropic_dilation_exact_(1-S)": dilation_exact,
    "isotropic_dilation_approx_(half_beta2)": dilation_approx,
    "dipole_amplitude_~beta": dipole_amp,
    "quadrupole_amplitude_~beta2": quadrupole_amp,
    "quadrupole_anisotropy_peak_to_peak": anisotropy_quad_p2p,
    "existing_bounds_fractional_freq_anisotropy": bounds,
    "tightest_bound": tightest_bound,
    "loosest_bound": loosest_bound,
    "OOM_quadrupole_above_tightest_bound": oom_above_tightest,
    "OOM_quadrupole_above_loosest_bound": oom_above_loosest,
    "OOM_dipole_above_tightest_bound": oom_dipole_above_tightest,
}

print("=" * 72)
print("Phase-16: velocity-anisotropic phase-transition threshold shift")
print("=" * 72)
print(f"v (solar system vs CMB)      = {v_cmb} km/s")
print(f"beta = v/c                   = {beta:.6e}")
print(f"gamma - 1                    = {gamma-1.0:.6e}")
print(f"isotropic dilation 1 - S     = {dilation_exact:.6e}  (== (1/2)beta^2 = {dilation_approx:.6e})")
print("-" * 72)
print("Angular structure of the boost-induced internal-frequency modulation:")
print(f"  MONOPOLE  (1/2)beta^2      = {monopole_amp:.3e}   <- isotropic time dilation")
print(f"  DIPOLE    ~beta            = {dipole_amp:.3e}   <- odd, direction of v")
print(f"  QUADRUPOLE ~beta^2         = {quadrupole_amp:.3e}   <- leading ANISOTROPY")
print(f"  quadrupole p2p across sky  = {anisotropy_quad_p2p:.3e}")
print("-" * 72)
print("Existing isotropy / LIV bounds on fractional frequency anisotropy:")
for k, vv in bounds.items():
    print(f"  {k:50s} ~ {vv:.0e}")
print("-" * 72)
print(f"Quadrupole anisotropy (~{quadrupole_amp:.1e}) sits")
print(f"   {oom_above_loosest:.1f} orders above the loosest  bound (~{loosest_bound:.0e})")
print(f"   {oom_above_tightest:.1f} orders above the tightest bound (~{tightest_bound:.0e})")
print(f"Even the DIPOLE (~{dipole_amp:.1e}) sits {oom_dipole_above_tightest:.1f} orders above the tightest bound.")
print("=" * 72)

resdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(resdir, exist_ok=True)
with open(os.path.join(resdir, "phase16_velocity_anisotropic_threshold_result.json"), "w") as f:
    json.dump(out, f, indent=2)
print(f"wrote results/phase16_velocity_anisotropic_threshold_result.json")
