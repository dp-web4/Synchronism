#!/usr/bin/env python3
"""
Session 669: The chemistry "r=0.98 validations" are the Debye model relabeled.
Runs the null comparison S651 recommended but never executed -- and corrects
S651's premise (the right null is the Debye model, not a polynomial in Z).

The framework (Research/Chemistry/Framework_Summary.md):
  - DEFINES gamma_phonon = 2*T/theta_D            (lines 179,184,428)
  - reports "Sound Velocity: v vs theta_D (r=0.984), v ~ 2/gamma"  (line 184)
  - reports BOTH "v_D vs theta_D: r=0.982" AND "v_D vs 1/gamma_phonon: r=0.982"
    as separate "EXCELLENT validations" (line 519)
  - "Heat Capacity: C_p/C_classical = gamma/2 (r=-0.988 Debye)"  (line 179)

Claim under test (S651): r=0.98 is compared against an implicit r=0 null; the
relevant null is a smooth monotonic baseline. S651 guessed "polynomial in Z" and
predicted a tie. This script shows the right null is the DEBYE MODEL (1912), and
the framework ties it not approximately but by an EXACT DEFINITIONAL IDENTITY.

Part 1 (exact, needs only theta_D): since 1/gamma_phonon = theta_D/(2T) is a
positive-linear function of theta_D at fixed T, Pearson r(X, 1/gamma_phonon) =
r(X, theta_D) for ANY property X, to machine precision. So "v vs 1/gamma" and
"v vs theta_D" are the SAME correlation; gamma_phonon adds zero information.

Part 2 (illustration): the Debye relation theta_D = (hbar/k_B)(6 pi^2 n)^(1/3) v_D
holds on real elements -> "v vs theta_D r=0.98" is the Debye model, not a
coherence discovery. So Delta_r(Synchronism - Debye) = 0, by construction.
"""
import numpy as np

hbar = 1.054571817e-34
kB   = 1.380649e-23
T    = 300.0

# Debye temperatures (K), Debye/longitudinal sound velocity (m/s, textbook),
# mass density (kg/m^3), molar mass (g/mol). Common elements; values are standard
# textbook figures (accuracy is not load-bearing -- Part 1 is exact algebra).
# columns: theta_D, v_sound, rho, M
elements = {
    "C(diam)": (2230, 13300, 3515, 12.011),
    "Be":      (1440,  8970, 1850,  9.012),
    "Si":      ( 645,  6300, 2330, 28.085),
    "Al":      ( 428,  5100, 2700, 26.982),
    "Fe":      ( 470,  4910, 7874, 55.845),
    "Ni":      ( 450,  4970, 8908, 58.693),
    "Cu":      ( 343,  3810, 8960, 63.546),
    "W":       ( 400,  4620,19250,183.84),
    "Ag":      ( 225,  2680,10490,107.868),
    "Au":      ( 165,  2030,19300,196.967),
    "Pb":      ( 105,  1260,11340,207.2),
    "Na":      ( 158,  2000,  971, 22.990),
    "K":       (  91,  1480,  862, 39.098),
}
names = list(elements)
theta = np.array([elements[e][0] for e in names], float)
vss   = np.array([elements[e][1] for e in names], float)
rho   = np.array([elements[e][2] for e in names], float)
M     = np.array([elements[e][3] for e in names], float)

# number density n (atoms / m^3)
NA = 6.02214076e23
n  = rho / (M * 1e-3) * NA

def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return np.corrcoef(a, b)[0, 1]

# Framework's coherence variable, by its own definition
gamma_phonon = 2 * T / theta
inv_gamma = 1.0 / gamma_phonon          # = theta_D / (2T)

print("=" * 72)
print("PART 1 (exact): gamma_phonon = 2T/theta_D  =>  r(X,1/gamma)=r(X,theta_D)")
print("=" * 72)
r_v_theta = pearson(vss, theta)
r_v_invgamma = pearson(vss, inv_gamma)
print(f"   r(sound velocity, theta_D)        = {r_v_theta:.6f}")
print(f"   r(sound velocity, 1/gamma_phonon) = {r_v_invgamma:.6f}")
print(f"   difference                        = {abs(r_v_theta - r_v_invgamma):.2e}")
print("   => identical to machine precision. The framework's line 519 reports both")
print("      as 0.982 because 1/gamma_phonon IS theta_D/(2T). gamma_phonon carries")
print("      ZERO information beyond theta_D. 'v vs 1/gamma' is 'v vs theta_D'.")
# verify the identity for an arbitrary unrelated property too
rng = np.random.default_rng(0)
X = rng.random(len(theta))
print(f"   sanity: random X -> r(X,theta_D)={pearson(X,theta):.6f}, "
      f"r(X,1/gamma)={pearson(X,inv_gamma):.6f} (equal for any X)")

print()
print("=" * 72)
print("PART 2 (illustration): 'v vs theta_D' IS the Debye model (1912)")
print("=" * 72)
# Debye relation: theta_D = (hbar/kB) (6 pi^2 n)^(1/3) v_D
debye_const = (hbar / kB) * (6 * np.pi**2 * n)**(1.0/3.0)   # multiplies v to give theta
theta_pred = debye_const * vss
print(f"   {'elem':>8} {'theta_D':>8} {'theta_D(Debye pred)':>20} {'ratio':>8}")
for i, e in enumerate(names):
    print(f"   {e:>8} {theta[i]:>8.0f} {theta_pred[i]:>20.0f} {theta[i]/theta_pred[i]:>8.2f}")
r_debye = pearson(theta, theta_pred)
print(f"   r(theta_D measured, theta_D from Debye formula v,n) = {r_debye:.4f}")
print("   => theta_D ~ v * n^(1/3): correlating sound velocity with theta_D recovers")
print("      Debye's 1912 relation. The 'coherence' result is the Debye model.")

print()
print("=" * 72)
print("PART 3: the corrected null and Delta_r")
print("=" * 72)
print("   S651 proposed the null 'polynomial in Z' and predicted a tie. But the high-r")
print("   properties are NOT monotonic in Z (atomic volume is the canonical PERIODIC")
print("   property -- Lothar Meyer 1870), so a poly-in-Z null would do poorly. The")
print("   correct null is the DEBYE MODEL, which gamma_phonon is a relabeling of.")
print(f"   Delta_r(Synchronism - Debye) = {r_v_invgamma:.4f} - {r_v_theta:.4f} "
      f"= {r_v_invgamma - r_v_theta:+.2e}  (ZERO, by definitional identity)")
print()
print("   Verdict: confirms S647 (self-correlation) with the exact mechanism;")
print("   corrects S651 (wrong null premise). The chemistry phonon-property network")
print("   re-derives the Debye model with theta_D relabeled as 2T/gamma_phonon.")
print("   Tension #3 answered: '89% validated against what?' -> against textbook")
print("   relations the coherence variables are definitional relabelings of.")
