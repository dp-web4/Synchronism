"""Publisher 2026-11-18: is the open-questions "virial cap" (x <= ~0.02 for every bound
system) a property of the density law, or of a length convention?

Checks back-annotation 5d2d684c (ii): "the virial-ratio 'knee unreachable in every sector'
argument holds only for A*V^2 with the same length in A and in the smoothing of rho".

Part 1 (identity controls, algebra):
  two-sided convention  rho_bar(l) = 3 V_c(l)^2 / (4 pi G l^2),  rho_crit = A(l) V^2,
  A(l) = 4 pi / (beta^2 G l^2)  ->  x = (3/16 pi^2) beta^2 (V_c/V)^2   (paper L45: 0.019 at beta=1)
  primary form (S53 L58-59, no 4 pi): A = 1/(G alpha^2 l^2) -> x = 3/(4 pi) alpha^2 (V_c/V)^2 (paper: 0.24)
  Jeans A at l = 8 kpc, 4 pi form: 4.57e-5 (paper conclusion L78 / exec L190: 4.56e-5)
Part 2 (data): S53's executed law rho_crit = 0.028 V^0.5 (V in km/s, M_sun/pc^3), i.e. the same
  V^2/(G alpha^2 R_half^2) with R_half = R0 V^0.75, against the midplane density the archive's
  fits use, rho = Sigma / 2h, h = 0.3 kpc. Sigma = 0.5 SBdisk + 0.7 SBbul (stars only; gas
  omitted, so Sigma and every x here are LOWER bounds). Count discs whose innermost point has x > 1.
  V = mean Vobs over the last 3 points (a V_flat proxy, not Lelli's Vflat column).
"""
import math, sys

G = 4.30091e-3  # pc (km/s)^2 / M_sun

# ---- Part 1: identity controls
x2 = 3 / (16 * math.pi**2)
xp = 3 / (4 * math.pi)
AJ = 4 * math.pi / (G * (8000.0) ** 2)
print(f"two-sided 4pi-form cap at beta=1: {x2:.4f}  (paper 0.019)")
print(f"two-sided primary-form cap at alpha=1: {xp:.4f}  (paper 0.24)")
print(f"Jeans A at 8 kpc, 4pi form: {AJ:.3e}  (paper 4.56e-5)")
ok = abs(x2 - 0.019) < 5e-4 and abs(xp - 0.24) < 5e-3 and abs(AJ - 4.56e-5) / 4.56e-5 < 0.01
print("IDENTITY CONTROLS", "PASS" if ok else "FAIL")
if not ok:
    sys.exit(1)

# S53 coefficient check: A = 1/(G alpha^2 R0^2), R0 = 0.088 kpc/(km/s)^0.75, alpha = 1.1
A53 = 1 / (G * 1.1**2 * (88.0) ** 2)
print(f"S53 A = 1/(G a^2 R0^2) at alpha=1.1, R0=0.088 kpc: {A53:.4f}  (S53 prints ~0.028)")

# ---- Part 2: SPARC
path = sys.argv[1]
gal = {}
for line in open(path, encoding="latin-1"):
    p = line.split()
    if len(p) != 10:
        continue
    try:
        vals = list(map(float, p[1:]))
    except ValueError:
        continue
    gal.setdefault(p[0], []).append(vals)

h = 300.0  # pc
n = cross = 0
xs = []
for g, rows in gal.items():
    rows.sort(key=lambda r: r[1])
    if len(rows) < 4:
        continue
    V = sum(r[2] for r in rows[-3:]) / 3
    r0 = rows[0]
    sigma = 0.5 * r0[7] + 0.7 * r0[8]
    rho = sigma / (2 * h)
    rc = 0.028 * V**0.5
    x = rho / rc
    xs.append((x, g, V))
    n += 1
    cross += x > 1
xs.sort()
print(f"galaxies: {n}; innermost-point x > 1 under 0.028 V^0.5 with rho = Sigma_*/2h: {cross}")
print("median x:", f"{xs[len(xs)//2][0]:.3f}", " max:", f"{xs[-1][0]:.1f} ({xs[-1][1]})")
for name in ("NGC2403", "NGC3198", "NGC7331", "DDO154", "UGC00128"):
    for x, g, V in xs:
        if g == name:
            print(f"  {g}: V={V:.0f} km/s x_inner={x:.2f}")
# the same discs under the two-sided virial convention can never exceed ~0.24 (primary) at V_c<=V
