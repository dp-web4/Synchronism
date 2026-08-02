import math
Om, OL = 0.315, 0.685
E   = lambda z: math.sqrt(Om*(1+z)**3 + OL)
dEdz= lambda z: 1.5*Om*(1+z)**2 / E(z)

# Ciocan et al. 2026: a0(z) = a0(0) + a1 z ; a1 = 1.59 +- 0.1 (95% CI => 1sigma = 0.051)
a1, a1_95 = 1.59, 0.10
a1_1s = a1_95/1.96

anchors = [("Ciocan+2026 fitted intercept", 1.00, 0.02),
           ("framework  cH0/2pi",           1.04, None),
           ("McGaugh+2016 SPARC (canon.)",  1.20, 0.26),
           ("Varasteanu+2025 MIGHTEE-HI",   1.69, 0.13)]

print("dE/dz|0 = 1.5*Om = %.4f     E(1) = %.4f     E(2) = %.4f" % (dEdz(0), E(1), E(2)))
print()
print("THE SLOPE RATIO I CALLED 'NORMALISATION-FREE' ON 08-01")
print("%-30s %8s %10s %10s" % ("anchor a0(0)", "a0(0)", "slope", "a1/slope"))
for n,a,_ in anchors:
    s = a*dEdz(0)
    print("%-30s %8.2f %10.3f %9.2fx" % (n, a, s, a1/s))
rs=[a1/(a*dEdz(0)) for _,a,_ in anchors]
print("  -> ratio spans %.2fx - %.2fx  (spread %.0f%% of the min; sign never flips)" %
      (min(rs), max(rs), 100*(max(rs)-min(rs))/min(rs)))
print()
print("THE z~1 SIGMA STATEMENT, with Ciocan's +-0.1 read as 95%% CI (1sigma = %.3f)" % a1_1s)
meas1 = 1.00 + a1*1.0   # Ciocan linear model at z=1
print("  Ciocan linear model at z=1: %.3f" % meas1)
print("%-30s %10s %10s %10s" % ("anchor", "brA(z=1)", "sig(meas)", "sig(both)"))
for n,a,ea in anchors:
    p = a*E(1)
    d = meas1 - p
    s_meas = d/ (a1_1s*1.0)
    tot = math.sqrt((a1_1s)**2 + (ea*E(1))**2) if ea else a1_1s
    print("%-30s %10.3f %+9.2f %+9.2f" % (n, p, s_meas, d/tot))
print()
print("z=2 GROWTH FACTOR (the LCDM-degeneracy comparison)")
print("  branch (A) E(2)              = %.3f" % E(2))
print("  Magneticum (Mayer+2023)      ~ 3    (1 sig fig, verified at source)")
print("  Ciocan linear extrapolation  = %.2f" % ((1.00+a1*2)/1.00))
