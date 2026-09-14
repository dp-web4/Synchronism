"""Is the Session-100 dark-energy sector a Cardassian model?  (maintainer, 2026-09-14)

Framework:  H^2 = (8 pi G / 3) rho_m / C(rho_m),  C = tanh(gamma ln(1+x)),  x = rho_m/rho_crit
            rho_DE = rho_m (1-C)/C
Cardassian (Freese & Lewis 2002):   H^2 = (8 pi G/3) g(rho_m)          [general class]
MP Cardassian (Gondolo & Freese):   g = rho_m [1 + (rho_car/rho_m)^(q(1-n))]^(1/q)

Checks
 1. identity  rho_m/(rho_m+rho_DE) == C  (so C(a) is the model's Omega_m(a))
 2. gamma = 1/2: 1/C == 1 + 2/x exactly  -> MP Cardassian q=1, n=0, rho_car=2 rho_crit (= LCDM)
 3. x >> 1: 1/C -> 1 + 2 x^(-2 gamma)     -> MP Cardassian q=1, n=1-2gamma, rho_car=2^(1/2gamma) rho_crit
 4. x << 1: rho_DE -> rho_crit/gamma const -> w -> -1 (NOT the MP-Cardassian low-density limit unless gamma=1/2)
 5. where today sits: C(x0) = Omega_m ; and w(x) from continuity, w = -1 - (1/3) dln rho_DE/dln a
"""
import numpy as np

Om = 0.315
def C(x, g): return np.tanh(g*np.log1p(x))
def rho_de(x, g): return x*2/np.expm1(2*g*np.log1p(x))   # = x(1-C)/C exactly, no 1-C underflow; units of rho_crit

print("1. identity rho_m/(rho_m+rho_DE) - C, max |diff| over x in [1e-3,1e4], gamma in {0.3,0.487,0.5,2}:")
xs = np.logspace(-3, 4, 400)
for g in (0.3, 0.487, 0.5, 2.0):
    d = np.max(np.abs(xs/(xs+rho_de(xs, g)) - C(xs, g)))
    print(f"   gamma={g:<5}  {d:.2e}")

print("2. gamma=1/2: max |1/C - (1+2/x)| =", f"{np.max(np.abs(1/C(xs,0.5)-(1+2/xs))):.2e}")

print("3. high-density asymptote, relative error of 1+2x^(-2g) vs 1/C:")
for g in (0.3, 0.487, 2.0):
    for x in (1, 10, 100, 1e3):
        ex = 1/C(x, g); asy = 1 + 2*x**(-2*g)
        print(f"   gamma={g:<5} x={x:<6g} 1/C={ex:.5f} asym={asy:.5f} rel.err={abs(asy-ex)/ex:.2e}")

print("4. low-density limit rho_DE(x->0) * gamma (expect -> 1):")
for g in (0.3, 0.487, 2.0):
    print(f"   gamma={g:<5} x=1e-4: {rho_de(1e-4,g)*g:.5f}")

print("5. today's x0 (C(x0)=Omega_m) and w(x) = -1 + (dln rho_DE/dln x)   [rho_m ~ a^-3 => dln/dln a = -3 dln/dln x]")
from math import atanh, expm1
for g in (0.3, 0.487, 0.5, 2.0):
    x0 = expm1(atanh(Om)/g)
    lx = np.log(xs); lr = np.log(rho_de(xs, g))
    w = -1 + np.gradient(lr, lx)
    w0 = np.interp(np.log(x0), lx, w)
    print(f"   gamma={g:<5} x0={x0:.3f}  w(x0)={w0:+.4f}  w(x=1e4)={w[-1]:+.4f} (asym -2g={-2*g:+.3f})  w(x=1e-3)={w[0]:+.4f}")
    # w monotone? (no -1 crossing)
    print(f"             w monotone in x: {bool(np.all(np.diff(w)<=1e-9) or np.all(np.diff(w)>=-1e-9))}, crosses -1: {bool(np.any(np.sign(w[1:]+1)!=np.sign(w[:-1]+1)) and g!=0.5)}")
