"""Is K1's null a property of (b,c), or of the equation's background stability?
Chartered equation: dA/dt = A + (1+ib) lap(A) - (1+ic)|A|^2 A   (1D periodic)
Linearise about A=0:  mode k grows at Re[1 - (1+ib)k^2] = 1 - k^2.
=> every mode with |k|<1 grows; k=0 grows at rate 1. The zero background is
   LINEARLY UNSTABLE for ALL (b,c). A localized pulse sits on a growing background.
Test it directly at four (b,c) from the arc's own grid, from a narrow Gaussian IC.
"""
import numpy as np
L, N = 256, 256
x = np.arange(N); dx = 1.0
def lap(A): return (np.roll(A,1)+np.roll(A,-1)-2*A)/dx**2
def run(b, c, T=2000, dt=0.02):
    A = (0.5*np.exp(-((x-N/2)**2)/(2*4.0**2))).astype(complex)
    A += 1e-6*(np.random.default_rng(0).standard_normal(N)+0j)
    def rhs(A): return A + (1+1j*b)*lap(A) - (1+1j*c)*np.abs(A)**2*A
    for _ in range(int(T/dt)):
        k1=rhs(A); k2=rhs(A+dt/2*k1); k3=rhs(A+dt/2*k2); k4=rhs(A+dt*k3)
        A = A + dt/6*(k1+2*k2+2*k3+k4)
        if not np.isfinite(A).all(): return None
    return A
def width(A):
    p = np.abs(A)**2
    if p.sum() <= 0: return np.nan
    return float((p>0.1*p.max()).sum())          # support above 10% of peak
print("%-6s %-6s %-14s %-14s %-14s" % ("b","c","|A| far-field","|A| peak","support(px)"))
for b,c in [(0.0,0.0),(1.5,1.5),(3.0,0.0),(7.0,7.0)]:
    A = run(b,c)
    if A is None: print("%-6s %-6s  blew up" % (b,c)); continue
    far = float(np.abs(A[:20]).mean())           # far from the initial pulse
    print("%-6s %-6s %-14.4f %-14.4f %-14.0f" % (b,c,far,float(np.abs(A).max()),width(A)))
print()
print("Prediction if the background is unstable: far-field |A| -> 1 (saturation of")
print("the cubic term), support -> full lattice (256), at EVERY (b,c). No localization")
print("is available anywhere in the chartered family.")
