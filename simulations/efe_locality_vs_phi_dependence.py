"""
EFE under a non-local C versus a Phi-dependent C, and the DF2 arithmetic.
Publisher lane, 2026-08-24.  Supports the [AMENDED 2026-08-24] note in
whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md

Question: does making C(rho) NON-LOCAL invalidate the EFE = 0 derivation?
Answer:   no.  The derivation turns on C being independent of Phi, not on C
          being local.  A fully non-local but Phi-independent C leaves the
          field equation linear in Phi, so superposition holds and EFE = 0
          exactly.  Only Phi-dependence (a grad-Phi keyed C) produces an EFE.

Measured here:  non-local fixed C  -> EFE = 5.6e-13  (solver precision)
                grad-Phi keyed C   -> EFE = 4.6e-02  (11 orders larger)

Part 2 recomputes the NGC 1052-DF2 numbers the repo circulated on 2026-08-23.
"""
import math
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

# 2D grid; solve  div[ C(x) grad Phi ] = 4 pi G rho   with Dirichlet Phi=0 on boundary.
N=161; L=40.0; h=L/(N-1)
x=np.linspace(-L/2,L/2,N); X,Y=np.meshgrid(x,x,indexing='ij')
def gauss(x0,y0,s,amp): return amp*np.exp(-((X-x0)**2+(Y-y0)**2)/(2*s*s))

rho_int = gauss(0,0,1.5,1.0)          # the dwarf, at origin
rho_ext = gauss(12.0,0,2.5,6.0)       # the host, 12 units away

def build(C):
    # face-centred conductivities, standard 5-point conservative discretisation
    idx=lambda i,j: i*N+j
    Cx=0.5*(C[1:,:]+C[:-1,:]); Cy=0.5*(C[:,1:]+C[:,:-1])
    rows=[];cols=[];vals=[]
    for i in range(N):
        for j in range(N):
            k=idx(i,j)
            if i==0 or j==0 or i==N-1 or j==N-1:
                rows.append(k);cols.append(k);vals.append(1.0); continue
            d=0.0
            for (di,dj,cf) in ((-1,0,Cx[i-1,j]),(1,0,Cx[i,j]),(0,-1,Cy[i,j-1]),(0,1,Cy[i,j])):
                rows.append(k);cols.append(idx(i+di,j+dj));vals.append(cf/h**2); d-=cf/h**2
            rows.append(k);cols.append(k);vals.append(d)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N*N,N*N))

def solve(A,f):
    b=(4*np.pi*f).ravel().copy()
    mask=np.zeros((N,N),bool); mask[0,:]=mask[-1,:]=mask[:,0]=mask[:,-1]=True
    b[mask.ravel()]=0.0
    return spl.spsolve(A,b).reshape(N,N)

def internal_field(Phi):
    # relative acceleration inside the dwarf: grad Phi minus its value at the dwarf centre
    gx,gy=np.gradient(Phi,h)
    c=N//2
    core=(np.abs(X)<4)&(np.abs(Y)<4)
    return np.stack([(gx-gx[c,c])[core],(gy-gy[c,c])[core]])

print("=== CASE A: C fixed, NON-LOCAL in origin (formation-history field), independent of Phi ===")
# a C field that is NOT a function of local rho at all: it 'remembers' a compact past state
C_hist = 0.6*np.exp(-((X)**2+(Y)**2)/(2*3.0**2)) + 0.05
A=build(C_hist)
Pi=solve(A,rho_int); Pe=solve(A,rho_ext); Pt=solve(A,rho_int+rho_ext)
print(f"  superposition residual |Phi_tot-(Phi_int+Phi_ext)|_max = {np.abs(Pt-(Pi+Pe)).max():.3e}")
fi=internal_field(Pi); ft=internal_field(Pt-Pe)
print(f"  internal relative field, isolated vs host-present: max diff = {np.abs(ft-fi).max():.3e}"
      f"  (scale {np.abs(fi).max():.3e})  -> EFE = {np.abs(ft-fi).max()/np.abs(fi).max():.2e}")

print()
print("=== CASE B: C depends on |grad Phi| (Phi-dependent, AQUAL-like) ===")
def solve_nl(f):
    Phi=np.zeros((N,N))
    for it in range(60):
        gx,gy=np.gradient(Phi,h); g=np.sqrt(gx**2+gy**2)
        C=np.tanh(0.5*np.log1p(g/0.02))+0.05
        Pn=solve(build(C),f)
        d=np.abs(Pn-Phi).max(); Phi=0.5*Phi+0.5*Pn
        if d<1e-10: break
    return Phi
Pi2=solve_nl(rho_int); Pe2=solve_nl(rho_ext); Pt2=solve_nl(rho_int+rho_ext)
fi2=internal_field(Pi2); ft2=internal_field(Pt2-Pe2)
print(f"  superposition residual = {np.abs(Pt2-(Pi2+Pe2)).max():.3e}")
print(f"  internal relative field, isolated vs host-present: max diff = {np.abs(ft2-fi2).max():.3e}"
      f"  (scale {np.abs(fi2).max():.3e})  -> EFE = {np.abs(ft2-fi2).max()/np.abs(fi2).max():.2e}")


print()
print("=== PART 2: the DF2 arithmetic ===")
G = 4.300917270e-6  # kpc (km/s)^2 / Msun

# NGC 1052-DF2 baryonic parameters (van Dokkum+2018; Danieli+2019)
M_star = 2.0e8      # Msun
R_e    = 2.2        # kpc  (half-light)
sig_obs = 8.5       # km/s (Danieli+2019)

# Wolf+2010 mass estimator inverted: M(<r_1/2) = 3 sigma^2 r_1/2 / G,  r_1/2 = 4/3 R_e
r_half = (4.0/3.0)*R_e
# Newtonian (baryons only) dispersion at r_half, half the stellar mass enclosed
for fenc in (0.5, 1.0):
    sig_N = math.sqrt(G*fenc*M_star/(3*r_half))
    B_req = (sig_obs/sig_N)**2
    print(f"f_enc={fenc}: sigma_N={sig_N:.2f} km/s   B_req=(obs/N)^2={B_req:.2f}   C_req=1/B={1/B_req:.3f}")

print()
# What does strict local C ~ 0.04 predict?  g = g_N/C  =>  sigma = sigma_N/sqrt(C)
for C in (0.04,):
    for fenc in (0.5, 1.0):
        sig_N = math.sqrt(G*fenc*M_star/(3*r_half))
        print(f"C={C}: predicted sigma = {sig_N/math.sqrt(C):.1f} km/s  (obs {sig_obs}) -> factor {sig_N/math.sqrt(C)/sig_obs:.1f}x")
print()
# What sigma_N would be needed for the draft's "sigma ~ 80" at C=0.04?
print(f"draft's sigma~80 at C=0.04 requires sigma_N = {80*math.sqrt(0.04):.1f} km/s")
print(f"draft's sigma~80 would need C = {(1/(80/6.0))**2:.4f} if sigma_N=6.0 ... i.e. C=(sig_N/sig)^2")
for sN in (6.0, 7.0, 8.0):
    print(f"   sigma_N={sN}: C needed for sigma=80 is {(sN/80)**2:.5f}")
print()
# C_formation 0.5-0.7 -> predicted sigma
for C in (0.5, 0.6, 0.7):
    for fenc in (0.5,):
        sig_N = math.sqrt(G*fenc*M_star/(3*r_half))
        print(f"C_formation={C}: sigma = {sig_N/math.sqrt(C):.2f} km/s (obs {sig_obs}), ratio {sig_N/math.sqrt(C)/sig_obs:.2f}")
