#!/usr/bin/env python3
"""
Phase-9 — the MATTER sector: does the inflow model reproduce GR's *other* classic
test (perihelion precession), and what does that say about the instrument-effect bet?

Context
-------
The substrate arc (Phase-3c / 7 / 8) showed the framework's gravity-as-inflow model
reproduces GR *light* bending (4GM/c^2 b) via the eikonal swimmer  H = c|k| + k.u  with
u = -sqrt(2GM/r) r-hat (the Gullstrand-Painleve / "river" profile). EVERY phase traced
LIGHT. None ever put MATTER in an orbit. So "reproduces GR" was, strictly, "reproduces
GR's NULL sector." The other canonical GR test -- perihelion precession (a TIMELIKE-sector
effect) -- was never computed.

The faithful matter version of the framework's OWN rule (Galilean advection of the
medium-frame dispersion, exactly as light uses H=c|k|+k.u in absolute time) is:

    H = sqrt(m^2 c^4 + c^2 |p|^2) + p.u           (massive swimmer; m=0 recovers light)

This script:
  (A) Integrates a bound orbit under that swimmer Hamiltonian and measures precession.
  (B) Integrates a Schwarzschild geodesic (reference GR) and measures precession the
      same way.
  (C) Compares both to the weak-field GR formula  dphi = 6*pi*GM / (c^2 a (1-e^2)).
  (D) Checks the ALGEBRAIC claim: squaring the swimmer dispersion regenerates the exact
      Gullstrand-Painleve matter dispersion (so the agreement should be EXACT, not
      approximate). If precession instead comes out WRONG, that is the finding (the model
      would then reproduce GR kinematics/light only -- the analog-gravity limit).
  (E) Clock universality probe -- the instrument-effect tension: in this model every
      pattern is advected by the SAME u(r); dtau/dt for a static clock is sqrt(1-u^2/c^2)
      with NO internal-mechanism parameter. So all clocks dilate identically. That is
      exactly the clock UNIVERSALITY the "time-dilation-as-instrument-effect" bet needs to
      FAIL (its discriminator is mechanism-dependence).

Units: c = 1, m = 1 (test particle). GM set per-run.
Author: CBP-Claude (Opus 4.8), autonomous.
"""

import numpy as np
import json, os
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# (A) Framework swimmer:  H = sqrt(1 + p_r^2 + p_phi^2/r^2) - w(r) p_r,  w=sqrt(2GM/r)
#     polar canonical coords (r, phi, p_r, p_phi); p_phi conserved.
# ----------------------------------------------------------------------
def swimmer_deriv(state, GM):
    r, phi, pr, pphi = state
    w  = np.sqrt(2.0*GM/r)
    S  = np.sqrt(1.0 + pr*pr + (pphi*pphi)/(r*r))
    rdot   = pr/S - w
    phidot = (pphi/(r*r))/S
    prdot  = (pphi*pphi)/(r**3 * S) - pr*w/(2.0*r)
    return np.array([rdot, phidot, prdot, 0.0])

# ----------------------------------------------------------------------
# (B) Schwarzschild reference: standard binet orbit equation in phi,
#     d2u/dphi2 + u = GM/L^2 + 3 GM u^2   (c=1, u=1/r, L=specific ang. mom.)
#     Integrate in phi; measure perihelion-to-perihelion angle.
# ----------------------------------------------------------------------
def schwarzschild_orbit_precession(GM, L, u0, dphi=1e-4, n_orbits=6):
    # state: u, du/dphi ; start at perihelion (du/dphi=0, u=u0=1/r_min)
    u, up = u0, 0.0
    phi = 0.0
    peri_phis = [0.0]
    prev_u = u
    going_down = False  # after perihelion u decreases
    # integrate; record each perihelion (u local max -> r local min)
    u_hist = [u]; phi_hist=[0.0]
    steps = int(n_orbits * 2*np.pi / dphi * 1.6)
    for _ in range(steps):
        # RK4 on (u, up)
        def f(s):
            uu, uup = s
            return np.array([uup, GM/L**2 + 3*GM*uu*uu - uu])
        s = np.array([u, up])
        k1=f(s); k2=f(s+0.5*dphi*k1); k3=f(s+0.5*dphi*k2); k4=f(s+dphi*k3)
        s = s + (dphi/6.0)*(k1+2*k2+2*k3+k4)
        u, up = s
        phi += dphi
        u_hist.append(u); phi_hist.append(phi)
    u_hist=np.array(u_hist); phi_hist=np.array(phi_hist)
    # find perihelia = local maxima of u
    peri=[]
    for i in range(1,len(u_hist)-1):
        if u_hist[i]>u_hist[i-1] and u_hist[i]>=u_hist[i+1]:
            # parabolic refine
            denom=(u_hist[i-1]-2*u_hist[i]+u_hist[i+1])
            shift=0.5*(u_hist[i-1]-u_hist[i+1])/denom if denom!=0 else 0.0
            peri.append(phi_hist[i]+shift*dphi)
    if len(peri)<2: return None,None
    gaps=np.diff(peri)
    prec=np.mean(gaps)-2*np.pi
    # a,e from u extremes
    umin=np.min(u_hist); umax=np.max(u_hist)
    rmin=1/umax; rmax=1/umin
    a=0.5*(rmin+rmax); e=(rmax-rmin)/(rmax+rmin)
    return prec, (a,e)

def EL_from_turning_points(GM, r_p, r_a):
    """Conserved (E,L) for a bound orbit with turning points r_p<r_a.
    Turning-point condition rdot=0 of the swimmer reduces EXACTLY to the
    Schwarzschild one:  E^2 = (1-2GM/r)(1+L^2/r^2).  Solve both ends."""
    fp=1-2*GM/r_p; fa=1-2*GM/r_a
    L2=(fa-fp)/(fp/r_p**2 - fa/r_a**2)
    E2=fp*(1+L2/r_p**2)
    return np.sqrt(E2), np.sqrt(L2)

def integrate_swimmer(GM, r_peri, r_apo, n_orbits=6):
    # state = [r, phi, pr]; pphi=L conserved. Start AT perihelion turning point:
    # rdot=0 => pr = w*S, with S=E/(1-w^2).
    E,L=EL_from_turning_points(GM,r_peri,r_apo)
    pphi=L
    w0=np.sqrt(2*GM/r_peri); S0=E/(1-w0*w0); pr0=w0*S0
    def rhs(t,s):
        r,phi,pr=s
        w=np.sqrt(2.0*GM/r)
        S=np.sqrt(1.0+pr*pr+(pphi*pphi)/(r*r))
        rdot=pr/S - w
        phidot=(pphi/(r*r))/S
        prdot=(pphi*pphi)/(r**3*S) - pr*w/(2.0*r)
        return [rdot,phidot,prdot]
    def peri_event(t,s):
        r,phi,pr=s
        w=np.sqrt(2.0*GM/r); S=np.sqrt(1.0+pr*pr+(pphi*pphi)/(r*r))
        return pr/S - w
    peri_event.direction=1.0
    rc=0.5*(r_peri+r_apo)
    Sc=np.sqrt(1+(pphi/rc)**2); phidot_c=(pphi/rc**2)/Sc
    T_est=2*np.pi/abs(phidot_c)
    t_max=T_est*(n_orbits+2)
    sol=solve_ivp(rhs,[0,t_max],[r_peri,0.0,pr0],events=peri_event,
                  rtol=1e-12,atol=1e-13,max_step=T_est/200.0)
    ev=sol.y_events[0]
    if ev is None or len(ev)<2: return None,None,(E,L)
    peri_phi=ev[:,1]
    gaps=np.diff(peri_phi)
    prec=float(np.mean(gaps)-2*np.pi)
    rmin=float(np.min(sol.y[0])); rmax=float(np.max(sol.y[0]))
    a=0.5*(rmin+rmax); e=(rmax-rmin)/(rmax+rmin)
    return prec,(a,e),(E,L)

def weakfield_precession(GM,a,e):
    return 6*np.pi*GM/(a*(1-e*e))   # c=1

# ----------------------------------------------------------------------
# (D) algebraic identity check (symbolic-numeric): for random (r,pr,pphi),
#     compute H_swim = sqrt(1+pr^2+pphi^2/r^2) - w*pr, then verify
#     H^2 + 2 H w pr  ==  1 + (1-w^2) pr^2 + pphi^2/r^2   (GP matter dispersion)
# ----------------------------------------------------------------------
def identity_check(GM, ntrials=100000, seed=1):
    rng=np.random.default_rng(seed)
    r=rng.uniform(3*GM,500*GM,ntrials)     # outside horizon
    pr=rng.uniform(-2,2,ntrials)
    pphi=rng.uniform(-50,50,ntrials)
    w=np.sqrt(2*GM/r)
    S=np.sqrt(1+pr*pr+pphi*pphi/(r*r))
    H=S - w*pr
    lhs=H*H + 2*H*w*pr
    rhs=1 + (1-w*w)*pr*pr + pphi*pphi/(r*r)
    return float(np.max(np.abs(lhs-rhs)))

# ----------------------------------------------------------------------
# (E) clock universality: dtau/dt for a STATIC clock at r, in this model.
#     static => rdot=0 => pr=wS; the proper-time rate is the medium-frame
#     internal advance per coordinate tick. For GP/Schwarzschild static
#     observer dtau/dt = sqrt(1-2GM/r)=sqrt(1-w^2). Show it depends ONLY on
#     r (via universal u), with no internal-mechanism parameter.
# ----------------------------------------------------------------------
def clock_rate_static(GM, r):
    w=np.sqrt(2*GM/r)
    return np.sqrt(1-w*w)   # = sqrt(1-2GM/r), Schwarzschild static redshift

# ----------------------------------------------------------------------
def main():
    GM=1.0
    print("="*72)
    print("Phase-9: MATTER sector of the inflow-gravity model (precession)")
    print("="*72)

    # (D) identity first -- if this holds, precession MUST match GR exactly
    maxdev=identity_check(GM)
    print(f"\n(D) Swimmer dispersion vs exact GP/Schwarzschild matter dispersion")
    print(f"    max |H^2+2Hw p_r  -  (1+(1-w^2)p_r^2+p_phi^2/r^2)|  over 1e5 random states")
    print(f"    = {maxdev:.3e}   ->  {'ALGEBRAICALLY IDENTICAL' if maxdev<1e-10 else 'DIFFERS'}")

    # (A,B,C) precession across a few field strengths, parametrized by turning points
    print(f"\n(A/B/C) Perihelion precession (GM=1, c=1); orbit set by turning points")
    print(f"  {'r_peri':>7} {'r_apo':>7} | {'a':>8} {'e':>6} | {'swimmer Δφ':>13} "
          f"{'Schw Δφ':>13} {'6πGM/a(1-e²)':>14} | {'swim/Schw':>9}")
    rows=[]
    for r_peri,r_apo in [(120,180),(120,240),(200,300),(60,120),(40,90)]:
        ps,(a,e),(E,L)=integrate_swimmer(GM,r_peri,r_apo,n_orbits=5)
        u0=1.0/r_peri  # exact perihelion
        pg,ae2=schwarzschild_orbit_precession(GM,L,u0,n_orbits=5)
        wf=weakfield_precession(GM,a,e)
        ratio=ps/pg if pg else float('nan')
        print(f"  {r_peri:>7} {r_apo:>7} | {a:>8.2f} {e:>6.3f} | {ps:>13.6f} "
              f"{(pg if pg else float('nan')):>13.6f} {wf:>14.6f} | {ratio:>9.6f}")
        rows.append(dict(r_peri=r_peri,r_apo=r_apo,a=a,e=e,E=E,L=L,
                         swim_prec=ps,schw_prec=pg,weakfield=wf,swim_over_schw=ratio))

    # (E) clock universality
    print(f"\n(E) Static-clock rate dtau/dt = sqrt(1-2GM/r)  (NO mechanism parameter)")
    for r in [10,20,50,100]:
        print(f"    r={r:>4}:  dtau/dt = {clock_rate_static(GM,r):.6f}   "
              f"(identical for ANY internal clock mechanism -> UNIVERSAL)")

    print("\nInterpretation")
    print("-"*72)
    print("If (D) is identical and (A)≈(B)≈(C): the swimmer rule reproduces GR's TIMELIKE")
    print("sector EXACTLY, not just the null sector -- because it IS Gullstrand-Painleve")
    print("(=Schwarzschild) in disguise, for matter as well as light. The model therefore")
    print("agrees with GR with ZERO wiggle room wherever u=sqrt(2GM/r). A novel (Bucket-0)")
    print("prediction can ONLY live where u DEPARTS from the GP profile (Phase-8 Faith-B")
    print("saturation at small r) or where DISCRETENESS breaks the continuum identity")
    print("(B7 Umklapp / LIV). And (E) shows the inflow model PREDICTS clock universality,")
    print("which is exactly what the 'time-dilation-as-instrument-effect' bet needs to FAIL.")

    out=dict(GM=GM, dispersion_max_dev=maxdev, rows=rows,
             clock_rates={int(r):clock_rate_static(GM,r) for r in [10,20,50,100]})
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/phase9_matter_sector_precession_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/phase9_matter_sector_precession_result.json")

if __name__=="__main__":
    main()
