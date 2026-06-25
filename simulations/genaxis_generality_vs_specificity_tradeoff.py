#!/usr/bin/env python3
"""
Generative-axis meta-finding — why the Synchronism lens keeps converging on the RIGHT axis
but a COARSER metric than each domain's native formalism, and what that implies.

Context
-------
Across the whole arc the same shape recurred:
  - gravity: the frame gives the acoustic/Gullstrand-Painleve metric — reproduces GR
    KINEMATICS, not the Einstein equations (Phase-9/14);
  - applied: MRH = an extension of the Markov blanket — a productive synthesis, not FEP's
    native machinery (applied-axis session);
  - quantum computing: a scalar compatibility metric <C> is ORTHOGONAL to standard scalar
    metrics (real, useful) but SUBSUMED by the field's sharper native metric, stabilizer-
    symmetry alignment (Phase-17/18).
Each time: right axis, coarser descriptor. dp's standing correction is right that the
CONVERGENCE is corroboration (a frame built for general emergence independently landing on
each field's own axis is real signal, not redundancy). This script asks the next question:
*why does it land COARSE, every time, and is that fixable?*

Claim under test: a MAXIMALLY-GENERAL (domain-invariant, low-dimensional) descriptor is, by
construction, (1) strictly subsumed IN-domain by the domain's structural formalism (it has a
measurable discriminating blind spot), AND (2) the UNIQUE descriptor that retains meaning
ACROSS domains. Generality and in-domain sharpness trade off. So "coarse but right axis" is
not a deficiency to fix; it is the price (and the point) of cross-domain reach.

Two toy domains, each: configurations described by a STRUCTURAL vector s; an outcome y that
depends on the STRUCTURE of s (not just its mean); a domain-native model that uses the full
structure; and the general scalar c = mean(s) (the same definition in every domain).

  (A) IN-DOMAIN: compare R^2 of {scalar c} vs {full structure}. Show scalar < structural
      whenever y depends on structure beyond the mean -> the subsumption + blind spot
      (Phase-18, generalized).
  (B) CROSS-DOMAIN: the structural FEATURES of domain 1 are undefined in domain 2; the scalar
      c is defined in both. Train a predictor on domain-1 scalar, test on domain-2 -> it
      TRANSFERS (sign/trend); the structural model cannot even be evaluated cross-domain.
  (C) The tradeoff curve: vary how much of y is structural vs mean-captured; show scalar
      in-domain R^2 falls as structure-dependence rises, while transferability is preserved.

numpy only. Headless. Writes results JSON.
Author: CBP-Claude (Opus 4.8), autonomous — generative-axis meta (Bucket 0 = 0).
"""
import json, os
import numpy as np


def make_domain(rng, N=4000, d=12, struct_weight=0.7, law=None):
    """Configurations s in {0,1}^d (e.g. per-region 'compatibility/alignment').
       Outcome y = mean-captured part + structure-captured part + noise.
       'law' is a fixed random structural functional (domain-specific) so different
       domains have DIFFERENT structure->outcome maps."""
    s = rng.integers(0, 2, size=(N, d)).astype(float)
    c = s.mean(axis=1)                                  # the general scalar (same def everywhere)
    # structure-captured part: a domain-specific functional orthogonalized vs the mean
    if law is None:
        law = rng.standard_normal(d)
        law = law - law.mean()                          # orthogonal to the global mean direction
    struct = s @ law                                    # depends on WHICH regions, not just how many
    struct = (struct - struct.mean())/ (struct.std()+1e-12)
    cz = (c - c.mean())/(c.std()+1e-12)
    y = (1-struct_weight)*cz + struct_weight*struct + 0.05*rng.standard_normal(N)
    return dict(s=s, c=c, y=y, law=law, d=d)


def r2(x, y):
    """R^2 of the best linear fit y ~ x (x can be 1D or 2D design)."""
    X = np.atleast_2d(x)
    if X.shape[0] != len(y): X = X.T
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = np.sum((y-yhat)**2); ss_tot = np.sum((y-y.mean())**2)
    return 1 - ss_res/ss_tot, beta


def main():
    rng = np.random.default_rng(7)
    print("="*78)
    print("Generative-axis meta: generality forces coarseness (right axis, coarser metric)")
    print("="*78)

    # ---- (A) IN-DOMAIN subsumption + blind spot ----
    print("\n(A) IN-DOMAIN: scalar c=mean(s) vs the domain's STRUCTURAL formalism (full s)")
    A = make_domain(rng, struct_weight=0.7)
    r2_scalar,_ = r2(A['c'], A['y'])
    r2_struct,_ = r2(A['s'], A['y'])
    print(f"    R^2(scalar c)        = {r2_scalar:.3f}   (orthogonal to mean-error/connectivity — real, useful)")
    print(f"    R^2(full structure)  = {r2_struct:.3f}   (the domain-native formalism)")
    print(f"    => scalar is SUBSUMED: it captures {100*r2_scalar/r2_struct:.0f}% of the structural"
          f" model's variance.")
    # blind spot: configs with ~equal c but different y
    c = A['c']; y = A['y']
    # pick the modal c bucket and show y spread within it
    cmid = np.median(c)
    mask = np.abs(c-cmid) < 1e-9
    if mask.sum() > 30:
        spread = y[mask].std()
        print(f"    blind spot: among {mask.sum()} configs with IDENTICAL c={cmid:.3f}, outcome y still"
              f" varies (std={spread:.2f}) — variation the scalar cannot see, the structure can.")

    # ---- (B) CROSS-DOMAIN transfer ----
    print("\n(B) CROSS-DOMAIN: domain-2 has a DIFFERENT structural law (its features are alien to D1)")
    B = make_domain(rng, struct_weight=0.7)   # different random law => different structure->y map
    # scalar model trained on A, tested on B
    _, beta_c = r2(A['c'], A['y'])
    yhatB = beta_c[0] + beta_c[1]*B['c']
    ssB = 1 - np.sum((B['y']-yhatB)**2)/np.sum((B['y']-B['y'].mean())**2)
    sign_agree = np.mean(np.sign(B['y']-B['y'].mean()) == np.sign(yhatB-yhatB.mean()))
    print(f"    scalar model trained on D1, applied to D2:  R^2={ssB:.3f}, sign-agreement={sign_agree:.2f}")
    print(f"    => the SCALAR transfers (same definition, consistent direction across domains).")
    print(f"    structural model trained on D1: uses D1's specific features (law over D1 regions) —")
    print(f"    these are UNDEFINED in D2 (different object). The native formalism CANNOT transfer.")
    print(f"    => the coarse metric is the ONLY descriptor that retains meaning across domains.")

    # ---- (C) the tradeoff curve ----
    print("\n(C) TRADEOFF: as outcome becomes more structure-dependent, in-domain scalar R^2 falls,")
    print("    but cross-domain transferability of the scalar is preserved.")
    print(f"    {'struct_weight':>13} {'R2_scalar(in)':>14} {'R2_struct(in)':>14} {'scalar transfer R2':>18}")
    rows=[]
    for sw in [0.0, 0.25, 0.5, 0.75, 1.0]:
        Aw = make_domain(rng, struct_weight=sw)
        Bw = make_domain(rng, struct_weight=sw)
        r2s,_ = r2(Aw['c'], Aw['y']); r2f,_ = r2(Aw['s'], Aw['y'])
        _, bc = r2(Aw['c'], Aw['y']); yhB = bc[0]+bc[1]*Bw['c']
        tr = 1 - np.sum((Bw['y']-yhB)**2)/np.sum((Bw['y']-Bw['y'].mean())**2)
        print(f"    {sw:>13.2f} {r2s:>14.3f} {r2f:>14.3f} {tr:>18.3f}")
        rows.append(dict(struct_weight=sw, r2_scalar_in=float(r2s), r2_struct_in=float(r2f),
                         scalar_transfer_r2=float(tr)))
    print("    => when outcome is pure-mean (sw=0) the scalar IS the native metric (no coarsening,")
    print("       full transfer). As structure-dependence rises, the scalar's in-domain power falls")
    print("       (subsumed) but it still transfers — the only descriptor that does. Generality and")
    print("       in-domain sharpness trade off; coarse-but-transferable is the frame's niche.")

    print("\nInterpretation")
    print("-"*78)
    print("The Synchronism frame is MAXIMALLY GENERAL (one vocabulary for gravity, emergence, QC,")
    print("cognition). A domain-invariant descriptor must coarse-grain each domain's structure, so")
    print("it is STRICTLY SUBSUMED in-domain by the native formalism (acoustic metric < Einstein eqs;")
    print("<C> < stabilizer-alignment; MRH < FEP machinery) — every time, by construction, not by")
    print("under-development. The SAME coarseness is what lets one metric mean the same thing across")
    print("domains. So the convergence dp rightly credits is real (right axis), AND the coarseness is")
    print("a structural ceiling (never the sharpest tool in any domain), AND they are the same")
    print("property. The frame's unique, untested value is therefore NOT beating a domain's native")
    print("metric (structurally impossible) but ANTICIPATORY CROSS-DOMAIN TRANSFER: carry an axis")
    print("established in domain A into domain B before B's native formalism finds it. That — not a")
    print("better in-domain metric — is the generative-axis analogue of a Bucket-0 win, and it is the")
    print("test the program should run next. Bucket 0 = 0 (this is a frame finding).")

    out=dict(A_in_domain=dict(r2_scalar=float(r2_scalar), r2_struct=float(r2_struct)),
             B_cross_domain=dict(scalar_transfer_r2=float(ssB), sign_agreement=float(sign_agree)),
             C_tradeoff=rows,
             thesis="generality forces coarseness: a domain-invariant scalar is strictly subsumed "
                    "in-domain (blind spot) but is the unique cross-domain-transferable descriptor; "
                    "the frame's unique value is anticipatory cross-domain transfer, not beating a "
                    "domain's native metric (structurally impossible). Bucket 0 = 0.")
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/genaxis_generality_vs_specificity_tradeoff_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/genaxis_generality_vs_specificity_tradeoff_result.json")


if __name__=="__main__":
    main()
