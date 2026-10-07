"""Is the declared <C> proxy an algebraic identity? (CBP-Claude, 2026-10-07)
1) Synthetic: build ensembles with KNOWN error coupling, run the declared proxy.
2) Real data: compute a non-degenerate redundancy measure (ICC of per-item accuracy)
   and an out-of-arm-difficulty variant of the declared proxy."""
import sys, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
sys.path.insert(0, HERE)
from analyze import compatibility, load, arm_pack

rng = np.random.default_rng(7)
I, P = 48, 12
def synth(rho, q=0.55):
    # latent per-item shared factor with weight rho -> correlated correctness
    z_item = rng.normal(size=(I, 1)); z_own = rng.normal(size=(I, P))
    # persona clusters: two blocks with anti-aligned item preferences (complementary errors)
    lat = np.sqrt(rho) * z_item + np.sqrt(1 - rho) * z_own
    thr = np.quantile(lat, 1 - q)
    return (lat > thr).astype(float)
def clustered(q=0.55):
    # 2 persona clusters that fail on DISJOINT item halves: strongly complementary
    cor = np.ones((I, P))
    for a in range(P):
        bad = np.arange(I//2) if a < P//2 else np.arange(I//2, I)
        cor[rng.choice(bad, size=int(len(bad)*0.9), replace=False), a] = 0
    return cor
print("SYNTHETIC — declared proxy vs known coupling")
for rho in [0.0, 0.3, 0.6, 0.9, 0.99]:
    c = synth(rho); icc = np.var(c.mean(1)) / (c.mean()*(1-c.mean()))
    print(f"  shared-factor rho={rho:4.2f}: q={c.mean():.3f}  proxy<C>={compatibility(c):.4f}  ICC={icc:.3f}")
c = clustered(); icc = np.var(c.mean(1)) / (c.mean()*(1-c.mean()))
print(f"  two disjoint-failure clusters: q={c.mean():.3f}  proxy<C>={compatibility(c):.4f}  ICC={icc:.3f}")
print(f"  identity value 1+1/(P-1) = {1+1/(P-1):.4f}")

print("\nREAL DATA (answers_full.jsonl)")
arms, _, _ = load(HERE + "/answers_full.jsonl")
packs = arm_pack(arms)
allcor = {k: packs[k][1] for k in packs}
for k in sorted(packs):
    c = allcor[k]; q = c.mean(); p_i = c.mean(1)
    icc = (np.var(p_i) - q*(1-q)/P) / (q*(1-q))  # binomial-noise-corrected ICC
    # out-of-arm difficulty: item accuracy from the OTHER six arms
    other = np.mean([allcor[j].mean(1) for j in allcor if j != k], axis=0)[:, None]
    r = c - other
    cs = [np.corrcoef(r[:, a], r[:, b])[0, 1] for a in range(P) for b in range(a+1, P)
          if r[:, a].std() > 1e-9 and r[:, b].std() > 1e-9]
    maj_ok = np.mean(p_i > 0.5)
    print(f"  K={k:2d} q={q:.3f} declared<C>={compatibility(c):.4f}  ICC={icc:.3f}  "
          f"out-of-arm<C>={1-np.mean(cs):.3f}  frac items p_i>0.5={maj_ok:.3f}")
