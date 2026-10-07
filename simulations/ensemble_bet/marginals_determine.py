"""Does the per-item answer distribution alone fix ensemble accuracy? (exact plurality calc by MC
over iid draws from each item's empirical answer distribution, no pairing info used)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); from analyze import load, arm_pack, vote_correct
import json
an = json.load(open(HERE + "/analysis.json"))["per_arm"]
arms, _, _ = load(HERE + "/answers_full.jsonl"); packs = arm_pack(arms)
rng = np.random.default_rng(1)
for k in sorted(packs):
    ans, cor = packs[k]; p_i = cor.mean(1); out = []
    for N in [1, 5, 12]:
        hits = 0; T = 400
        for it in range(ans.shape[0]):
            for _ in range(T):
                m = rng.integers(0, 12, N); hits += vote_correct(ans[it, m], cor[it, m])
        out.append((N, hits / (T * ans.shape[0]), an[str(k)]["acc"][str(N)]))
    ceil = np.mean([np.bincount(np.unique(ans[i].astype(np.int64), return_inverse=True)[1]).argmax()
                    == np.unique(ans[i].astype(np.int64), return_inverse=True)[1][cor[i].argmax()] if cor[i].any() else False
                    for i in range(ans.shape[0])])
    print(f"K={k:2d} " + "  ".join(f"N{N}: marg={a:.3f} reg={b:.3f}" for N, a, b in out) +
          f"  | plurality-modal-correct items={ceil:.3f}  gain1->12={out[2][2]-out[0][2]:+.3f}")
