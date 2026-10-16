"""Does count compensate in an INTERACTING ensemble? Sign of d(post-interaction acc)/d(ln N) vs discrimination.
Same model as positive_control.py (uniform verification V, false-acceptance F), N varied. (CBP-Claude, 2026-10-16)"""
import numpy as np, positive_control as pc
out = []
for level in [0.05, 0.1, 0.2, 0.4]:
    accs, votes = [], []
    for N in [3, 5, 9, 15, 25]:
        pc.N = N
        pc.skill = np.full((N, pc.K), 0.35)
        for a in range(N): pc.skill[a, a % pc.K] = 0.80
        post, vote = [], []
        for rep in range(4):
            cats = pc.rng.integers(0, pc.K, 200); truth, A = pc.answer_matrix(200, cats)
            V = np.full((N, N, pc.K), level); rr = np.random.default_rng(rep)
            vote.append(pc.accuracy(A, truth)); post.append(pc.accuracy(pc.interact(A, truth, cats, V, rr), truth))
        accs.append(np.mean(post)); votes.append(np.mean(vote))
    slope_int = np.polyfit(np.log([3, 5, 9, 15, 25]), accs, 1)[0]
    slope_vote = np.polyfit(np.log([3, 5, 9, 15, 25]), votes, 1)[0]
    print(f"V={level:.2f}  post acc by N(3,5,9,15,25)={np.round(accs,3).tolist()}  slope/lnN interact={slope_int:+.3f}  vote={slope_vote:+.3f}")
