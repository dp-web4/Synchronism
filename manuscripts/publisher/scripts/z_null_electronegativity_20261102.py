"""Test S651's premise: is Pauling electronegativity 'near-monotonic in Z' so that a smooth
monotonic / polynomial-in-Z null reaches r >= 0.95 'by construction'? Pauling values (CRC), Z=1..54,
noble gases without a Pauling value omitted. Typed from standard tables, not loaded from a database."""
import numpy as np
EN = {1:2.20,3:0.98,4:1.57,5:2.04,6:2.55,7:3.04,8:3.44,9:3.98,11:0.93,12:1.31,13:1.61,14:1.90,15:2.19,
      16:2.58,17:3.16,19:0.82,20:1.00,21:1.36,22:1.54,23:1.63,24:1.66,25:1.55,26:1.83,27:1.88,28:1.91,
      29:1.90,30:1.65,31:1.81,32:2.01,33:2.18,34:2.55,35:2.96,36:3.00,37:0.82,38:0.95,39:1.22,40:1.33,
      41:1.60,42:2.16,43:1.90,44:2.20,45:2.28,46:2.20,47:1.93,48:1.69,49:1.78,50:1.96,51:2.05,52:2.10,
      53:2.66,54:2.60}
Z = np.array(sorted(EN)); y = np.array([EN[z] for z in Z])
from scipy.stats import spearmanr
print(f"n={len(Z)}  Spearman(EN, Z) = {spearmanr(Z, y).correlation:+.3f}")
for deg in (1, 2, 3, 5):
    fit = np.polyval(np.polyfit(Z, y, deg), Z)
    print(f"polynomial-in-Z degree {deg}: r = {np.corrcoef(fit, y)[0,1]:.3f}")
# best monotonic function of Z = isotonic regression (upper bound for ANY monotone null)
def iso(y):
    blocks = [[v, 1] for v in y]
    out = []
    for b in blocks:
        out.append(b)
        while len(out) > 1 and out[-2][0] > out[-1][0]:
            v2, n2 = out.pop(); v1, n1 = out.pop(); out.append([(v1*n1+v2*n2)/(n1+n2), n1+n2])
    return np.concatenate([[v]*n for v, n in out])
for name, yy in (("increasing", y), ("decreasing", -y)):
    f = iso(yy); print(f"best monotone-{name} fit of Z: r = {np.corrcoef(f, yy)[0,1]:.3f}")
