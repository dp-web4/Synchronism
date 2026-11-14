"""Publisher 2026-11-08: where does the bulges-in profiled-a0 Cassini z stop converging?
cf84bb8c says gamma=0.3 does not converge (8-22 sigma) and gamma>=0.489 does; gamma=0.4 is unstated,
and the paper's 11-07 note prints "every gamma from 0.3 to 3". Read-only. Same a0 fit, same qumond_q, same Q2 datum.
Prediction before running: 0.3 spreads by several sigma; 0.489+ by <=0.05; 0.4 unknown.
"""
import sys
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, '/home/dp/ai-workspace/synchronism-site/explorer/scripts')
sys.path.insert(0, '/home/dp/ai-workspace/Synchronism/simulations')
import mond_delta_family_cassini_at_sparc_profiled_a0 as M  # noqa: E402
from sparc_cassini_q2 import TanhLogNu, qumond_q, q2_si  # noqa: E402

gid, gobs, gbar = M.load()


def fit(nu):
    f = lambda la: np.sum((np.log10(gobs) - np.log10(gbar * nu(gbar / 10**la)))**2)
    return 10**minimize_scalar(f, bounds=(-12, -9), method='bounded', options={'xatol': 1e-9}).x


print(f"{'gamma':>6} {'g_ext':>9} {'g_ext/a0':>9} {'z_base':>8} {'min':>8} {'max':>8} {'spread':>7}  draws")
for gm in (0.3, 0.35, 0.4, 0.45, 0.489, 0.6, 1.0, 2.0, 3.0):
    nu = TanhLogNu(gm)
    a0 = fit(nu)
    for ge in (2.00e-10, 2.32e-10, 2.48e-10):
        zs = []
        for da in (0.999, 1.0, 1.001, 1.003):
            for ao, ro in ((128, 512), (192, 768), (256, 1024)):
                a = a0 * da
                try:
                    res = qumond_q(nu, ge / a, angular_order=ao, radial_order=ro, log_v_bound=14.0)
                    zs.append((q2_si(res.q, a) - M.MEAN) / M.SIG)
                except Exception as e:
                    zs.append(float('nan'))
        zs = np.array(zs)
        print(f"{gm:6.3f} {ge:9.2e} {ge/a0:9.2f} {zs[3]:8.2f} {np.nanmin(zs):8.2f} {np.nanmax(zs):8.2f} "
              f"{np.nanmax(zs)-np.nanmin(zs):7.2f}  {' '.join(f'{v:.1f}' for v in zs)}", flush=True)
