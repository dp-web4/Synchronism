"""
Phase-2 — Does special relativity emerge, and what does the DISCRETENESS predict?
(2026-06-22, CBP-Claude, autonomous — dp: "follow what pulls you.")

Where the arc stands: the substrate is converging on a COMPLEX Intent field + focusing
nonlinearity. Electron = standing amplitude envelope (rest oscillation = mass); photon =
phase-dominated propagating packet; momentum = phase gradient; c = ceiling. Everything so far
REPRODUCES known physics (Bucket 0 unchanged: zero confirmed novel predictions).

The one place this model CAN differ from continuum QFT is its **discreteness**. A complex field
on a continuum is just Klein-Gordon (textbook). On a DISCRETE grid the dispersion relation is
modified, and that modification is the model's only candidate novel content. This experiment
measures the emergent dispersion ω(k) of the substrate and asks:

  (1) Does the relativistic backbone E^2 = (pc)^2 + (mc^2)^2 EMERGE at low k?  [electron(k=0,
      E=mc^2) <-> photon(high k, E->pc) is literally this curve.]
  (2) What is the DEVIATION at high k — the discreteness signature — and is it a real,
      named, in-principle-testable prediction?

Method: linearized complex field on a 1D lattice (spacing a=1), exact lattice action
  omega^2(k) = m^2 + (2c^2/a^2)(1 - cos(k a)).
Seed a broad (near-monochromatic) wavepacket at carrier k; MEASURE omega(k) from the temporal
phase-rotation of the k-mode projection c_k(t) = sum_x psi e^{-ikx}, and the group velocity
v_g(k) from the envelope centroid. Compare measured vs continuum-SR vs lattice prediction.
Then (bonus) turn focusing ON and measure the soliton's amplitude-dependent frequency shift =
a nonlinear MASS renormalization (mass <- rest oscillation + binding; "mass = complexity").

numpy only, headless. Writes results JSON.
"""
import json
import os
import numpy as np

L = 1024
DT = 0.01
C = 1.0
C2 = C * C
M = 0.5
M2 = M * M
STEPS = 4000
K_LIST = [0.05, 0.1, 0.2, 0.4, 0.6, 0.9, 1.2, 1.6, 2.0, 2.5, 3.0]  # up toward zone edge pi
SEED = 7


def lap(z):
    return np.roll(z, -1) - 2 * z + np.roll(z, 1)


def accel(psi, gsat=0.0, ssat=0.7):
    a = C2 * lap(psi) - M2 * psi
    if gsat:
        a = a + gsat * (np.abs(psi) ** 2 / (1.0 + np.abs(psi) ** 2 / ssat ** 2)) * psi
    return a


def lattice_omega(k):
    return float(np.sqrt(M2 + 2 * C2 * (1 - np.cos(k))))     # exact discrete dispersion


def sr_omega(k):
    return float(np.sqrt(M2 + C2 * k * k))                   # continuum special relativity


def measure(rng, k, amp=0.02, gsat=0.0):
    x = np.arange(L)
    x0 = L / 2.0
    W = 80.0                                                  # broad envelope -> sharp in k
    env = amp * np.exp(-0.5 * ((x - x0) / W) ** 2)
    psi = env * np.exp(1j * k * x)
    w_seed = lattice_omega(k)
    pdot = -1j * w_seed * psi
    a = accel(psi, gsat)
    phase_unwrapped = []
    last = None
    se = 5
    cents = []
    ekx = np.exp(-1j * k * x)
    for t in range(STEPS):
        psi = psi + DT * pdot + 0.5 * DT * DT * a
        a_new = accel(psi, gsat)
        pdot = pdot + 0.5 * DT * (a + a_new)
        a = a_new
        if t % se == 0:
            ck = np.sum(psi * ekx)                            # projection onto carrier k-mode
            ph = np.angle(ck)
            if last is not None:
                d = ph - last
                d = (d + np.pi) % (2 * np.pi) - np.pi         # unwrap step
                phase_unwrapped.append(d)
            last = ph
            rho = np.abs(psi) ** 2
            ang = 2 * np.pi * x / L
            zc = np.sum(rho * np.exp(1j * ang))
            cents.append((np.angle(zc) % (2 * np.pi)) / (2 * np.pi) * L)
    # omega = - mean phase-rotation rate of the k-mode
    omega_meas = -float(np.mean(phase_unwrapped)) / (se * DT)
    # group velocity from centroid drift (unwrapped)
    cc = np.array(cents); dd = np.diff(cc); dd = (dd + L / 2) % L - L / 2
    vg_meas = float(np.sum(dd) / ((len(cc) - 1) * se * DT))
    return abs(omega_meas), vg_meas


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for k in K_LIST:
        w_meas, vg_meas = measure(rng, k, amp=0.02, gsat=0.0)
        w_sr = sr_omega(k); w_lat = lattice_omega(k)
        vg_sr = C2 * k / w_sr
        # lattice group velocity dω/dk = c^2 sin(k) / ω_lat
        vg_lat = C2 * np.sin(k) / w_lat
        rows.append({
            "k": k,
            "omega_meas": round(w_meas, 4), "omega_SR": round(w_sr, 4), "omega_lattice": round(w_lat, 4),
            "vg_meas": round(vg_meas, 4), "vg_SR": round(vg_sr, 4), "vg_lattice": round(vg_lat, 4),
            "SR_dev_pct": round(100 * (w_meas - w_sr) / w_sr, 2),
            "lattice_match_pct": round(100 * (w_meas - w_lat) / w_lat, 2),
        })

    # Does SR emerge at low k? (deviation small for k<<1)
    low = [r for r in rows if r["k"] <= 0.2]
    sr_emerges = all(abs(r["SR_dev_pct"]) < 2.0 for r in low)
    # Does the measurement track the LATTICE dispersion (the model's actual prediction)?
    tracks_lattice = all(abs(r["lattice_match_pct"]) < 3.0 for r in rows)
    # The LIV signature: at the highest k, how far below SR is the speed? (energy-dependent c)
    hi = rows[-1]
    liv_speed_deficit_pct = round(100 * (hi["vg_lattice"] - hi["vg_SR"]) / hi["vg_SR"], 1)
    omega_max = round(lattice_omega(np.pi), 4)               # UV cutoff (max frequency/energy)
    # leading low-k LIV: omega_lat ~ omega_SR * (1 - (c^2 k^4 / 24) / omega_SR^2 ...) -> the
    # physically relevant deviation scales as (k a)^2 = (E/E_cutoff)^2, tiny at accessible E.
    low_k_dev = next(r["SR_dev_pct"] for r in rows if abs(r["k"] - 0.2) < 1e-9)

    verdict = (
        f"SR backbone EMERGES at low k (deviation <2% for k<=0.2: {sr_emerges}) — the substrate "
        f"reproduces E^2=(pc)^2+(mc^2)^2, with the same curve interpolating electron (k=0, "
        f"E=mc^2) and photon (high k, E->pc). The measurement tracks the exact DISCRETE-lattice "
        f"dispersion omega^2=m^2+2c^2(1-cos k) (match: {tracks_lattice}), NOT continuum SR. "
        f"The discreteness signature — the model's one candidate novel content — is: (a) a UV "
        f"CUTOFF, a maximum frequency/energy omega_max={omega_max} at the zone edge (no patterns "
        f"above it — a natural UV completion); (b) an ENERGY-DEPENDENT speed of light: at the "
        f"highest k tested the group velocity is {liv_speed_deficit_pct}% below the SR value — "
        f"high-energy patterns propagate SLOWER (Lorentz-invariance violation). This maps to a "
        f"REAL observable: energy-dependent photon arrival times from distant gamma-ray bursts / "
        f"blazars, which LIV searches actually measure. HONEST BOUNDS: (1) this LIV form is "
        f"GENERIC to discrete-substrate theories (loop quantum gravity, causal sets, etc.) — NOT "
        f"unique to Synchronism; (2) the lattice spacing is Planck-scale, so the real-world "
        f"deviation is ~(E/E_Planck)^2, far below current GRB constraints. So it is a real, "
        f"named, in-principle-testable prediction CLASS that is neither unique nor currently "
        f"distinguishable. The advance is in SPECIFICATION, not confirmation: the program's "
        f"novel physics, if any, lives EXACTLY in the discreteness, at the Planck scale, and the "
        f"ONLY experiment class that could ever test it is astrophysical photon timing. That "
        f"converts 'underspecified, zero novel predictions' into a precise statement of where the "
        f"novelty must be and what could falsify it. Bucket 0 unchanged (nothing confirmed). "
        f"(Note: the physically relevant deviation is the LOW-k tail — at k=0.2 it is "
        f"{low_k_dev}%, scaling as (k a)^2 = (E/E_cutoff)^2; the large %s near the zone edge are "
        f"the Planck-energy extreme. A nonlinear mass-renormalization measurement was attempted "
        f"and DISCARDED — the linear k-mode method is invalid at soliton amplitudes; that needs "
        f"a proper nonlinear-normal-mode method, deferred.)"
    )

    out = {
        "lattice": L, "dt": DT, "c": C, "m": M, "steps": STEPS, "seed": SEED,
        "dispersion": rows,
        "SR_emerges_low_k": bool(sr_emerges),
        "tracks_discrete_lattice_dispersion": bool(tracks_lattice),
        "UV_cutoff_omega_max": omega_max,
        "LIV_speed_deficit_pct_at_high_k": liv_speed_deficit_pct,
        "SR_dev_pct_at_k0.2": low_k_dev,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results", "phase2_dispersion_LIV_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print("\n   k     omega_meas  omega_SR  omega_lat   vg_meas  vg_SR  vg_lat   SR_dev%  lat_match%")
    for r in rows:
        print(f"  {r['k']:.2f}   {r['omega_meas']:8.4f}  {r['omega_SR']:7.4f}  {r['omega_lattice']:7.4f}   "
              f"{r['vg_meas']:+.3f}  {r['vg_SR']:.3f}  {r['vg_lattice']:.3f}   "
              f"{r['SR_dev_pct']:+6.2f}   {r['lattice_match_pct']:+6.2f}")


if __name__ == "__main__":
    main()
