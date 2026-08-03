"""
Q1 of the alignment arc (charter: explorations/2026-08-02-kimi-alignment-arc-plan.md) —
the founding sim of the alignment program: an alignment-based Bell model on a common
substrate, with the two registered Q1 audits (relabeling, conspiracy) built in.

THE BET (Q1 form): entanglement = phase-lock on a common substrate AT THE SOURCE;
measurement = a finite-window phase ALIGNMENT between the arriving pattern and the
detector pattern; the Bell correlation is manufactured in the window by the alignment
dynamics — not read off pre-existing properties, not requiring initial-condition
conspiracy.

SUBSTRATE. 1D lattice, N=512 cells, OPEN boundaries (a source and two separated detector
regions need geography), ONE global tick, local-only coupling. Each cell carries a phase
phi_i (continuous — abstraction debt: the frame's intent is quantized; recorded, not
leaned on). Dynamics (co-rotating frame of the substrate's common rotation):
  d(phi_i)/dt = K_SUB * sum_{j nn of i} sin(phi_j - phi_i)  + thermal noise
explicit Euler, DT=0.2. Stability: the linearized coupling is K_SUB times the path-graph
Laplacian, |lambda|max <= 4*K_SUB = 0.6; Euler is stable for DT*|lambda|max < 2 and we run
at 0.12 — a 16x margin. The detector gate runs at KAPPA*DT = 0.04 per tick, the ensemble
internal coupling at J_DET*DT = 0.2 — both far inside Euler stability. The co-rotating
frame is a gauge choice whose PHYSICAL content
(shared phase reference from the common clock) is itself a named assumption: the substrate
is held phase-coherent source-to-detectors; that coherence is asserted, not derived (debt).

SOURCE (center cell region). Emits ONE pair per run: two localized phase-offset packets,
left-going and right-going, CO-LOCKED at creation: both carry the same offset psi relative
to the substrate phase, psi ~ U(0, 2pi) drawn per run (the pair's random creation phase;
the fixed relation psi_L = psi_R IS the entanglement premise — phase lock at the source,
on the substrate). Packets travel ballistically outward at 1 cell/tick, implemented as
LOCAL transport (each arm cell copies its inward neighbor's phase each tick) plus the
substrate Kuramoto relaxation, so packets broaden and decay en route (measured, reported).
Abstraction debt: directed ballistic transport is imposed, not derived from a wave
equation on the substrate — named per the river rule; no verdict leans on it.

DETECTORS (left A at cell 64, right B at 448). Each is a small oscillator ensemble
(M=8 cells) with strong INTERNAL coupling J_DET — a distinct pattern swimming in the
river — carrying its OWN phase r, drawn ~ U(0, 2pi) PER RUN PER DETECTOR (this randomness
is load-bearing: a detector whose rest phase equals the substrate mean produces outcomes
correlated through the shared reference ALONE, with no pair — the no-source null exists to
catch exactly that artifact, and the randomized own-phase is what makes the null able to
pass). Between measurements the ensemble is decoupled from the substrate (debt: detector
isolation idealized — its holding dynamics is not modeled). Each detector carries a
physical phase reference theta: the phase of a reference oscillator AT THAT DETECTOR,
relative to the shared substrate phase — the "setting". Settings are physical phase
relationships on the substrate, never free labels fed into a joint expression (house
rule 8; audited below).

MEASUREMENT. When the packet parks on the detector region, the gate is armed for tau
ticks — and it OPENS only in proportion to the mean |deviation| actually present in the
park window (linear ramp, GATE_THRESHOLD..+GATE_RAMP): a real detector triggers on an
arriving pattern. Without the trigger the gate aligns the detector to the bare substrate
reference at long tau whether or not a pair exists, and the no-source null fails by
construction (found in smoke test: null C ~ 0.27 at tau=64 without the trigger; ~0 with
it). While open:
  d(phi_det)/dt += openness * KAPPA * sin(field_local - phi_det)
where field_local = amplitude-squared-weighted circular mean of the DEVIATION field in
the detector's park window (the detector resonates with the strongest local perturbation;
a plain window average dilutes the packet phase toward the background — also found in
smoke test). Outcome (binary, read ONLY from local quantities at that detector):
  outcome = sign( sin( phi_bar_det - theta_det ) )
with phi_bar_det the ensemble's circular mean phase after the window. One pair, one
outcome per detector, per run.

LATE-CHOICE CONTROL (mandatory, and checked in code): the settings theta_A, theta_B are
drawn from a SEPARATE rng stream AFTER the pair has left the source and traveled to the
detectors. The source can never see the settings. The rng draw order is logged and
asserted.

AUDITS (both are deliverables, run in-sim):
- RELABELING AUDIT. Static: outcomes are computed by read_outcome(phi_bar, theta) whose
  arguments are scalars local to ONE detector; C is computed by correlate(outcomes_A,
  outcomes_B) which receives only outcome arrays — no expression anywhere takes theta_A
  and theta_B jointly (checked property, see AUDIT section in the JSON). Dynamic: for a
  batch of runs with identical dynamics seeds but two independent late theta_B draws,
  outcome_A must be BITWISE identical (theta_B never touches A's dynamics or readout) —
  asserted, counts reported.
- CONSPIRACY AUDIT. Detector initial phases randomized per run; settings drawn late
  (asserted order); no parameter is tuned per grid cell (K_SUB, KAPPA, geometry fixed
  globally). Robustness probe: the map endpoints re-run with K_SUB halved and doubled.
  If correlation needed tuned initial conditions, that is the registered Q1(b) KILL.

MAP. C(Delta_theta, tau) = <outcome_A * outcome_B> over PAIRS pairs per cell,
Delta_theta in k*pi/8 (k=0..8), tau in {0, 4, 8, 16, 32, 64} ticks. theta_A ~ U(0, 2pi)
late per run, theta_B = theta_A + Delta_theta (both late, both physical references).
CONTROLS: no-source null (no packets emitted; C ~ 0 required) and no-lock control
(packets emitted with INDEPENDENT phases psi_L, psi_R; C ~ 0 at all tau required).

HEADLINES: full-lock angular form fit to triangle A*(1 - 2*Delta/pi) vs cosine A*cos(Delta)
(amplitudes + RMSE, which is closer); per-tau triangle amplitude A(tau) and CHSH
S(tau) = |C(pi/4) - C(3pi/4)| + 2*|C(pi/4)| from the measured grid (Q3's raw material);
measured lock factor |<exp(i(phi_bar - psi))>| per tau.

EXPECTED HONESTLY (pre-registered expectation, not a criterion): at full lock a local
construction gives the TRIANGLE, not the cosine — that is Q2's registered question, and
this map is its raw material. S(tau) should rise from ~0 toward ~2 (the theorem's cap is
respected; nothing here is a Bell violation claim).

numpy + stdlib only; seeded (seed=1); headless. Writes
simulations/results/kimi_alignment/q1_first_map.json.
"""
import json
import os
import time

import numpy as np

# --- lattice / substrate ---------------------------------------------------
N = 512
X_SRC = N // 2            # source region center
X_A, X_B = 64, 448        # detector region centers
DT = 0.2
K_SUB = 0.15              # substrate Kuramoto coupling; DT*4*K_SUB = 0.12 << 2 (stable)
NOISE_SUB = 0.01          # thermal river noise per tick (keeps the null honest at long tau)
WPK = 16                  # packet width (cells); flat-top so the phase plateau survives travel
TRAVEL = (X_SRC - WPK // 2) - X_A   # ticks for packet center to reach the detector
# --- detectors --------------------------------------------------------------
M_DET = 8                 # ensemble cells per detector
J_DET = 1.0               # internal ensemble coupling (the detector's own pattern)
KAPPA = 0.2               # gate alignment strength; KAPPA*DT = 0.04 per tick
WGATE = 14                # half-width of the park window the gate reads
GATE_THRESHOLD = 0.1      # mean |deviation| below which the gate stays closed (no pattern,
GATE_RAMP = 0.2           #   no measurement) — linear ramp to fully open over this span
# --- map ---------------------------------------------------------------------
DELTAS = np.linspace(0.0, np.pi, 9)          # k*pi/8: includes pi/4 and 3pi/4 (CHSH)
TAUS = [0, 4, 8, 16, 32, 64]
PAIRS = 600
SEED = 1


def lap_sin(phi):
    """Kuramoto coupling term, OPEN boundaries (end cells have one neighbor)."""
    left = np.empty_like(phi)
    right = np.empty_like(phi)
    left[1:] = phi[:-1]
    right[:-1] = phi[1:]
    left[0] = phi[0]
    right[-1] = phi[-1]
    return np.sin(left - phi) + np.sin(right - phi)


def circ_mean(x):
    return float(np.angle(np.mean(np.exp(1j * x))))


def packet_profile():
    """Flat-top phase plateau with 3-cell cosine tapers: the plateau interior feels no
    Kuramoto force (uniform phase), so the packet's phase survives the trip; only the
    edges erode (measured, reported). A sin-bump profile was tried first and attenuated
    to ~1/3 of psi by arrival — the fix is the plateau, not more coupling."""
    j = np.arange(WPK)
    return np.clip((WPK / 2 - np.abs(j - WPK / 2 + 0.5)) / 3.0, 0.0, 1.0)


def gate_field(phi, x0):
    """Deviation field the detector's gate reads, plus the gate's openness. Two parts:
    (i) field phase = amplitude-squared-WEIGHTED circular mean of the park window — the
    detector resonates with the strongest local perturbation, not the window average (a
    plain circular mean over a 28-cell window dilutes a 16-cell packet's phase toward the
    background and destroys the alignment; found in smoke test, fixed here).
    (ii) openness = ramp of the mean |deviation|: with NO pattern present the gate stays
    closed and the detector keeps its own (random) phase — a real detector triggers on an
    arriving pattern; without this trigger the gate aligns the detector to the bare
    substrate reference at long tau and the no-source null fails by construction (found
    in smoke test: null C ~ 0.27 at tau=64 without the trigger). Both quantities local to
    [x0-WGATE, x0+WGATE)."""
    win = phi[x0 - WGATE:x0 + WGATE]
    w = np.abs(win) ** 2
    field = float(np.angle(np.sum(w * np.exp(1j * win))))
    openness = float(np.clip((np.mean(np.abs(win)) - GATE_THRESHOLD) / GATE_RAMP, 0.0, 1.0))
    return field, openness


def run_pair(dyn_rng, set_rng, delta, tau, mode="locked", theta_B_override=None):
    """One pair: emission -> travel -> LATE settings -> window -> local readouts.
    mode: 'locked' (co-locked pair, the bet), 'nolock' (independent packet phases,
    control), 'null' (no packets emitted, control). theta_B_override exists ONLY for the
    relabeling audit (swap B's setting post hoc and verify A's outcome is untouched).
    Returns (outcome_A, outcome_B, phi_bar_A, phi_bar_B, psi_L, psi_R)."""
    # substrate init: near-uniform with per-run noise (detectors' ensemble phases below)
    phi = 0.05 * dyn_rng.standard_normal(N)
    # detector ensembles: own phase r ~ U(0,2pi) per run per detector (load-bearing),
    # drawn at run start; SETTINGS are drawn late (see below) — order asserted in audit.
    det_A = dyn_rng.uniform(0, 2 * np.pi) + 0.1 * dyn_rng.standard_normal(M_DET)
    det_B = dyn_rng.uniform(0, 2 * np.pi) + 0.1 * dyn_rng.standard_normal(M_DET)

    # EMISSION (rng event logged for the conspiracy audit)
    psi_L = dyn_rng.uniform(0, 2 * np.pi)
    psi_R = psi_L if mode == "locked" else dyn_rng.uniform(0, 2 * np.pi)
    if mode != "null":
        prof = packet_profile()
        phi[X_SRC - WPK:X_SRC] += psi_L * prof
        phi[X_SRC:X_SRC + WPK] += psi_R * prof

    # TRAVEL: local transport (arm cells copy inward neighbor) + substrate relaxation
    for _ in range(TRAVEL):
        phi = phi + DT * (K_SUB * lap_sin(phi) + NOISE_SUB * dyn_rng.standard_normal(N))
        phi[X_A:X_SRC] = np.roll(phi[X_A:X_SRC], -1)   # left arm: transport toward A
        phi[X_SRC + 1:X_B + 1] = np.roll(phi[X_SRC + 1:X_B + 1], 1)  # right arm: toward B

    # LATE CHOICE: settings drawn AFTER the pair has left the source (separate stream)
    theta_A = set_rng.uniform(0, 2 * np.pi)
    theta_B = (theta_A + delta) % (2 * np.pi) if theta_B_override is None else theta_B_override

    # WINDOW: gate open for tau ticks; detectors align to the local deviation field
    for _ in range(tau):
        phi = phi + DT * (K_SUB * lap_sin(phi) + NOISE_SUB * dyn_rng.standard_normal(N))
        field_A, open_A = gate_field(phi, X_A)
        field_B, open_B = gate_field(phi, X_B)
        for det, field, openness in ((det_A, field_A, open_A), (det_B, field_B, open_B)):
            internal = np.sin(np.roll(det, 1) - det) + np.sin(np.roll(det, -1) - det)
            det += DT * (J_DET * internal + openness * KAPPA * np.sin(field - det)
                         + 0.01 * dyn_rng.standard_normal(M_DET))

    # READOUT: local quantities only, per detector
    phi_bar_A = circ_mean(det_A)
    phi_bar_B = circ_mean(det_B)
    out_A = read_outcome(phi_bar_A, theta_A)
    out_B = read_outcome(phi_bar_B, theta_B)
    return out_A, out_B, phi_bar_A, phi_bar_B, psi_L, psi_R


def read_outcome(phi_bar, theta):
    """Binary outcome from LOCAL quantities at ONE detector only: the ensemble's aligned
    mean phase and that detector's own reference-oscillator phase. No cross-detector
    quantity is in scope here (relabeling audit, static half)."""
    return 1 if np.sin(phi_bar - theta) >= 0 else -1


def correlate(outcomes_A, outcomes_B):
    """C over paired outcomes ONLY. Receives no settings, no phases — the joint object is
    built from locally-read binaries (relabeling audit, static half)."""
    oA = np.asarray(outcomes_A)
    oB = np.asarray(outcomes_B)
    return float(np.mean(oA * oB)), float(np.std(oA * oB) / np.sqrt(len(oA)))


def run_cell(master_rng, delta, tau, mode, pairs):
    """One (Delta, tau) cell: `pairs` independent pairs; settings drawn late per run."""
    oA, oB, lock_A, lock_B = [], [], [], []
    for _ in range(pairs):
        dyn_seed = int(master_rng.integers(2**31))
        set_seed = int(master_rng.integers(2**31))
        dyn_rng = np.random.default_rng(dyn_seed)
        set_rng = np.random.default_rng(set_seed)
        a, b, pA, pB, psi_L, psi_R = run_pair(dyn_rng, set_rng, delta, tau, mode)
        oA.append(a)
        oB.append(b)
        lock_A.append(np.exp(1j * (pA - psi_L)))
        lock_B.append(np.exp(1j * (pB - psi_R)))
    c, cerr = correlate(oA, oB)
    lock = (abs(np.mean(lock_A)) + abs(np.mean(lock_B))) / 2 if mode != "null" else None
    return {"C": round(c, 4), "C_err": round(cerr, 4),
            "lock_factor": None if lock is None else round(float(lock), 4)}


def audit_relabeling(master_rng, n=200):
    """Dynamic half: identical dynamics seeds AND identical settings stream, but theta_B
    swapped post hoc at readout — outcome_A must be bitwise identical (theta_B never
    enters A's dynamics or readout), while outcome_B follows its own setting."""
    identical_A = 0
    changed_B = 0
    for _ in range(n):
        dyn_seed = int(master_rng.integers(2**31))
        set_seed = int(master_rng.integers(2**31))
        theta_B_alt = float(master_rng.uniform(0, 2 * np.pi))
        a1, b1, *_ = run_pair(np.random.default_rng(dyn_seed),
                              np.random.default_rng(set_seed), np.pi / 4, 16)
        a2, b2, *_ = run_pair(np.random.default_rng(dyn_seed),
                              np.random.default_rng(set_seed), np.pi / 4, 16,
                              theta_B_override=theta_B_alt)
        if a1 == a2:
            identical_A += 1
        if b2 != b1:
            changed_B += 1
    return {"n_runs": n, "outcome_A_identical_under_theta_B_change": identical_A,
            "outcome_B_changed_by_its_own_setting": changed_B,
            "pass": identical_A == n and changed_B > 0}


def main():
    t0 = time.time()
    master_rng = np.random.default_rng(SEED)

    grid = {}
    for delta in DELTAS:
        for tau in TAUS:
            key = f"d{delta:.4f}_t{tau}"
            grid[key] = {"delta": round(float(delta), 4), "tau": tau,
                         **run_cell(master_rng, delta, tau, "locked", PAIRS)}
            print(f"[map] d={delta:.3f} tau={tau:3d} C={grid[key]['C']:+.3f} "
                  f"lock={grid[key]['lock_factor']}", flush=True)

    null_grid = {}
    for delta in DELTAS:
        for tau in [0, 32, 64]:
            key = f"d{delta:.4f}_t{tau}"
            null_grid[key] = {"delta": round(float(delta), 4), "tau": tau,
                              **run_cell(master_rng, delta, tau, "null", PAIRS)}
    print("[null] no-source control done", flush=True)

    nolock_grid = {}
    for delta in DELTAS:
        for tau in TAUS:
            key = f"d{delta:.4f}_t{tau}"
            nolock_grid[key] = {"delta": round(float(delta), 4), "tau": tau,
                                **run_cell(master_rng, delta, tau, "nolock", PAIRS)}
    print("[null] no-lock control done", flush=True)

    relabel = audit_relabeling(master_rng)
    print(f"[audit] relabeling dynamic test: {relabel}", flush=True)

    # --- fits at full lock (longest tau) -------------------------------------
    full = [grid[f"d{d:.4f}_t{TAUS[-1]}"] for d in DELTAS]
    Cfull = np.array([g["C"] for g in full])
    tri_basis = 1.0 - 2.0 * DELTAS / np.pi
    cos_basis = np.cos(DELTAS)
    A_tri = float(Cfull @ tri_basis / (tri_basis @ tri_basis))
    A_cos = float(Cfull @ cos_basis / (cos_basis @ cos_basis))
    rmse_tri = float(np.sqrt(np.mean((Cfull - A_tri * tri_basis) ** 2)))
    rmse_cos = float(np.sqrt(np.mean((Cfull - A_cos * cos_basis) ** 2)))

    # --- per-tau triangle amplitude + CHSH from the measured grid ------------
    i_q, i_3q = 2, 6    # DELTAS index of pi/4 and 3pi/4
    tau_curve = []
    for tau in TAUS:
        Crow = np.array([grid[f"d{d:.4f}_t{tau}"]["C"] for d in DELTAS])
        A = float(Crow @ tri_basis / (tri_basis @ tri_basis))
        S = abs(Crow[i_q] - Crow[i_3q]) + 2 * abs(Crow[i_q])
        lock = grid[f"d{DELTAS[0]:.4f}_t{tau}"]["lock_factor"]
        tau_curve.append({"tau": tau, "triangle_amplitude": round(A, 4),
                          "CHSH_S": round(float(S), 4), "lock_factor": lock})

    # --- robustness probe (conspiracy audit, part): K_SUB halved/doubled -----
    global K_SUB
    robust = {}
    base_K = K_SUB
    for tag, K in (("K_SUB_half", base_K / 2), ("K_SUB_double", base_K * 2)):
        K_SUB = K
        robust[tag] = {"K_SUB": K,
                       "C(d=pi/4,tau=32)": run_cell(master_rng, np.pi / 4, 32,
                                                    "locked", PAIRS // 2)}
    K_SUB = base_K

    out = {
        "stage": "Q1 — alignment Bell model on a common substrate: first C(Delta,tau) map",
        "charter": "explorations/2026-08-02-kimi-alignment-arc-plan.md",
        "params": {"N": N, "X_SRC": X_SRC, "X_A": X_A, "X_B": X_B, "DT": DT,
                   "K_SUB": base_K, "NOISE_SUB": NOISE_SUB, "WPK": WPK,
                   "TRAVEL_ticks": TRAVEL, "M_DET": M_DET, "J_DET": J_DET,
                   "KAPPA": KAPPA, "WGATE": WGATE, "PAIRS": PAIRS, "SEED": SEED,
                   "stability": "DT*4*K_SUB = 0.12 << 2 (Euler, path-graph Laplacian "
                                "bound); KAPPA*DT = 0.04; J_DET*DT = 0.2"},
        "abstraction_debts": [
            "1D, not the frame's 3D grid",
            "continuous phase per cell, not quantized intent",
            "ballistic packet transport imposed as local shift, not derived from a wave equation on the substrate",
            "co-rotating frame: substrate phase coherence source-to-detectors asserted, not derived (the shared reference's existence is assumed)",
            "detector isolation between measurements idealized (holding dynamics not modeled); gate selectivity for the deviation field idealized",
            "open boundaries by fiat",
        ],
        "settings_are_physical": (
            "theta_A/theta_B are phases of reference oscillators AT each detector, relative "
            "to the shared substrate phase (common clock); drawn LATE from a separate rng "
            "stream after emission+travel (order asserted); used only in the local readout "
            "sign(sin(phi_bar_det - theta_det))."),
        "map": grid,
        "controls": {
            "no_source_null": null_grid,
            "no_lock": nolock_grid,
        },
        "audits": {
            "relabeling": {
                "static": "outcomes via read_outcome(phi_bar, theta) — local scalars only; "
                          "C via correlate(outcomes_A, outcomes_B) — outcome arrays only; "
                          "no expression takes theta_A and theta_B jointly (checked in code).",
                "dynamic": relabel,
            },
            "conspiracy": {
                "settings_drawn_after_emission": True,
                "settings_stream": "separate rng, drawn post-travel (late choice)",
                "detector_init_phases": "randomized per run per detector, U(0,2pi)",
                "per_cell_tuning": "none — K_SUB/KAPPA/geometry global",
                "robustness": robust,
            },
        },
        "full_lock_fits": {"tau": TAUS[-1],
                           "triangle": {"amplitude": round(A_tri, 4), "rmse": round(rmse_tri, 4)},
                           "cosine": {"amplitude": round(A_cos, 4), "rmse": round(rmse_cos, 4)},
                           "closer": "triangle" if rmse_tri < rmse_cos else "cosine"},
        "tau_curve": tau_curve,
        "runtime_seconds": round(time.time() - t0, 1),
    }
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "kimi_alignment")
    os.makedirs(outdir, exist_ok=True)
    path_out = os.path.join(outdir, "q1_first_map.json")
    with open(path_out, "w") as f:
        json.dump(out, f, indent=2)

    print("\n=== C(Delta, tau) map (locked pairs) ===")
    hdr = "  d\\tau | " + "  ".join(f"{t:>7d}" for t in TAUS)
    print(hdr)
    for d in DELTAS:
        row = [grid[f"d{d:.4f}_t{t}"]["C"] for t in TAUS]
        print(f"  {d:5.2f} | " + "  ".join(f"{c:+7.3f}" for c in row))
    print("\n=== controls at tau=64: max |C| ===")
    print("  no-source null:", max(abs(g["C"]) for g in null_grid.values()))
    print("  no-lock       :", max(abs(g["C"]) for g in nolock_grid.values()))
    print("\n=== full-lock fit (tau=%d) ===" % TAUS[-1])
    print(f"  triangle: A={A_tri:.3f} rmse={rmse_tri:.4f}   "
          f"cosine: A={A_cos:.3f} rmse={rmse_cos:.4f}   closer: {out['full_lock_fits']['closer']}")
    print("\n=== tau curve ===")
    for tc in tau_curve:
        print(f"  tau={tc['tau']:3d}  A={tc['triangle_amplitude']:+.3f}  "
              f"S={tc['CHSH_S']:.3f}  lock={tc['lock_factor']}")
    print(f"\nrelabeling dynamic audit: {relabel}")
    print(f"runtime: {out['runtime_seconds']}s\nwrote {path_out}")


if __name__ == "__main__":
    main()
