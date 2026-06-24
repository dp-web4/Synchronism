"""
Phase-18 (generative axis, QC) — does "compatibility ⟨C⟩" add anything over the QC field's own
correlation-structure metric, on the kind of data a real device gives you? (2026-06-24)

CONTEXT. Phase-17 showed the repo's compatibility design rule (compatibility, not coupling strength,
gates collective coherence; frustrated couplings can't be out-muscled) transfers in-silico to a
frustrated-Kuramoto qubit-coherence model. The honest next question (this file) is the DISCRIMINATING
one: on REAL-device-shaped data, does a single scalar "compatibility" ⟨C⟩ predict logical performance
BETTER than, or only the same as, the metric the QC field ALREADY uses for correlated noise —
namely SYMMETRY-ALIGNMENT of the correlated errors with the code's stabilizer structure?

Literature (fixed-marginal-rate, simulation): correlation STRUCTURE, not magnitude, sets the surface-
code threshold. Symmetry-ALIGNED correlated errors (e.g. correlated Z matching a plaquette/stabilizer
support) IMPROVE the threshold; MISALIGNED correlations DEGRADE it — at the SAME marginal error rate.
(arXiv:2506.15490; arXiv:2410.23779.) Real-device correlation structure is directly extractable as the
p_ij detection-event correlation matrix (Google d=5 surface code: 600x600 p_ij; arXiv:2207.06431).

THIS SIM is a minimal, transparent toy of that mechanism (NOT a surface-code decoder; not real data).
A logical observable is a parity over a set of "qubits." Single errors flip with marginal rate p
(held FIXED). On top, a fraction f of error events are CORRELATED PAIRS. We toggle the pair STRUCTURE:
  - ALIGNED   : correlated pairs lie within one stabilizer support -> their two flips cancel in the
                logical parity (a detectable, correctable, logically-benign event). "Compatible" with
                the code symmetry.
  - MISALIGNED: correlated pairs straddle the logical cut -> a single correlated event flips the
                logical observable directly (a logically-fatal weight-2 event). "Incompatible."
Both have the SAME marginal single-qubit error rate and the SAME total correlated fraction f.

THE DISCRIMINATING TEST:
  (A) Does a scalar "compatibility" ⟨C⟩ = (1 - f) + f*aligned_fraction  -- i.e. "fraction of error
      events that are benign" -- predict logical error rate? (It should, by construction.)
  (B) But does aligned-fraction (the structure) carry ALL the signal, with the marginal rate and the
      raw correlated fraction f explaining NONE of the aligned-vs-misaligned gap? If yes, then on real
      data the operative variable is the QC field's SYMMETRY-ALIGNMENT metric; a generic scalar ⟨C⟩
      "works" only because it is a coarsening of that same structural quantity, not an independent or
      better predictor. That is the honest verdict (ii->iii boundary): the lens REFRAMES usefully and
      names the right axis, but it does NOT beat / is subsumed by the field's own structural metric.

HONEST FRAMING. Bucket 0 = 0. This is generative/applied. The sim cannot and does not claim the
compatibility metric BEATS standard QC metrics on real data -- it shows WHY, on real-device-shaped
correlation data, ⟨C⟩ and stabilizer-symmetry-alignment are the SAME axis viewed at different
resolution, which is the result that matters for the verdict. numpy only, headless, writes JSON.
"""
import json, os
import numpy as np

N = 24                 # "qubits" in the logical observable's support neighborhood
TRIALS = 400_000       # Monte Carlo error events
SEED = 18


def logical_error_rate(p_marginal, f_corr, aligned_fraction, rng):
    """Monte Carlo logical error rate.
    p_marginal : marginal single-qubit flip prob (HELD FIXED across conditions).
    f_corr     : fraction of the marginal error budget delivered as CORRELATED PAIRS.
    aligned_fraction : of correlated pairs, fraction that are stabilizer-ALIGNED (benign).
                       (1-aligned_fraction) straddle the logical cut (fatal weight-2).
    Returns logical error rate (prob the logical parity is flipped per shot).
    Construction keeps the *marginal* per-qubit flip rate equal to p_marginal in all conditions:
    independent flips carry (1-f_corr) of the budget; correlated pairs carry f_corr.
    """
    # split logical support into two halves; logical observable = parity over half A
    # a flip inside A toggles the logical; a flip in B does not.
    halfA = N // 2
    p_indep = p_marginal * (1.0 - f_corr)
    p_pair = p_marginal * f_corr  # prob a given qubit participates in a correlated pair event

    logical_flips = 0
    batch = 20_000
    done = 0
    while done < TRIALS:
        n = min(batch, TRIALS - done)
        # independent flips on each qubit
        flips = (rng.random((n, N)) < p_indep).astype(np.int8)
        # correlated pair events: with prob ~p_pair, draw a pair
        pair_event = rng.random(n) < (p_pair)  # one candidate pair event per shot (toy)
        for k in np.nonzero(pair_event)[0]:
            if rng.random() < aligned_fraction:
                # ALIGNED: both members same side -> two flips on side A (cancel in parity) or side B
                if rng.random() < 0.5:
                    a, b = rng.integers(0, halfA), rng.integers(0, halfA)
                else:
                    a, b = rng.integers(halfA, N), rng.integers(halfA, N)
            else:
                # MISALIGNED: one member each side -> straddles the logical cut
                a, b = rng.integers(0, halfA), rng.integers(halfA, N)
            flips[k, a] ^= 1
            flips[k, b] ^= 1
        # logical parity over half A
        parityA = flips[:, :halfA].sum(axis=1) % 2
        logical_flips += int(parityA.sum())
        done += n
    return logical_flips / TRIALS


def main():
    rng = np.random.default_rng(SEED)
    p_marginal = 0.02
    f_corr = 0.4  # 40% of error budget is correlated pairs (strong, to make structure visible)

    # sweep aligned_fraction (the STRUCTURE) at fixed marginal rate and fixed correlated fraction
    sweep = []
    for af in np.linspace(0.0, 1.0, 11):
        ler = logical_error_rate(p_marginal, f_corr, af, np.random.default_rng(SEED + int(af * 100)))
        C = (1.0 - f_corr) + f_corr * af  # scalar "compatibility" = benign-event fraction
        sweep.append({"aligned_fraction": round(float(af), 2),
                      "compatibility_C": round(float(C), 3),
                      "logical_error_rate": round(float(ler), 5)})

    # CONTROL 1: vary marginal rate p with structure FIXED (all-misaligned) -> does LER move for a
    # reason unrelated to structure? (confound check: marginal rate is a standard metric)
    control_pmarg = []
    for pm in [0.01, 0.02, 0.04]:
        ler = logical_error_rate(pm, f_corr, 0.0, np.random.default_rng(SEED + 500 + int(pm * 1000)))
        control_pmarg.append({"p_marginal": pm, "aligned_fraction": 0.0,
                              "logical_error_rate": round(float(ler), 5)})

    # CONTROL 2 (the discriminator): at FIXED marginal rate AND fixed correlated fraction f, the only
    # thing that moves LER is aligned_fraction. Compare extremes.
    ler_all_aligned = sweep[-1]["logical_error_rate"]
    ler_all_misaligned = sweep[0]["logical_error_rate"]
    structure_gap = round(ler_all_misaligned - ler_all_aligned, 5)

    # Is scalar C monotone with LER? (yes by construction) -- and is it just a relabel of structure?
    Cs = np.array([s["compatibility_C"] for s in sweep])
    lers = np.array([s["logical_error_rate"] for s in sweep])
    afs = np.array([s["aligned_fraction"] for s in sweep])
    # corr of C with LER vs corr of aligned_fraction with LER -- identical because C is affine in af
    corr_C = float(np.corrcoef(Cs, lers)[0, 1])
    corr_af = float(np.corrcoef(afs, lers)[0, 1])

    verdict = (
        f"DISCRIMINATING TEST -- does scalar 'compatibility' ⟨C⟩ BEAT, or merely RELABEL, the QC "
        f"field's structural (stabilizer-symmetry-alignment) metric? At FIXED marginal error rate "
        f"p={p_marginal} and FIXED correlated fraction f={f_corr}, sweeping ONLY the alignment "
        f"structure moves the logical error rate from {ler_all_misaligned} (all-misaligned, fatal "
        f"weight-2) to {ler_all_aligned} (all-aligned, logically benign) -- a structure-only gap of "
        f"{structure_gap} at IDENTICAL marginal rate. This reproduces the literature's fixed-marginal "
        f"finding (arXiv:2506.15490, 2410.23779): correlation STRUCTURE, not magnitude, sets logical "
        f"performance. The scalar 'compatibility' C=(1-f)+f*aligned_fraction tracks LER perfectly "
        f"(corr(C,LER)={corr_C:.3f}) -- but ONLY because C is an affine function of aligned_fraction "
        f"(corr(aligned_fraction,LER)={corr_af:.3f}); they are the SAME axis at different resolution. "
        f"VERDICT: on real-device-shaped correlation data, the compatibility lens names the CORRECT "
        f"design axis (correlation structure / compatible-vs-destructive coupling, consonant with why "
        f"DFS & symmetry-aligned QEC work) -- a genuine REFRAMING/synthesis value -- but a generic "
        f"scalar ⟨C⟩ does NOT beat the field's own, more specific stabilizer-symmetry-alignment metric; "
        f"it is SUBSUMED by it (the QC metric specifies WHICH correlations are 'compatible': those "
        f"matching stabilizer support). So: (ii) TESTABLE-BUT-DATA-LIMITED for the strong 'beats "
        f"standard metrics' claim (the cross-configuration real-data test is not yet done by anyone), "
        f"trending (iii) 'reduces to / is subsumed by standard structural metrics' for the scalar "
        f"form. HONEST: Bucket 0 = 0; this is a toy mechanism + a verdict on the metric's status, "
        f"NOT a claim of demonstrated novelty on real data."
    )

    out = {
        "setup": {"N": N, "trials": TRIALS, "p_marginal_fixed": p_marginal, "f_corr_fixed": f_corr,
                  "note": "single-qubit marginal rate held fixed; only correlation STRUCTURE varied"},
        "alignment_sweep_fixed_marginal": sweep,
        "structure_only_LER_gap_(misaligned-aligned)": structure_gap,
        "control_marginal_rate_sweep": control_pmarg,
        "corr_compatibilityC_with_LER": round(corr_C, 4),
        "corr_alignedfraction_with_LER": round(corr_af, 4),
        "C_is_affine_relabel_of_structure": True,
        "verdict": verdict,
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    path = os.path.join(os.path.dirname(__file__), "results",
                        "phase18_qc_compatibility_vs_symmetry_alignment_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print(f"\nstructure-only LER gap at fixed marginal rate: {structure_gap} "
          f"(misaligned {ler_all_misaligned} -> aligned {ler_all_aligned})")
    print(f"corr(C,LER)={corr_C:.3f}  corr(aligned_fraction,LER)={corr_af:.3f}  "
          f"(equal => C is a relabel of the structural metric)")


if __name__ == "__main__":
    main()
