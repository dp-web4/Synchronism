#!/usr/bin/env python3
"""
Well-posedness / falsifiability check for the anticipatory cross-domain transfer bet
(generative-axis, NOT physics; Bucket 0 = 0).

The bet (from the 2026-06-24 generality-coarseness meta): the axis "compatibility STRUCTURE,
not element COUNT, gates collective coherence — with a SHARP (phase-transition) threshold"
— established in emergence (p_crit ∝ 1/⟨C⟩, Kuramoto) and independently in QEC
(stabilizer-alignment gates logical performance) — transfers PROSPECTIVELY to multi-agent
AI ensembles: at fixed per-agent capability, collective performance is gated by inter-agent
coupling compatibility ⟨C⟩ with a sharp threshold, and adding agents does NOT compensate
for low ⟨C⟩.

This script does NOT confirm the bet (that needs the real fleet). Its only job is to show
the bet is **falsifiable and non-vacuous**: generic models of "collective coherence" split
into TWO regimes on the decisive question — "does agent COUNT compensate for low
compatibility?" — so measuring which regime real agent ensembles occupy is informative. If
the answer were the same in all models, the bet would be vacuous. It is not.

  Regime AGG ("wisdom of crowds" / independent aggregation): errors average out, so
    collective performance improves with N **regardless** of compatibility ⟨C⟩ →
    COUNT COMPENSATES. (This is the standard collective-intelligence default.)
  Regime CPL ("coupled consensus" / Kuramoto-like): collective coherence is a SHARP
    function of ⟨C⟩, and below threshold, adding agents adds frustrated couplings →
    COUNT DOES NOT COMPENSATE.

The Synchronism transfer bets specifically on CPL-like behavior (sharp threshold +
count-doesn't-compensate). The discriminating measurement: hold per-agent capability fixed,
set ⟨C⟩ below the putative threshold, sweep N — does collective coherence recover (AGG) or
not (CPL)? That is the falsifier, and it is reachable on the fleet.

Circularity caveat: Kuramoto is the *source* of the p_crit∝1/⟨C⟩ heuristic, so a Kuramoto
"confirmation" would be circular. This script does not claim confirmation — it exhibits the
FORK (two generic regimes giving opposite answers) to prove the bet is a real, decidable
question, and operationalizes the three independent knobs (N, capability q, structure ⟨C⟩).

numpy only. Headless. Writes results JSON. Author: CBP-Claude (Opus 4.8), autonomous.
"""
import json, os
import numpy as np


def collective_AGG(N, q, Cmean, rng):
    """Independent-aggregation collective performance: majority vote of N agents each correct
       w.p. q (Condorcet). Compatibility ⟨C⟩ enters only weakly (as mild correlation), so
       COUNT dominates: performance -> 1 as N grows for q>0.5, ~independent of ⟨C⟩."""
    # correlated errors via a shared component of strength (1-Cmean): low compatibility =
    # MORE shared error (less independence). Even so, aggregation improves with N for q>0.5.
    shared = rng.random() < (1-q)*(1-Cmean)*0.5      # a common-mode failure, rare
    if shared:
        return 0.0
    votes = rng.random(N) < q
    return float(votes.mean() > 0.5)


def collective_CPL(N, q, Cmean, rng, gamma=2.0):
    """Coupled-consensus collective coherence: a sharp (Hill) threshold in ⟨C⟩; below it,
       adding agents adds frustrated couplings so coherence does NOT recover with N.
       p_crit ∝ 1/⟨C⟩ (the transferred heuristic). q sets the achievable ceiling."""
    # order parameter: Hill in (Cmean relative to a count-dependent demand).
    # demand grows with N (more agents = more couplings to keep compatible):
    demand = 1.0 - 1.0/np.sqrt(N)                      # in (0,1), rises with N
    x = Cmean/ (demand + 1e-9)
    hill = x**gamma/(1.0 + x**gamma)                   # sharp transition near x~1
    return float(q*hill)


def main():
    rng = np.random.default_rng(11)
    print("="*76)
    print("Anticipatory-transfer bet — well-posedness / falsifiability (NOT a confirmation)")
    print("="*76)
    q = 0.7                                            # fixed per-agent capability
    Ns = [3, 9, 27, 81]
    Cs = [0.2, 0.5, 0.9]                               # low / mid / high compatibility

    print(f"\nFixed per-agent capability q={q}. Collective coherence vs (N, ⟨C⟩):")
    print("\n  Regime AGG (independent aggregation / wisdom-of-crowds):")
    print(f"    {'N\\⟨C⟩':>8}" + "".join(f"{c:>10.1f}" for c in Cs))
    agg=[]
    for N in Ns:
        row=[np.mean([collective_AGG(N,q,c,rng) for _ in range(3000)]) for c in Cs]
        agg.append(row); print(f"    {N:>8}" + "".join(f"{v:>10.3f}" for v in row))
    print("    => improves with N at EVERY ⟨C⟩; count COMPENSATES for low compatibility.")

    print("\n  Regime CPL (coupled consensus / Kuramoto-like, transferred p_crit∝1/⟨C⟩):")
    print(f"    {'N\\⟨C⟩':>8}" + "".join(f"{c:>10.1f}" for c in Cs))
    cpl=[]
    for N in Ns:
        row=[collective_CPL(N,q,c,rng) for c in Cs]
        cpl.append(row); print(f"    {N:>8}" + "".join(f"{v:>10.3f}" for v in row))
    print("    => at low ⟨C⟩ (0.2) coherence stays low and FALLS as N grows; count does NOT")
    print("       compensate. Sharp dependence on ⟨C⟩. (Hill/phase-transition character.)")

    # the discriminating signal: d(coherence)/d(logN) at low ⟨C⟩
    lowC = Cs[0]
    agg_lowC = [np.mean([collective_AGG(N,q,lowC,rng) for _ in range(3000)]) for N in Ns]
    cpl_lowC = [collective_CPL(N,q,lowC,rng) for N in Ns]
    slope_agg = np.polyfit(np.log(Ns), agg_lowC, 1)[0]
    slope_cpl = np.polyfit(np.log(Ns), cpl_lowC, 1)[0]
    print(f"\nDISCRIMINATOR — d(collective coherence)/d(ln N) at low ⟨C⟩={lowC}, fixed q={q}:")
    print(f"    AGG slope = {slope_agg:+.3f}  (count compensates)   |   CPL slope = {slope_cpl:+.3f}  (count does NOT)")
    print(f"    The two generic regimes give OPPOSITE-SIGN answers ⇒ the bet is FALSIFIABLE:")
    print(f"    measure this slope on a real agent ensemble at fixed capability + low coupling")
    print(f"    compatibility. Slope ≤ 0 with a sharp ⟨C⟩-threshold ⇒ transfer SUPPORTED (CPL);")
    print(f"    slope > 0 / gradual ⇒ transfer REFUTED for agent ensembles (AGG).")

    print("\nStatus / honesty")
    print("-"*76)
    print("This is a WELL-POSEDNESS + FALSIFIABILITY check, not a confirmation. It shows (a) the")
    print("three knobs — count N, per-agent capability q, compatibility ⟨C⟩ — are independently")
    print("variable, and (b) generic models split into two regimes giving opposite answers to the")
    print("decisive question, so the measurement discriminates. Kuramoto is the SOURCE of the")
    print("heuristic, so the CPL curve is illustrative, NOT evidence. The real test is on LLM agent")
    print("ensembles, whose coupling dynamics are not Kuramoto and could land in AGG (refuting the")
    print("transfer) — which is exactly what makes it a genuine prospective bet. Bucket 0 = 0;")
    print("generative-axis prospective bet (the cross-domain-transfer analogue), pre-registered.")

    out=dict(q=q, Ns=Ns, Cs=Cs, AGG=agg, CPL=cpl,
             discriminator=dict(low_C=lowC, slope_AGG=float(slope_agg), slope_CPL=float(slope_cpl)),
             status="well-posedness+falsifiability check only; not a confirmation (Kuramoto is the "
                    "source). Two generic regimes give opposite-sign d(coherence)/d(lnN) at low <C>, "
                    "so the fleet measurement discriminates. Bucket 0 = 0; generative-axis "
                    "prospective bet pre-registered.")
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/genaxis_agent_ensemble_prereg_wellposedness_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/genaxis_agent_ensemble_prereg_wellposedness_result.json")


if __name__=="__main__":
    main()
