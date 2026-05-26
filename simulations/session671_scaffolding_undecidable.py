#!/usr/bin/env python3
"""
Session 671: Is "productive scaffolding" a real defense, or non-discriminating?

S670 gave a Tier-1 discriminator separating reparametrization (central object ===
a pre-existing quantity, no new prediction) from confirmed discovery (new confirmed
prediction). But there is a THIRD category the program leans on as its last defense
(S614/S615/S627: "wrong theories motivate right questions"; S663 option C
"productive scaffolding"): a GENERATIVE REFORMULATION that makes no new prediction
at proposal time but later cashes out into confirmed novel predictions.

Question: can ANY test, applied AT PROPOSAL TIME, distinguish a sterile
reparametrization from a generative-but-not-yet-cashed-out reformulation?

Method: take a historical reference class with KNOWN outcomes and score each on the
only properties observable at proposal time, plus the retrospective outcome.
  P1: did it motivate new questions / research?         (observable at proposal)
  P2: did it make a NEW confirmed prediction at proposal? (observable at proposal)
  P3: did it EVENTUALLY cash out into confirmed novelty?  (retrospective only)
If sterile and generative cases are identical on (P1,P2) and differ only on P3,
then no proposal-time test discriminates them -- including S670's Tier-1 criterion
(which keys on P2) and the "motivates questions" defense (which keys on P1).
"""

# name, class, P1 motivated_questions, P2 new_prediction_at_proposal,
#       P3 eventually_cashed_out, years_to_cashout (None if never), note
cases = [
    ("Ptolemaic epicycles",      "sterile",   True,  False, False, None,
     "fit data ever better; never produced a confirmed novel prediction; superseded"),
    ("Phlogiston",               "sterile",   True,  False, False, None,
     "organized 18thC chemistry questions; no confirmed novelty; replaced by O2"),
    ("Caloric theory of heat",   "sterile",   True,  False, False, None,
     "motivated calorimetry; no confirmed novelty; replaced by kinetic theory"),
    ("Heliocentrism (1543)",     "generative",True,  False, True,  64,
     "Copernicus: NOT more accurate than Ptolemy at proposal (still circular orbits "
     "+ epicycles); cashed out via Kepler(1609)/Newton(1687) -> new confirmed preds"),
    ("Lagrangian mechanics",     "generative",True,  False, True,  130,
     "1788: same predictions as Newton; cashed out via Noether(1918)/QFT -> "
     "conservation theorems, gauge theory"),
    ("Hamilton-Jacobi theory",   "generative",True,  False, True,  90,
     "1830s: same predictions as Newton/Lagrange; cashed out as the scaffold "
     "Schrodinger(1926) built wave mechanics on -> new confirmed predictions"),
    ("Synchronism (1995-2026)",  "UNKNOWN",   True,  False, False, None,
     "motivates questions (S614/627); 0 confirmed novel predictions in 670+ "
     "sessions / ~30 yr; cash-out status: none so far"),
]

def col(i): return [c[i] for c in cases]

print("=" * 78)
print("Reference class: properties observable AT PROPOSAL vs retrospective outcome")
print("=" * 78)
print(f"{'case':<26}{'class':<11}{'P1 ask?':<8}{'P2 pred?':<9}{'P3 cashed?':<11}{'yrs':<5}")
for c in cases:
    yrs = "" if c[5] is None else str(c[5])
    print(f"{c[0]:<26}{c[1]:<11}{str(c[2]):<8}{str(c[3]):<9}{str(c[4]):<11}{yrs:<5}")

# Can (P1,P2) -- the proposal-time observables -- separate sterile from generative?
sterile = [c for c in cases if c[1] == "sterile"]
generative = [c for c in cases if c[1] == "generative"]
def signature(c): return (c[2], c[3])     # (motivated_questions, new_prediction_at_proposal)
ster_sigs = set(signature(c) for c in sterile)
gen_sigs  = set(signature(c) for c in generative)
overlap = ster_sigs & gen_sigs

print()
print("=" * 78)
print("Discriminability at proposal time")
print("=" * 78)
print(f"  sterile    (P1,P2) signatures: {ster_sigs}")
print(f"  generative (P1,P2) signatures: {gen_sigs}")
print(f"  overlap: {overlap}")
if ster_sigs == gen_sigs:
    print("  => IDENTICAL. No proposal-time test on (motivated questions, new prediction)")
    print("     separates sterile from generative. Both 'motivate questions' (P1=True)")
    print("     and both make NO new prediction at proposal (P2=False).")
print()
print("  - S670 Tier-1 criterion keys on P2 (new prediction): fires 'reparametrization'")
print("    on ALL of them at proposal time, including heliocentrism and Bohr-Sommerfeld.")
print("  - The 'motivates right questions' defense keys on P1: TRUE for phlogiston and")
print("    caloric too. Non-discriminating -- it is satisfied by sterile theories.")
print("  - The ONLY separating variable is P3 (eventual cash-out), knowable only in")
print("    retrospect. Generativity is not a proposal-time property; it is a verdict")
print("    history returns later.")

print()
print("=" * 78)
print("Where Synchronism sits")
print("=" * 78)
syn = cases[-1]
print(f"  (P1,P2) signature = {signature(syn)} -- identical to BOTH reference classes.")
print("  So at proposal time Synchronism is in the reference class but its subset")
print("  (sterile vs generative) is UNDECIDABLE by any current test.")
print("  The only evidence that bears on P3 is the track record: 0 confirmed novel")
print("  predictions in 670+ sessions / ~30 years.")
print()
print("  Heliocentrism cashed out in ~64 yr, Hamilton-Jacobi ~90, Lagrangian ~130 --")
print("  all eventually produced a CONFIRMED NOVEL PREDICTION. Synchronism has produced")
print("  none in ~30 yr. The 'scaffolding' label claims generative status that only a")
print("  cash-out confers; it is unearned until a Tier-1 prediction appears.")
print()
print("  NOTE on the Bohr-Sommerfeld analogy the defense invokes: B-S is a POOR analogy")
print("  because it made new confirmed predictions AT PROPOSAL (hydrogen spectrum /")
print("  Rydberg 1913; fine structure 1916) -- P2=True, a genuine discovery, not a")
print("  P2=False scaffold. Synchronism (P2=False, 0 predictions) is LESS like")
print("  Bohr-Sommerfeld than the defense implies. The analogy is doubly unearned.")
print()
print("  HONEST VERDICT: not 'refuted' (undecidable), not 'vindicated as scaffolding'")
print("  (that label is unearned). Status = UNDECIDED, with a 0/670 null as the only")
print("  evidence -- which is compatible with sterile and provides no support for")
print("  generative. A single future Tier-1 confirmed prediction would flip it.")
