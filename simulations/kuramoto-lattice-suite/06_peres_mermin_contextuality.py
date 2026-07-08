#!/usr/bin/env python3
"""
06 — Peres–Mermin square vs the CRT scanning model (explorer 2026-07-08).

The CRT temporal-scanning picture (site: /two-reframes) models superposition as a
fast deterministic cycle through definite states, sampled by the observer. At any
instant, every observable of the system has a definite value fixed by "where the
scan is" — independent of which *other* compatible observables are co-measured.
That is precisely a NON-CONTEXTUAL value assignment.

Kochen–Specker (1967) rules such assignments out for dim ≥ 3. The Peres–Mermin
square (Peres 1990, Mermin 1990) is the minimal state-independent proof in dim 4:
nine two-qubit observables in a 3×3 grid,

        A11 = X⊗I   A12 = I⊗X   A13 = X⊗X
        A21 = I⊗Z   A22 = Z⊗I   A23 = Z⊗Z
        A31 = X⊗Z   A32 = Z⊗X   A33 = Y⊗Y

Each row and each column is a mutually commuting (co-measurable) set. QM fixes the
product of each row to +I and of each column to +I, EXCEPT the third column, whose
product is −I. A scanning model must supply v(A) ∈ {+1, −1} for all nine at each
instant, satisfying all six product constraints — for EVERY instant of the scan,
because each constraint is verifiable on any state at any time (state-independent).

This script:
  1. DERIVES the six constraints from the operators (numpy) — nothing asserted.
  2. Exhaustively enumerates all 2^9 = 512 non-contextual assignments.
  3. Reports how many satisfy all six constraints (theorem: zero), and the maximum
     number of constraints any assignment satisfies (known: 5 of 6).
  4. States the consequence for the scanning model: no instant of any deterministic
     cycle can carry a consistent value set; time-averaging over the cycle cannot
     help, because each single co-measurement already verifies its constraint with
     certainty. The only escape is to make the scanned value depend on which
     compatible set is being co-measured — i.e., contextuality — which surrenders
     the "just sampling timing" content of the analogy.

Companion to 05 (CHSH): the same non-contextual real-valued ontology fails two
independent theorems — Bell (statistical, needs space-like separation) and
Kochen–Specker (algebraic, needs only one lab).

Output: results/06_peres_mermin_contextuality.json
"""

import itertools
import json
import os

import numpy as np

# ---------------------------------------------------------------- operators
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(a, b):
    return np.kron(a, b)


# The Peres–Mermin square, indexed [row][col]
SQUARE = [
    [kron(X, I2), kron(I2, X), kron(X, X)],
    [kron(I2, Z), kron(Z, I2), kron(Z, Z)],
    [kron(X, Z), kron(Z, X), kron(Y, Y)],
]
NAMES = [
    ["X⊗I", "I⊗X", "X⊗X"],
    ["I⊗Z", "Z⊗I", "Z⊗Z"],
    ["X⊗Z", "Z⊗X", "Y⊗Y"],
]

ID4 = np.eye(4, dtype=complex)


def check_qm_structure():
    """Derive (not assert) the co-measurability and product constraints."""
    facts = {"observables_squared_identity": True, "contexts": []}

    # every observable is ±1-valued: A² = I
    for r in range(3):
        for c in range(3):
            if not np.allclose(SQUARE[r][c] @ SQUARE[r][c], ID4):
                facts["observables_squared_identity"] = False

    contexts = []
    for r in range(3):  # rows
        ops = [SQUARE[r][c] for c in range(3)]
        contexts.append((f"row{r+1}", [NAMES[r][c] for c in range(3)], ops))
    for c in range(3):  # columns
        ops = [SQUARE[r][c] for r in range(3)]
        contexts.append((f"col{c+1}", [NAMES[r][c] for r in range(3)], ops))

    for label, names, ops in contexts:
        commuting = all(
            np.allclose(a @ b, b @ a) for a, b in itertools.combinations(ops, 2)
        )
        prod = ops[0] @ ops[1] @ ops[2]
        if np.allclose(prod, ID4):
            sign = +1
        elif np.allclose(prod, -ID4):
            sign = -1
        else:
            sign = None  # would falsify the construction
        facts["contexts"].append(
            {"context": label, "observables": names, "commuting": commuting,
             "product_sign": sign}
        )
    return facts


def exhaust_noncontextual(facts):
    """Try all 512 assignments v: observable -> ±1 against the derived constraints."""
    constraints = [
        (ctx["context"], ctx["product_sign"]) for ctx in facts["contexts"]
    ]
    # map context -> the 3 flat indices (r*3+c) it constrains
    idx = {
        "row1": [0, 1, 2], "row2": [3, 4, 5], "row3": [6, 7, 8],
        "col1": [0, 3, 6], "col2": [1, 4, 7], "col3": [2, 5, 8],
    }
    n_valid = 0
    best = 0
    best_examples = []
    for values in itertools.product([+1, -1], repeat=9):
        satisfied = sum(
            1 for label, sign in constraints
            if values[idx[label][0]] * values[idx[label][1]] * values[idx[label][2]] == sign
        )
        if satisfied == 6:
            n_valid += 1
        if satisfied > best:
            best = satisfied
            best_examples = [values]
        elif satisfied == best and len(best_examples) < 3:
            best_examples.append(values)
    return n_valid, best, best_examples


def parity_argument(facts):
    """The 3-line proof, checked numerically: product of all six constraint signs."""
    signs = [ctx["product_sign"] for ctx in facts["contexts"]]
    required = int(np.prod(signs))          # QM: (+1)^5 · (−1) = −1
    # each observable appears in exactly one row and one column, so any
    # assignment contributes v(A)² = +1 twice-over: the product of all six
    # constraint left-hand sides is +1 for EVERY assignment.
    achievable = +1
    return {"required_product_of_signs": required,
            "achievable_by_any_assignment": achievable,
            "contradiction": required != achievable}


def main():
    facts = check_qm_structure()
    n_valid, best, best_examples = exhaust_noncontextual(facts)
    parity = parity_argument(facts)

    result = {
        "experiment": "Peres–Mermin square vs non-contextual (CRT-scanning) value assignments",
        "qm_structure_derived": facts,
        "exhaustive_search": {
            "assignments_tried": 512,
            "assignments_satisfying_all_6_constraints": n_valid,
            "max_constraints_satisfiable": best,
        },
        "parity_argument": parity,
        "consequence_for_scanning_model": (
            "A deterministic scan through definite states is a non-contextual value "
            "assignment at every instant. Zero of 512 assignments satisfy the six "
            "QM-derived product constraints (best achievable: 5 of 6), so no instant "
            "of any cycle reproduces QM's dim-4 predictions; time-averaging cannot "
            "repair per-context certainties. Escape requires the scanned value to "
            "depend on the co-measured context — surrendering 'just sampling timing'."
        ),
    }

    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    out = os.path.join(os.path.dirname(__file__), "results",
                       "06_peres_mermin_contextuality.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print("Peres–Mermin square — QM structure (derived, not asserted):")
    for ctx in facts["contexts"]:
        print(f"  {ctx['context']}: {' · '.join(ctx['observables'])} "
              f"commuting={ctx['commuting']} product={'+I' if ctx['product_sign']==1 else '−I'}")
    print(f"\nExhaustive non-contextual search: {n_valid}/512 assignments satisfy all 6 "
          f"constraints; best achievable = {best}/6")
    print(f"Parity argument: constraints require product {parity['required_product_of_signs']}, "
          f"any assignment yields {parity['achievable_by_any_assignment']} → "
          f"contradiction = {parity['contradiction']}")
    print(f"\nWritten: {out}")


if __name__ == "__main__":
    main()
