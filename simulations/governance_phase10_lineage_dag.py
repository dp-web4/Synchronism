#!/usr/bin/env python3
"""Phase 10: fork/merge lineage and provenance-preserving trust.

This is no longer primarily a Markov experiment. It is a provenance-DAG sanity check.

A society begins with a set of unique witnessed evidence receipts. Repeatedly:
1. it forks into two descendants;
2. both descendants inherit the same ancestral evidence;
3. each branch produces some new unique evidence;
4. the branches merge.

Compare two trust-accounting strategies:
- naive scalar inheritance: each child copies a parent score and a merge adds scores;
- provenance accounting: descendants carry receipt identities and a merge takes the set
  union before deriving any score.

The first double-counts shared ancestry and grows exponentially under repeated diamonds.
The second grows only when genuinely new evidence is created.

This is a governance/provenance toy, not a trust formula recommendation.
"""

from dataclasses import dataclass


@dataclass
class Node:
    name: str
    receipts: frozenset[str]
    naive_score: int
    parents: tuple[str, ...]


def make_receipts(prefix: str, count: int):
    return frozenset(f"{prefix}:{i:04d}" for i in range(count))


def fork(parent: Node, generation: int, branch: str, new_count: int):
    fresh = make_receipts(f"g{generation}-{branch}", new_count)
    return Node(
        name=f"g{generation}-{branch}",
        receipts=parent.receipts | fresh,
        naive_score=parent.naive_score + new_count,
        parents=(parent.name,),
    )


def merge(left: Node, right: Node, generation: int):
    return Node(
        name=f"g{generation}-merge",
        receipts=left.receipts | right.receipts,
        # Deliberately wrong control: scalar scores are added without provenance dedupe.
        naive_score=left.naive_score + right.naive_score,
        parents=(left.name, right.name),
    )


def main():
    initial_evidence = 100
    new_per_branch = 10
    generations = 8

    current = Node(
        name="root",
        receipts=make_receipts("root", initial_evidence),
        naive_score=initial_evidence,
        parents=(),
    )

    print(
        "generation  naive_inherited_score  unique_receipts  inflation_factor  "
        "new_evidence_this_generation"
    )
    print(
        f"{0:>10}  {current.naive_score:>21}  {len(current.receipts):>15}  "
        f"{current.naive_score / len(current.receipts):>16.3f}  {0:>28}"
    )

    for generation in range(1, generations + 1):
        left = fork(current, generation, "L", new_per_branch)
        right = fork(current, generation, "R", new_per_branch)
        current = merge(left, right, generation)

        unique = len(current.receipts)
        inflation = current.naive_score / unique
        print(
            f"{generation:>10}  {current.naive_score:>21}  {unique:>15}  "
            f"{inflation:>16.3f}  {2 * new_per_branch:>28}"
        )

    print("\nFinal DAG lesson:")
    print(
        "A merge may inherit multiple valid lineages, but shared ancestral evidence must "
        "be deduplicated before any trust derivation."
    )


if __name__ == "__main__":
    main()
