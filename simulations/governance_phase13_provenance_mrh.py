#!/usr/bin/env python3
"""Governance Phase 13: query-specific provenance proof slices.

This is a deliberately transparent dependency-graph toy. It does NOT discover a novel
minimal provenance representation. Instead it illustrates an established provenance /
program-slicing idea in the Web4 continuity setting:

    full provenance graph -> claim-specific sufficient dependency closure

The synthetic graph contains:
- one 80-transition continuity path from a relying party's admission anchor;
- authorization / law / witness objects for each transition;
- a few policy-relevance probe receipts;
- trust evidence receipts;
- scoped registry/fencing proof objects for an exclusivity claim;
- many unrelated provenance subgraphs.

For each claim we compute a simple backward dependency closure. This is sufficient by
construction, not guaranteed mathematically minimal. A full-audit control requires the
whole graph.

The exclusivity proof includes a sparse authenticated-registry path of ceil(log2(N))
proof objects to illustrate compact non-membership/fencing evidence. No cryptographic
security proof is made here.

Formalism/governance toy only; no physics result.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict


def build_graph(
    noise_events: int = 10_000,
    main_steps: int = 80,
    trust_receipts: int = 60,
    registry_members: int = 1_000_000,
):
    deps: dict[str, set[str]] = defaultdict(set)
    kinds: dict[str, str] = {}
    nodes: set[str] = set()

    def add(node: str, kind: str, *requires: str) -> str:
        nodes.add(node)
        kinds[node] = kind
        for requirement in requires:
            nodes.add(requirement)
            deps[node].add(requirement)
        return node

    anchor = add("token:anchor", "token")
    previous = anchor
    lineage_links: list[str] = []
    policy_checks: list[str] = []

    # One historically valid continuity path from an admission anchor to current state.
    relevant_policy_steps = {5, 15, 25, 35, 45, 55, 65, 75}
    for i in range(1, main_steps + 1):
        authorization = add(f"auth:{i}", "authorization")
        witness_a = add(f"witness:{i}:a", "witness")
        witness_b = add(f"witness:{i}:b", "witness")
        law = add(f"law:{i // 10}", "law")
        event = add(
            f"event:{i}",
            "continuity_event",
            authorization,
            witness_a,
            witness_b,
            law,
        )
        token = add(f"token:{i}", "token")
        link = add(f"lineage:{i}", "lineage_link", previous, token, event)
        lineage_links.append(link)
        previous = token

        # Only a few amendments touch dimensions relevant to this relying party.
        if i in relevant_policy_steps:
            probe = add(f"probe:{i}", "policy_probe")
            check = add(f"policy-check:{i}", "policy_relevance", event, probe)
            policy_checks.append(check)

    # Positive descent from the admission anchor to the current token.
    descent = add("claim:descent", "claim", *lineage_links)

    # Scoped exclusivity evidence. The authenticated path stands in for a Merkle-like
    # membership/non-membership proof over an authoritative active-lineage registry.
    registry_root = add("registry:root", "registry_root")
    fence = add("fence:parent", "fence_receipt", registry_root)
    auth_path = [
        add(f"registry:path:{i}", "registry_path")
        for i in range(math.ceil(math.log2(max(2, registry_members))))
    ]
    exclusive_evidence = add(
        "evidence:exclusive", "exclusive_evidence", fence, *auth_path
    )
    exclusive = add(
        "claim:exclusive-continuity", "claim", descent, exclusive_evidence
    )

    # Policy compatibility: lineage plus only policy changes/probes relevant to R.
    compatibility = add(
        "claim:policy-compatibility", "claim", descent, *policy_checks
    )

    # Trust: lineage plus unique evidence actually consumed by one derivation model.
    evidence = [
        add(f"trust-evidence:{i}", "trust_evidence")
        for i in range(trust_receipts)
    ]
    trust_model = add("trust:model", "trust_model")
    trust_projection = add(
        "evidence:trust-projection", "trust_projection", trust_model, *evidence
    )
    trust = add("claim:trust-current", "claim", descent, trust_projection)

    # Large amount of unrelated historical provenance. It is intentionally disconnected
    # from the four claims above; a full audit must retain it, targeted claims need not.
    for i in range(noise_events):
        source = add(f"noise:{i}:source", "noise")
        event = add(f"noise:{i}:event", "noise", source)
        add(f"noise:{i}:result", "noise", event)
        if i % 5 == 0:
            add(f"noise:{i}:extra", "noise", event)

    # Full-audit control. Build this last so it explicitly depends on every pre-existing
    # provenance object.
    audit_inputs = sorted(nodes)
    full_audit = add("claim:full-audit", "claim", *audit_inputs)

    claims = {
        "descent": descent,
        "exclusive": exclusive,
        "compatibility": compatibility,
        "trust": trust,
        "full_audit": full_audit,
    }
    return deps, kinds, nodes, claims


def dependency_closure(deps: dict[str, set[str]], root: str) -> set[str]:
    """Simple backward slice: all objects recursively required by root."""
    seen: set[str] = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(deps.get(node, ()))
    return seen


def print_claim_table(deps, kinds, nodes, claims):
    total = len(nodes)
    print(f"full synthetic provenance objects: {total}")
    print("\nclaim            proof objects   fraction of full graph")
    print("---------------- -------------   ----------------------")
    for name in ("descent", "exclusive", "compatibility", "trust", "full_audit"):
        proof = dependency_closure(deps, claims[name])
        print(f"{name:16s} {len(proof):13d}   {len(proof)/total:21.6%}")

    print("\nProof-object kinds by claim")
    for name in ("descent", "exclusive", "compatibility", "trust"):
        proof = dependency_closure(deps, claims[name])
        counts = Counter(kinds.get(node, "unknown") for node in proof)
        common = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
        print(f"{name:16s} {common}")


def noise_sweep(registry_members: int):
    print("\nUnrelated-history sweep")
    print("noise_events   full_graph   descent   exclusive   compatibility   trust")
    for noise in (0, 100, 1_000, 10_000, 50_000):
        deps, kinds, nodes, claims = build_graph(
            noise_events=noise, registry_members=registry_members
        )
        sizes = {
            name: len(dependency_closure(deps, root))
            for name, root in claims.items()
        }
        print(
            f"{noise:12d}   {len(nodes):10d}   {sizes['descent']:7d}   "
            f"{sizes['exclusive']:9d}   {sizes['compatibility']:13d}   "
            f"{sizes['trust']:5d}"
        )


def registry_sweep(noise_events: int):
    print("\nAuthenticated-registry-size sweep")
    print("registry entries   path objects   exclusive proof objects")
    for registry_members in (10, 100, 1_000, 10_000, 1_000_000, 1_000_000_000):
        deps, kinds, nodes, claims = build_graph(
            noise_events=noise_events, registry_members=registry_members
        )
        proof = dependency_closure(deps, claims["exclusive"])
        path_len = math.ceil(math.log2(max(2, registry_members)))
        print(f"{registry_members:16d}   {path_len:12d}   {len(proof):23d}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--noise-events", type=int, default=10_000)
    parser.add_argument("--registry-members", type=int, default=1_000_000)
    args = parser.parse_args()

    deps, kinds, nodes, claims = build_graph(
        noise_events=args.noise_events,
        registry_members=args.registry_members,
    )
    print_claim_table(deps, kinds, nodes, claims)
    noise_sweep(args.registry_members)
    registry_sweep(args.noise_events)


if __name__ == "__main__":
    main()
