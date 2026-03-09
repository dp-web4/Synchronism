"""
Compatibility-Synthon Experiment (Phase 2)
Tests synthon formation with heterogeneous agents and compatibility structure.

Phase 1 (coupling_coherence_experiment.py) showed:
  - Clear Hill-function phase transition in coherence vs coupling frequency p
  - p_crit derivation fails because compatibility = 1.0 was assumed (all agents identical)
  - Hypothesis: p_crit ∝ 1/<compatibility>

Phase 2 tests:
  1. Does p_crit scale inversely with mean compatibility?
  2. Do heterogeneous agents spontaneously specialize?
  3. Can two-agent synthons achieve cross-type inference neither can alone?
  4. Does compatibility structure (uniform / block-diagonal / hierarchical) change the Hill exponent?

Key design choice: agents observe different SUBSETS of edge types, creating genuine
informational complementarity. A type-specialist agent can learn edges of its type
accurately but remains ignorant of other types. Synthon formation is the moment
this specialization enables capabilities neither agent has alone.

2026-03-08
"""

import numpy as np
import json
import os
from itertools import combinations
from datetime import datetime


# ---------------------------------------------------------------------------
# World
# ---------------------------------------------------------------------------

class World:
    """Random directed typed knowledge graph — same as Phase 1."""

    def __init__(self, n_nodes=12, n_edges=30, n_types=4, seed=None):
        self.n_nodes = n_nodes
        self.n_edges = n_edges
        self.n_types = n_types
        self.rng = np.random.default_rng(seed)
        self.ground_truth = self._generate()

    def _generate(self):
        gt = np.zeros((self.n_nodes, self.n_nodes, self.n_types), dtype=np.float64)
        possible = [
            (i, j, t)
            for i in range(self.n_nodes)
            for j in range(self.n_nodes)
            for t in range(self.n_types)
            if i != j
        ]
        chosen = self.rng.choice(len(possible), size=min(self.n_edges, len(possible)), replace=False)
        for idx in chosen:
            i, j, t = possible[idx]
            gt[i, j, t] = 1.0
        return gt

    def edges_for_type(self, t):
        """Ground truth restricted to edge type t."""
        return self.ground_truth[:, :, t]


# ---------------------------------------------------------------------------
# Heterogeneous Agent
# ---------------------------------------------------------------------------

class HeterogeneousAgent:
    """
    Bayesian agent that can observe only a SUBSET of edge types (its specialty).
    Compatibility with other agents is encoded in the coupling step.
    """

    def __init__(self, n_nodes, n_types, noise_rate, obs_types, agent_id, rng):
        """
        obs_types: list of type indices this agent can directly observe.
                   Other types are invisible to this agent's sensors.
        """
        self.n_nodes = n_nodes
        self.n_types = n_types
        self.noise_rate = noise_rate
        self.obs_types = set(obs_types)
        self.agent_id = agent_id
        self.rng = rng
        self.log_odds = np.zeros((n_nodes, n_nodes, n_types), dtype=np.float64)

    @property
    def beliefs(self):
        return 1.0 / (1.0 + np.exp(-self.log_odds))

    def observe(self, world, n_observations):
        """Observe random edges, restricted to observable types."""
        obs_type_list = sorted(self.obs_types)
        if not obs_type_list:
            return
        for _ in range(n_observations):
            i = self.rng.integers(0, self.n_nodes)
            j = self.rng.integers(0, self.n_nodes)
            while j == i:
                j = self.rng.integers(0, self.n_nodes)
            t = self.rng.choice(obs_type_list)
            true_val = world.ground_truth[i, j, t]
            if self.rng.random() < self.noise_rate:
                observed = 1.0 - true_val
            else:
                observed = true_val
            if observed == 1.0:
                lr = (1.0 - self.noise_rate) / self.noise_rate
            else:
                lr = self.noise_rate / (1.0 - self.noise_rate)
            self.log_odds[i, j, t] += np.log(lr)

    def receive_beliefs(self, other_beliefs, compatibility, self_weight=0.7):
        """
        Integrate another agent's beliefs, weighted by compatibility.
        compatibility ∈ [0,1]: how useful the other agent's compression is to this one.
        When compatibility = 0, the beliefs are ignored.
        When compatibility = 1, full averaging (as Phase 1).
        """
        if compatibility <= 0:
            return
        effective_weight = self_weight + (1.0 - self_weight) * (1.0 - compatibility)
        my_beliefs = self.beliefs
        merged = effective_weight * my_beliefs + (1.0 - effective_weight) * other_beliefs
        merged = np.clip(merged, 1e-10, 1.0 - 1e-10)
        self.log_odds = np.log(merged / (1.0 - merged))

    def infer_graph(self, threshold=0.5):
        return (self.beliefs > threshold).astype(np.float64)

    def f1_score(self, world, type_mask=None):
        """F1 vs ground truth, optionally restricted to specific types."""
        pred = self.infer_graph()
        gt = world.ground_truth
        if type_mask is not None:
            # Only score specified types
            mask = np.zeros_like(gt)
            for t in type_mask:
                mask[:, :, t] = 1.0
            pred = pred * mask
            gt = gt * mask
        tp = np.sum(pred * gt)
        fp = np.sum(pred * (1.0 - gt))
        fn = np.sum((1.0 - pred) * gt)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        if precision + recall == 0:
            return 0.0
        return 2 * precision * recall / (precision + recall)

    def specialization_score(self):
        """
        How differentiated is this agent's belief distribution across types?
        High score = agent has strong opinions on its specialty types, less on others.
        Returns per-type mean belief confidence (distance from 0.5 prior).
        """
        b = self.beliefs
        per_type = {}
        for t in range(self.n_types):
            confidence = np.mean(np.abs(b[:, :, t] - 0.5))
            per_type[t] = float(confidence)
        return per_type


def jaccard_similarity(beliefs_a, beliefs_b, threshold=0.5):
    pred_a = (beliefs_a > threshold).astype(np.float64)
    pred_b = (beliefs_b > threshold).astype(np.float64)
    intersection = np.sum(pred_a * pred_b)
    union = np.sum(np.maximum(pred_a, pred_b))
    if union == 0:
        return 1.0
    return float(intersection / union)


# ---------------------------------------------------------------------------
# Compatibility Matrix Factories
# ---------------------------------------------------------------------------

def make_compatibility_uniform(n_agents, mean_compat):
    """All pairs have the same compatibility."""
    C = np.full((n_agents, n_agents), mean_compat)
    np.fill_diagonal(C, 1.0)
    return C


def make_compatibility_block(n_agents, within_compat, between_compat):
    """
    Block-diagonal: high within community, low between.
    Assumes n_agents is even; two equal-size communities.
    """
    half = n_agents // 2
    C = np.full((n_agents, n_agents), between_compat)
    C[:half, :half] = within_compat
    C[half:, half:] = within_compat
    np.fill_diagonal(C, 1.0)
    return C


def make_compatibility_random(n_agents, mean_compat, spread, rng):
    """Random compatibility matrix symmetric around mean."""
    raw = rng.uniform(
        max(0, mean_compat - spread),
        min(1, mean_compat + spread),
        (n_agents, n_agents)
    )
    # Symmetrize
    C = (raw + raw.T) / 2
    np.fill_diagonal(C, 1.0)
    return np.clip(C, 0, 1)


# ---------------------------------------------------------------------------
# Agent Observation Type Assignments
# ---------------------------------------------------------------------------

def assign_obs_types_specialist(n_agents, n_types):
    """
    Each agent specializes in one or two edge types.
    Creates genuine informational complementarity — to learn the full world,
    agents must share across type boundaries.
    """
    obs_types = []
    for i in range(n_agents):
        # Each agent observes 1 or 2 types (overlap allowed but primary is unique)
        primary = i % n_types
        secondary = (i + 1) % n_types
        obs_types.append([primary, secondary])
    return obs_types


def assign_obs_types_generalist(n_agents, n_types):
    """All agents observe all types — same as Phase 1."""
    return [list(range(n_types)) for _ in range(n_agents)]


# ---------------------------------------------------------------------------
# Synthon Detection
# ---------------------------------------------------------------------------

def cross_type_f1(agents, world, type_a, type_b):
    """
    Measure collective ability to infer cross-type patterns.
    Collective belief = mean of all agent beliefs on type_b,
    compared to how well any single type_a-specialist could do on type_b alone.

    Returns (collective_f1_on_b, best_individual_f1_on_b).
    A synthon is forming when collective >> best individual.
    """
    # Collective prediction on type_b: majority vote across agents
    all_beliefs_b = np.stack([a.beliefs[:, :, type_b] for a in agents], axis=0)
    collective_beliefs_b = np.mean(all_beliefs_b, axis=0)
    collective_pred = (collective_beliefs_b > 0.5).astype(np.float64)
    gt_b = world.ground_truth[:, :, type_b]

    tp = np.sum(collective_pred * gt_b)
    fp = np.sum(collective_pred * (1 - gt_b))
    fn = np.sum((1 - collective_pred) * gt_b)
    p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    collective_f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0

    # Best individual on type_b
    individual_f1s = [a.f1_score(world, type_mask=[type_b]) for a in agents]
    best_individual_f1 = max(individual_f1s)

    return float(collective_f1), float(best_individual_f1)


def specialization_index(agents):
    """
    How much do agents differentiate by type?
    High index = agents have developed distinct expertise profiles.
    Returns mean (max_type_confidence - min_type_confidence) across agents.
    """
    scores = []
    for agent in agents:
        per_type = agent.specialization_score()
        confidences = list(per_type.values())
        scores.append(max(confidences) - min(confidences))
    return float(np.mean(scores))


# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------

class CompatibilityExperiment:
    """
    Phase 2 experiment: heterogeneous agents + compatibility matrix.
    """

    def __init__(self, n_nodes=12, n_edges=30, n_types=4, n_agents=6,
                 noise_rate=0.15, obs_per_round=8, n_rounds=80,
                 coupling=0.0, compatibility_matrix=None,
                 obs_type_assignment='specialist', seed=None):
        self.rng = np.random.default_rng(seed)
        self.world = World(n_nodes, n_edges, n_types,
                           seed=int(self.rng.integers(0, 2**31)))
        self.n_agents = n_agents
        self.n_types = n_types
        self.obs_per_round = obs_per_round
        self.n_rounds = n_rounds
        self.coupling = coupling

        # Compatibility matrix (n_agents x n_agents)
        if compatibility_matrix is None:
            self.compat = make_compatibility_uniform(n_agents, 1.0)
        else:
            self.compat = compatibility_matrix

        # Observation type assignments
        if obs_type_assignment == 'specialist':
            obs_types_list = assign_obs_types_specialist(n_agents, n_types)
        else:
            obs_types_list = assign_obs_types_generalist(n_agents, n_types)

        self.agents = [
            HeterogeneousAgent(
                n_nodes, n_types, noise_rate,
                obs_types=obs_types_list[i],
                agent_id=i,
                rng=np.random.default_rng(int(self.rng.integers(0, 2**31)))
            )
            for i in range(n_agents)
        ]

    def measure_convergence(self):
        if self.n_agents < 2:
            return 1.0
        sims = [
            jaccard_similarity(self.agents[a].beliefs, self.agents[b].beliefs)
            for a, b in combinations(range(self.n_agents), 2)
        ]
        return float(np.mean(sims))

    def measure_correctness(self):
        return float(np.mean([a.f1_score(self.world) for a in self.agents]))

    def run(self):
        history = []
        for r in range(self.n_rounds):
            # 1. Observe (each agent observes only its specialist types)
            for agent in self.agents:
                agent.observe(self.world, self.obs_per_round)

            # 2. Coupling with compatibility weighting
            if self.coupling > 0:
                current_beliefs = [a.beliefs.copy() for a in self.agents]
                for i in range(self.n_agents):
                    for j in range(self.n_agents):
                        if i != j and self.rng.random() < self.coupling:
                            self.agents[i].receive_beliefs(
                                current_beliefs[j],
                                compatibility=self.compat[i, j]
                            )

            # 3. Measure
            c_conv = self.measure_convergence()
            c_corr = self.measure_correctness()
            c = float(np.sqrt(max(0, c_conv * c_corr)))
            spec_idx = specialization_index(self.agents)

            history.append({
                'round': r + 1,
                'C_conv': round(c_conv, 4),
                'C_corr': round(c_corr, 4),
                'C': round(c, 4),
                'spec_idx': round(spec_idx, 4),
            })

        # End-of-run synthon metrics
        synthon_metrics = self._measure_synthon()

        final = history[-1] if history else {}
        return {
            'final': final,
            'synthon': synthon_metrics,
            'last_rounds': history[-10:],
            'world_entropy': round(self.world.entropy() if hasattr(self.world, 'entropy') else 0, 2),
        }

    def _measure_synthon(self):
        """Measure synthon formation indicators at end of run."""
        # Cross-type inference: collective vs best individual on each type
        cross_type_results = {}
        for t in range(self.n_types):
            # Find which agents specialize in type t
            t_specialists = [
                i for i, a in enumerate(self.agents) if t in a.obs_types
            ]
            other_types = [tt for tt in range(self.n_types) if tt != t]
            if other_types:
                cross_t = other_types[0]  # test cross-inference to first other type
                coll_f1, best_ind_f1 = cross_type_f1(self.agents, self.world, t, cross_t)
                cross_type_results[f't{t}_to_t{cross_t}'] = {
                    'collective_f1': round(coll_f1, 4),
                    'best_individual_f1': round(best_ind_f1, 4),
                    'emergence_ratio': round(coll_f1 / (best_ind_f1 + 1e-6), 3),
                }

        # Specialization index
        spec_idx = specialization_index(self.agents)

        # Per-agent type expertise
        agent_profiles = {}
        for i, agent in enumerate(self.agents):
            profile = agent.specialization_score()
            agent_profiles[f'agent_{i}'] = {
                'obs_types': sorted(agent.obs_types),
                'type_confidences': {str(k): round(v, 4) for k, v in profile.items()},
            }

        return {
            'cross_type_inference': cross_type_results,
            'specialization_index': round(spec_idx, 4),
            'agent_profiles': agent_profiles,
        }


# ---------------------------------------------------------------------------
# Suite Runners
# ---------------------------------------------------------------------------

def run_compatibility_sweep(coupling_levels, compatibility_values, n_repetitions=15,
                             n_nodes=12, n_edges=30, n_types=4, n_agents=6,
                             noise_rate=0.15, obs_per_round=8, n_rounds=80,
                             obs_type_assignment='specialist',
                             compat_structure='uniform'):
    """
    Main sweep: vary both coupling (p) and mean compatibility (<C>).
    Tests the prediction p_crit ∝ 1/<C>.
    """
    results = []
    total = len(coupling_levels) * len(compatibility_values) * n_repetitions
    count = 0
    rng = np.random.default_rng(42)

    for mean_compat in compatibility_values:
        for p in coupling_levels:
            for rep in range(n_repetitions):
                seed = int(mean_compat * 1000) * 10000 + int(p * 10000) + rep * 100 + 7

                # Build compatibility matrix based on structure
                if compat_structure == 'uniform':
                    compat = make_compatibility_uniform(n_agents, mean_compat)
                elif compat_structure == 'block':
                    within = min(1.0, mean_compat * 1.5)
                    between = max(0.0, mean_compat * 0.5)
                    compat = make_compatibility_block(n_agents, within, between)
                elif compat_structure == 'random':
                    compat = make_compatibility_random(n_agents, mean_compat, 0.2,
                                                        np.random.default_rng(seed))
                else:
                    compat = make_compatibility_uniform(n_agents, mean_compat)

                exp = CompatibilityExperiment(
                    n_nodes=n_nodes, n_edges=n_edges, n_types=n_types,
                    n_agents=n_agents, noise_rate=noise_rate,
                    obs_per_round=obs_per_round, n_rounds=n_rounds,
                    coupling=p, compatibility_matrix=compat,
                    obs_type_assignment=obs_type_assignment,
                    seed=seed,
                )
                result = exp.run()
                results.append({
                    'coupling': round(p, 4),
                    'mean_compatibility': round(mean_compat, 3),
                    'compat_structure': compat_structure,
                    'obs_assignment': obs_type_assignment,
                    'repetition': rep,
                    'final': result['final'],
                    'synthon': result['synthon'],
                    'last_rounds': result['last_rounds'],
                })
                count += 1
                if count % 50 == 0:
                    print(f"  {count}/{total} runs complete", flush=True)

    return results


def run_replacement_resilience(coupling, mean_compat, n_reps=10,
                                n_nodes=12, n_edges=30, n_types=4, n_agents=6,
                                noise_rate=0.15, obs_per_round=8, n_rounds=80):
    """
    Test identity persistence: does the synthon survive agent replacement?
    Run to convergence, replace one agent with a fresh one, measure recovery.
    """
    results = []
    for rep in range(n_reps):
        seed = 99999 + rep * 1000
        compat = make_compatibility_uniform(n_agents, mean_compat)

        # Phase 1: run to convergence
        exp = CompatibilityExperiment(
            n_nodes=n_nodes, n_edges=n_edges, n_types=n_types,
            n_agents=n_agents, noise_rate=noise_rate,
            obs_per_round=obs_per_round, n_rounds=n_rounds,
            coupling=coupling, compatibility_matrix=compat,
            obs_type_assignment='specialist', seed=seed,
        )
        result_before = exp.run()
        C_before = result_before['final'].get('C', 0)

        # Phase 2: replace agent 0 with a fresh agent of the same type
        fresh_rng = np.random.default_rng(seed + 50000)
        obs_types_list = assign_obs_types_specialist(n_agents, n_types)
        exp.agents[0] = HeterogeneousAgent(
            n_nodes=n_nodes, n_types=n_types, noise_rate=noise_rate,
            obs_types=obs_types_list[0], agent_id=0, rng=fresh_rng,
        )
        exp.n_rounds = n_rounds // 2  # Recovery window: half the original run

        result_after = exp.run()
        C_after = result_after['final'].get('C', 0)

        results.append({
            'C_before_replacement': round(C_before, 4),
            'C_after_replacement': round(C_after, 4),
            'recovery_ratio': round(C_after / (C_before + 1e-6), 3),
        })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    output_dir = os.path.join(os.path.dirname(__file__), 'results')
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M')

    # -----------------------------------------------------------------------
    # Experiment A: p_crit ∝ 1/<C> test
    # Primary sweep: 6 coupling levels × 5 compatibility values × 15 reps
    # -----------------------------------------------------------------------
    print("=== Experiment A: p_crit vs compatibility (uniform structure) ===")
    coupling_levels = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.70, 1.0]
    compatibility_values = [0.2, 0.4, 0.6, 0.8, 1.0]

    results_a = run_compatibility_sweep(
        coupling_levels=coupling_levels,
        compatibility_values=compatibility_values,
        n_repetitions=15,
        compat_structure='uniform',
        obs_type_assignment='specialist',
    )
    fname_a = os.path.join(output_dir, f'compatibility_uniform_{timestamp}.json')
    with open(fname_a, 'w') as f:
        json.dump({'experiment': 'A_uniform_compatibility', 'results': results_a}, f, indent=2)
    print(f"  Saved: {fname_a} ({len(results_a)} runs)")

    # -----------------------------------------------------------------------
    # Experiment B: Compatibility structure (block vs random vs uniform)
    # Fixed mean_compat=0.5, vary structure
    # -----------------------------------------------------------------------
    print("\n=== Experiment B: Compatibility structure at mean=0.5 ===")
    results_b = {}
    for structure in ['uniform', 'block', 'random']:
        print(f"  Structure: {structure}")
        res = run_compatibility_sweep(
            coupling_levels=coupling_levels,
            compatibility_values=[0.5],
            n_repetitions=15,
            compat_structure=structure,
            obs_type_assignment='specialist',
        )
        results_b[structure] = res

    fname_b = os.path.join(output_dir, f'compatibility_structure_{timestamp}.json')
    with open(fname_b, 'w') as f:
        json.dump({'experiment': 'B_compatibility_structure', 'results': results_b}, f, indent=2)
    print(f"  Saved: {fname_b}")

    # -----------------------------------------------------------------------
    # Experiment C: Specialist vs generalist at high compatibility
    # Does specialization enable more efficient synthon formation?
    # -----------------------------------------------------------------------
    print("\n=== Experiment C: Specialist vs generalist agents (compat=1.0) ===")
    results_c = {}
    for assignment in ['specialist', 'generalist']:
        print(f"  Assignment: {assignment}")
        res = run_compatibility_sweep(
            coupling_levels=coupling_levels,
            compatibility_values=[1.0],
            n_repetitions=15,
            compat_structure='uniform',
            obs_type_assignment=assignment,
        )
        results_c[assignment] = res

    fname_c = os.path.join(output_dir, f'specialist_vs_generalist_{timestamp}.json')
    with open(fname_c, 'w') as f:
        json.dump({'experiment': 'C_specialist_vs_generalist', 'results': results_c}, f, indent=2)
    print(f"  Saved: {fname_c}")

    # -----------------------------------------------------------------------
    # Experiment D: Replacement resilience
    # High coupling + high compat: does the synthon persist through agent swap?
    # -----------------------------------------------------------------------
    print("\n=== Experiment D: Replacement resilience ===")
    results_d = run_replacement_resilience(
        coupling=0.3, mean_compat=0.8, n_reps=20,
    )
    fname_d = os.path.join(output_dir, f'replacement_resilience_{timestamp}.json')
    with open(fname_d, 'w') as f:
        json.dump({'experiment': 'D_replacement_resilience', 'results': results_d}, f, indent=2)
    print(f"  Saved: {fname_d} ({len(results_d)} resilience runs)")

    # Summary
    print("\n=== Summary ===")
    # Quick check of A results at high vs low compatibility
    high_compat = [r for r in results_a if r['mean_compatibility'] >= 0.9 and r['coupling'] >= 0.05]
    low_compat  = [r for r in results_a if r['mean_compatibility'] <= 0.25 and r['coupling'] >= 0.05]
    if high_compat:
        mean_C_high = np.mean([r['final'].get('C', 0) for r in high_compat])
        print(f"  High compat (≥0.9), p≥0.05: mean C = {mean_C_high:.3f}")
    if low_compat:
        mean_C_low = np.mean([r['final'].get('C', 0) for r in low_compat])
        print(f"  Low compat (≤0.25), p≥0.05: mean C = {mean_C_low:.3f}")

    # Synthon emergence check
    synthon_examples = [
        r for r in results_a
        if r['coupling'] >= 0.1 and r['mean_compatibility'] >= 0.6
    ]
    if synthon_examples:
        emergence_ratios = []
        for r in synthon_examples:
            for key, vals in r.get('synthon', {}).get('cross_type_inference', {}).items():
                emergence_ratios.append(vals.get('emergence_ratio', 1.0))
        if emergence_ratios:
            print(f"  Mean collective/individual emergence ratio: {np.mean(emergence_ratios):.3f}")
            print(f"  Max emergence ratio: {np.max(emergence_ratios):.3f}")

    print("\nAll experiments complete.")
