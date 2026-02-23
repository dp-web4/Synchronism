"""
Coupling-Coherence Experiment
Tests whether C(p) = tanh(γ · log(p/p_crit + 1)) describes coherence emergence
in multi-agent knowledge discovery with controllable coupling.

Designed in response to Andrei's AI challenge (2026-02-22).
"""

import numpy as np
import json
import os
import sys
from itertools import combinations
from datetime import datetime


class World:
    """Random directed knowledge graph serving as ground-truth oracle."""

    def __init__(self, n_nodes=20, n_edges=60, n_types=4, seed=None):
        self.n_nodes = n_nodes
        self.n_edges = n_edges
        self.n_types = n_types
        self.rng = np.random.default_rng(seed)
        self.ground_truth = self._generate()

    def _generate(self):
        """Generate random directed typed graph as binary tensor [N x N x T]."""
        gt = np.zeros((self.n_nodes, self.n_nodes, self.n_types), dtype=np.float64)
        # No self-loops
        possible = [(i, j, t)
                     for i in range(self.n_nodes)
                     for j in range(self.n_nodes)
                     for t in range(self.n_types)
                     if i != j]
        chosen = self.rng.choice(len(possible), size=min(self.n_edges, len(possible)), replace=False)
        for idx in chosen:
            i, j, t = possible[idx]
            gt[i, j, t] = 1.0
        return gt

    def entropy(self):
        """Information content of the world in bits."""
        total_slots = self.n_nodes * (self.n_nodes - 1) * self.n_types
        p_edge = self.n_edges / total_slots
        if p_edge <= 0 or p_edge >= 1:
            return 0.0
        return total_slots * (-p_edge * np.log2(p_edge) - (1 - p_edge) * np.log2(1 - p_edge))


class Agent:
    """Bayesian belief agent that observes the world and updates beliefs."""

    def __init__(self, n_nodes, n_types, noise_rate, rng):
        self.n_nodes = n_nodes
        self.n_types = n_types
        self.noise_rate = noise_rate
        self.rng = rng
        # Uninformative prior: log-odds = 0 (probability 0.5)
        self.log_odds = np.zeros((n_nodes, n_nodes, n_types), dtype=np.float64)

    @property
    def beliefs(self):
        """Convert log-odds to probabilities."""
        return 1.0 / (1.0 + np.exp(-self.log_odds))

    def observe(self, world, n_observations):
        """Observe random edges from the world with noise."""
        # Sample random (i, j, t) positions to observe
        indices = []
        for _ in range(n_observations):
            i = self.rng.integers(0, self.n_nodes)
            j = self.rng.integers(0, self.n_nodes)
            while j == i:
                j = self.rng.integers(0, self.n_nodes)
            t = self.rng.integers(0, self.n_types)
            indices.append((i, j, t))

        for i, j, t in indices:
            true_val = world.ground_truth[i, j, t]
            # Apply noise
            if self.rng.random() < self.noise_rate:
                observed = 1.0 - true_val
            else:
                observed = true_val

            # Bayesian update via log-odds
            # LR = P(obs|edge=1) / P(obs|edge=0)
            if observed == 1.0:
                lr = (1.0 - self.noise_rate) / self.noise_rate
            else:
                lr = self.noise_rate / (1.0 - self.noise_rate)
            self.log_odds[i, j, t] += np.log(lr)

    def receive_beliefs(self, other_beliefs, self_weight=0.7):
        """Update beliefs by averaging with another agent's beliefs."""
        my_beliefs = self.beliefs
        merged = self_weight * my_beliefs + (1.0 - self_weight) * other_beliefs
        # Convert back to log-odds, clipping to avoid inf
        merged = np.clip(merged, 1e-10, 1.0 - 1e-10)
        self.log_odds = np.log(merged / (1.0 - merged))

    def infer_graph(self, threshold=0.5):
        """Threshold beliefs into binary graph."""
        return (self.beliefs > threshold).astype(np.float64)

    def f1_score(self, world):
        """F1 score of inferred graph vs ground truth."""
        pred = self.infer_graph()
        gt = world.ground_truth
        tp = np.sum(pred * gt)
        fp = np.sum(pred * (1.0 - gt))
        fn = np.sum((1.0 - pred) * gt)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        if precision + recall == 0:
            return 0.0
        return 2 * precision * recall / (precision + recall)


def jaccard_similarity(beliefs_a, beliefs_b, threshold=0.5):
    """Jaccard similarity of predicted edge sets.

    Only counts edges where at least one agent has formed an opinion
    (moved away from the 0.5 prior). This avoids the saturation problem
    where vast numbers of unobserved edges (both at 0.5) inflate similarity.
    """
    pred_a = (beliefs_a > threshold).astype(np.float64)
    pred_b = (beliefs_b > threshold).astype(np.float64)
    intersection = np.sum(pred_a * pred_b)
    union = np.sum(np.maximum(pred_a, pred_b))
    if union == 0:
        return 1.0  # Both predict nothing → trivially agree
    return float(intersection / union)


class Experiment:
    """Single experiment run with fixed parameters."""

    def __init__(self, n_nodes=20, n_edges=60, n_types=4, n_agents=5,
                 noise_rate=0.15, obs_per_round=5, n_rounds=100,
                 coupling=0.0, seed=None):
        self.rng = np.random.default_rng(seed)
        self.world = World(n_nodes, n_edges, n_types, seed=self.rng.integers(0, 2**31))
        self.agents = [
            Agent(n_nodes, n_types, noise_rate, np.random.default_rng(self.rng.integers(0, 2**31)))
            for _ in range(n_agents)
        ]
        self.n_agents = n_agents
        self.obs_per_round = obs_per_round
        self.n_rounds = n_rounds
        self.coupling = coupling

    def measure_convergence(self):
        """Mean pairwise Jaccard similarity of agents' predicted edge sets.

        Uses Jaccard instead of JSD because JSD saturates when most beliefs
        are at the 0.5 prior (unobserved edges). Jaccard only measures
        agreement on edges the agents actually predict.
        """
        if self.n_agents < 2:
            return 1.0
        sims = []
        for a, b in combinations(range(self.n_agents), 2):
            sim = jaccard_similarity(self.agents[a].beliefs, self.agents[b].beliefs)
            sims.append(sim)
        return float(np.mean(sims))

    def measure_correctness(self):
        """Mean F1 of agents vs ground truth."""
        f1s = [a.f1_score(self.world) for a in self.agents]
        return float(np.mean(f1s))

    def run(self):
        """Execute the experiment, return per-round measurements."""
        history = []
        for r in range(self.n_rounds):
            # 1. Each agent observes the world
            for agent in self.agents:
                agent.observe(self.world, self.obs_per_round)

            # 2. Coupling: agents share beliefs
            if self.coupling > 0:
                # Collect current beliefs before sharing
                current_beliefs = [a.beliefs.copy() for a in self.agents]
                for i in range(self.n_agents):
                    for j in range(self.n_agents):
                        if i != j and self.rng.random() < self.coupling:
                            self.agents[i].receive_beliefs(current_beliefs[j])

            # 3. Measure
            c_conv = self.measure_convergence()
            c_corr = self.measure_correctness()
            c = float(np.sqrt(max(0, c_conv * c_corr)))

            history.append({
                'round': r + 1,
                'C_conv': round(c_conv, 4),
                'C_corr': round(c_corr, 4),
                'C': round(c, 4),
            })

        final = history[-1] if history else {}
        return {
            'final': final,
            'history': history,
            'world_entropy': round(self.world.entropy(), 2),
        }


def run_suite(coupling_levels, n_repetitions, progress=True, **world_params):
    """Run the full experimental suite."""
    results = []
    total = len(coupling_levels) * n_repetitions
    count = 0

    for p in coupling_levels:
        for rep in range(n_repetitions):
            seed = int(p * 10000) + rep * 100 + 42
            exp = Experiment(coupling=p, seed=seed, **world_params)
            result = exp.run()
            results.append({
                'coupling': round(p, 4),
                'repetition': rep,
                'final': result['final'],
                'world_entropy': result['world_entropy'],
                # Only store final 10 rounds to keep file size manageable
                'last_rounds': result['history'][-10:],
            })
            count += 1
            if progress and count % 50 == 0:
                print(f"  {count}/{total} runs complete ({100*count/total:.0f}%)")

    return results


def run_primary_experiment():
    """Run the primary experiment suite."""
    print("=" * 60)
    print("COUPLING-COHERENCE EXPERIMENT")
    print(f"Started: {datetime.now().isoformat()}")
    print("=" * 60)

    # Parameters chosen so individual agents CANNOT learn the world alone
    # but collectively with coupling they CAN:
    #   World: 12×11×3 = 396 possible edge positions
    #   Individual budget: 8 obs/round × 80 rounds = 640 obs (1.6 per position)
    #   Collective (K=5): 3200 obs (8.1 per position) — sufficient with coupling
    params = {
        'n_nodes': 12,
        'n_edges': 30,
        'n_types': 3,
        'n_agents': 5,
        'noise_rate': 0.15,
        'obs_per_round': 8,
        'n_rounds': 80,
    }

    # Fine resolution near 0 where the transition happens, coarser above 0.20
    coupling_levels = (
        [round(i * 0.005, 3) for i in range(21)] +  # 0.000 to 0.100, step 0.005
        [round(0.12 + i * 0.02, 2) for i in range(10)] +  # 0.12 to 0.30
        [round(0.35 + i * 0.05, 2) for i in range(14)]  # 0.35 to 1.00
    )
    n_repetitions = 20

    print(f"\nParameters: {params}")
    print(f"Coupling levels: {len(coupling_levels)} ({coupling_levels[0]} to {coupling_levels[-1]})")
    print(f"Repetitions per level: {n_repetitions}")
    print(f"Total runs: {len(coupling_levels) * n_repetitions}")
    print()

    results = run_suite(coupling_levels, n_repetitions, **params)

    output = {
        'experiment': 'coupling_coherence',
        'timestamp': datetime.now().isoformat(),
        'params': params,
        'coupling_levels': coupling_levels,
        'n_repetitions': n_repetitions,
        'results': results,
    }

    # Save
    out_dir = os.path.join(os.path.dirname(__file__), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'coupling_coherence_results.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")
    print(f"Finished: {datetime.now().isoformat()}")

    return output


def run_variations():
    """Run variation experiments for derived p_crit testing."""
    print("\n" + "=" * 60)
    print("VARIATION EXPERIMENTS")
    print("=" * 60)

    base = {
        'n_nodes': 12, 'n_edges': 30, 'n_types': 3,
        'n_agents': 5, 'noise_rate': 0.15, 'obs_per_round': 8, 'n_rounds': 80,
    }

    # Fewer coupling levels and repetitions for variations
    coupling_levels = [round(i * 0.05, 2) for i in range(21)]  # 0.00 to 1.00, step 0.05
    n_reps = 10

    variations = []

    # Vary noise
    for eta in [0.05, 0.15, 0.30]:
        params = {**base, 'noise_rate': eta}
        label = f"eta={eta}"
        print(f"\nRunning variation: {label}")
        results = run_suite(coupling_levels, n_reps, progress=False, **params)
        variations.append({'label': label, 'params': params, 'results': results})

    # Vary agents
    for k in [3, 5, 10]:
        params = {**base, 'n_agents': k}
        label = f"K={k}"
        print(f"Running variation: {label}")
        results = run_suite(coupling_levels, n_reps, progress=False, **params)
        variations.append({'label': label, 'params': params, 'results': results})

    # Vary world size (keep edge density ~ 30/396 ≈ 7.6% of possible edges)
    for n in [8, 12, 20]:
        possible = n * (n - 1) * base['n_types']
        e = max(10, int(0.076 * possible))
        params = {**base, 'n_nodes': n, 'n_edges': e}
        label = f"N={n}"
        print(f"Running variation: {label}")
        results = run_suite(coupling_levels, n_reps, progress=False, **params)
        variations.append({'label': label, 'params': params, 'results': results})

    # Vary observation rate
    for m in [2, 5, 10]:
        params = {**base, 'obs_per_round': m}
        label = f"m={m}"
        print(f"Running variation: {label}")
        results = run_suite(coupling_levels, n_reps, progress=False, **params)
        variations.append({'label': label, 'params': params, 'results': results})

    output = {
        'experiment': 'coupling_coherence_variations',
        'timestamp': datetime.now().isoformat(),
        'coupling_levels': coupling_levels,
        'n_repetitions': n_reps,
        'variations': [{'label': v['label'], 'params': v['params'],
                         'results': v['results']} for v in variations],
    }

    out_dir = os.path.join(os.path.dirname(__file__), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'coupling_coherence_variations.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nVariations saved to {out_path}")

    return output


if __name__ == '__main__':
    primary = run_primary_experiment()

    if '--variations' in sys.argv:
        run_variations()
    else:
        print("\nRun with --variations to also run parameter variation experiments.")
