import numpy as np

# MRH similarity: cosine similarity across flattened tensor vectors
def mrh_similarity(mrh1: Dict[str, Dict[str, float]], mrh2: Dict[str, Dict[str, float]]) -> float:
    def flatten_mrh(mrh: Dict[str, Dict[str, float]]) -> np.ndarray:
        flat = []
        for category in sorted(mrh.keys()):
            for key in sorted(mrh[category].keys()):
                flat.append(mrh[category][key])
        return np.array(flat)

    vec1 = flatten_mrh(mrh1)
    vec2 = flatten_mrh(mrh2)

    if np.linalg.norm(vec1) == 0 or np.linalg.norm(vec2) == 0:
        return 0.0  # Avoid division by zero

    similarity = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return round(float(similarity), 4)

# Create example MRH tensors
example_mrh_1 = generate_mrh_schema()
example_mrh_2 = generate_mrh_schema()

# Assign example weights to differentiate them
example_mrh_1["fractal_scale"]["biospheric"] = 0.8
example_mrh_1["informational_scope"]["technical"] = 0.6
example_mrh_1["geographic_scope"]["global"] = 0.7
example_mrh_1["action_scope"]["execution"] = 0.9
example_mrh_1["temporal_scope"]["days"] = 0.6

example_mrh_2["fractal_scale"]["biospheric"] = 0.7
example_mrh_2["informational_scope"]["technical"] = 0.5
example_mrh_2["geographic_scope"]["global"] = 0.6
example_mrh_2["action_scope"]["execution"] = 0.85
example_mrh_2["temporal_scope"]["days"] = 0.7

# Compute similarity
similarity_score = mrh_similarity(example_mrh_1, example_mrh_2)
similarity_score
