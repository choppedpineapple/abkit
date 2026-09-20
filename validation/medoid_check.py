#!/usr/bin/env python3

import numpy as np
from scipy import sparse
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from sklearn.preprocessing import normalize


def make_cluster(n, n_features, seed, minority_frac, separation):
    # Two lobes: a majority plus a smaller group sitting further out
    rng = np.random.default_rng(seed)
    core_a = rng.random(n_features) < 0.1
    core_b = rng.random(n_features) < 0.1

    n_minor = int(n * minority_frac)
    rows = [core_a * rng.random(n_features) for _ in range(n - n_minor)]
    rows += [core_b * rng.random(n_features) * separation for _ in range(n_minor)]
    return sparse.csr_matrix(np.array(rows))


def distances(A, B, metric):
    return cosine_distances(A, B) if metric == "cosine" else euclidean_distances(A, B)


def exact_medoid(X, metric):
    # O(N^2)
    D = distances(X, X, metric)
    return int(D.sum(axis=1).argmin()), D


def projected_medoid(X, metric):
    # O(N)
    centroid = sparse.csr_matrix(X.mean(axis=0))
    return int(distances(X, centroid, metric).ravel().argmin())


def compare(X, metric, unit_norm):
    if unit_norm:
        X = normalize(X)
    true_i, D = exact_medoid(X, metric)
    approx_i = projected_medoid(X, metric)
    return true_i == approx_i, D[approx_i].sum() / D[true_i].sum()


def sweep(metric, unit_norm, seeds=20):
    print(f"\nmetric={metric}  unit_norm={unit_norm}")
    for frac in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        results = [
            compare(make_cluster(200, 500, s, frac, separation=3.0), metric, unit_norm)
            for s in range(seeds)
        ]
        agree = sum(a for a, _ in results) / len(results)
        ratios = [r for _, r in results]
        print(
            f"  frac={frac:.1f}  agree={agree:>4.0%}  "
            f"median={np.median(ratios):.3f}  worst={max(ratios):.3f}"
        )


if __name__ == "__main__":
    sweep("cosine", unit_norm=True)
    sweep("cosine", unit_norm=False)
    sweep("euclidean", unit_norm=True)
