from pathlib import Path

import polars as pl
import pytest
from sklearn.metrics.pairwise import cosine_distances

from abkit.cluster import cluster_cdr3
from abkit.io import load_data

DATA_PATH = Path("data/PRJEB26509_IGH_100_PRJEB26509_IGH.tsv")


@pytest.fixture(scope="module")
def loaded():
    return load_data(DATA_PATH, "\t")


@pytest.fixture(scope="module")
def clustered(loaded):
    return cluster_cdr3(loaded, min_cluster_size=8)


def test_load_data(loaded):
    assert loaded.height > 0
    assert loaded["cdr3_aa"].n_unique() == loaded.height


def test_one_medoid_per_cluster(clustered):
    df_out, _ = clustered
    clusters = df_out.filter(pl.col("cluster_id") != -1)["cluster_id"].n_unique()
    medoids = df_out.filter(pl.col("is_medoid")).height
    assert clusters > 0
    assert clusters == medoids


def test_read_counts_preserved(loaded, clustered):
    df_out, _ = clustered
    assert df_out["read_count"].sum() == loaded["read_count"].sum()


def test_medoid_is_exact(clustered):
    df_out, X = clustered
    df_out = df_out.with_row_index("i")

    for cid in df_out.filter(pl.col("cluster_id") != -1)["cluster_id"].unique():
        members = df_out.filter(pl.col("cluster_id") == cid)
        costs = cosine_distances(X[members["i"].to_numpy()]).sum(axis=1)
        chosen = costs[members["is_medoid"].to_numpy()][0]
        assert chosen == pytest.approx(costs.min())
