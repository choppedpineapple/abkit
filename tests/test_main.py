#!/usr/bin/env python3

from pathlib import Path
import polars as pl
from abkit.cluster import cluster_cdr3
from abkit.io import load_data

SAMPLE_PATH = Path("data/sample_airr.tsv")


def test_load_data():
    df = load_data(SAMPLE_PATH, "\t")
    assert df.height > 0
    assert df["cdr3_aa"].n_unique() == df.height


def test_cluster_cdr3():
    df = load_data(SAMPLE_PATH, "\t")
    df_out, _ = cluster_cdr3(df, min_cluster_size=5)
    n_clusters = df_out.filter(pl.col("cluster_id") != -1)["cluster_id"].n_unique()
    n_medoids = df_out.filter(pl.col("is_medoid")).height
    assert n_clusters == n_medoids
