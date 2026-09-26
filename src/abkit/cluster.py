from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl
from hdbscan import HDBSCAN
from scipy.sparse import csr_matrix
from sklearn.feature_extraction.text import TfidfVectorizer


@dataclass(slots=True)
class CDR3Clusterer:
    min_cluster_size: int = 5
    ngram_range: tuple[int, int] = (3, 3)

    def fit_predict(self, df: pl.DataFrame) -> tuple[pl.DataFrame, csr_matrix]:
        vectorizer = TfidfVectorizer(
            analyzer="char_wb",
            ngram_range=self.ngram_range,
            norm="l2",
        )
        X: csr_matrix = vectorizer.fit_transform(df["cdr3_aa"].to_numpy())

        if df.height < self.min_cluster_size:
            print(
                f"Warning: dataset size ({df.height}) is smaller than min_cluster_size. "
                "Assigning all to noise."
            )
            cluster_ids = np.full(df.height, -1, dtype=int)
        else:
            # rows are L2-normalised, hence
            # euclidean distance ranks neighbours the same way as cosine.
            clusterer = HDBSCAN(
                min_cluster_size=self.min_cluster_size,
                metric="euclidean",
                cluster_selection_method="leaf",
            )
            cluster_ids = clusterer.fit_predict(X)

        df_out = df.with_columns(
            pl.Series("cluster_id", cluster_ids),
            pl.int_range(0, df.height).alias("row_idx"),
        )

        clustered_groups = (
            df_out.filter(
                (pl.col("cluster_id") != -1) & pl.col("cluster_id").is_not_null()
            )
            .group_by("cluster_id")
            .agg(pl.col("row_idx"))
        )

        medoid_row_indices: list[int] = []
        for row in clustered_groups.iter_rows(named=True):
            indices = np.asarray(row["row_idx"])
            if len(indices) == 1:
                medoid_row_indices.append(indices[0])
                continue

            sub = X[indices]
            centroid_sum = np.asarray(sub.sum(axis=0)).ravel()
            medoid_local_idx = int(np.argmax(sub.dot(centroid_sum)))
            medoid_row_indices.append(indices[medoid_local_idx])

        is_medoid = np.zeros(df_out.height, dtype=bool)
        if medoid_row_indices:
            is_medoid[medoid_row_indices] = True

        df_out = df_out.with_columns(pl.Series("is_medoid", is_medoid)).drop("row_idx")

        return df_out, X


def cluster_cdr3(
    df: pl.DataFrame, min_cluster_size: int = 5
) -> tuple[pl.DataFrame, csr_matrix]:
    return CDR3Clusterer(min_cluster_size=min_cluster_size).fit_predict(df)
