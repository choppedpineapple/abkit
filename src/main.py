#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import polars as pl
from hdbscan import HDBSCAN
from scipy.sparse import csr_matrix
from sklearn.feature_extraction.text import TfidfVectorizer


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract, cluster, and plot antibody CDR3 repertoires."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=check_file_exist,
        help="Path to input igblastn AIRR report in TSV or CSV format.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=False,
        default="clustered_output.tsv",
        help="Path to output clustered file.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        required=False,
        default="clusters.pdf",
        help="Path to the output clusters plot file.",
    )
    parser.add_argument(
        "-m",
        "--min_cluster_size",
        required=False,
        default=5,
        type=int,
        help="Minimum cluster size for HDBSCAN.",
    )
    return parser.parse_args()


def check_file_exist(filepath: str) -> str:
    try:
        with open(filepath, "r") as _:
            return filepath
    except FileNotFoundError:
        raise argparse.ArgumentTypeError(f"File {filepath!r} not found!")


def load_data(input_path: Path) -> pl.DataFrame:
    schema = pl.scan_csv(input_path, separator="\t").collect_schema().names()
    cdr3_col = "cdr3_aa" if "cdr3_aa" in schema else "junction_aa"
    if cdr3_col not in schema:
        raise ValueError(
            f"Input file must contain 'cdr3_aa' or 'junction_aa'. Found: {schema}"
        )

    has_productive = "productive" in schema

    count_cols = ["duplicate_count", "consensus_count", "read_count", "count"]
    detected_count_col = next((c for c in count_cols if c in schema), None)

    cols_to_select = [cdr3_col]
    if has_productive:
        cols_to_select.append("productive")
    if detected_count_col:
        cols_to_select.append(detected_count_col)

    query = pl.scan_csv(input_path, separator="\t").select(cols_to_select)

    if has_productive:
        query = query.filter(
            pl.col("productive")
            .cast(pl.String)
            .str.to_uppercase()
            .is_in(["T", "TRUE", "1"])
        )

    agg_expr = (
        pl.col(detected_count_col).cast(pl.UInt32).sum().alias("read_count")
        if detected_count_col
        else pl.len().alias("read_count")
    )

    query = (
        query.filter(
            pl.col(cdr3_col).is_not_null()
            & (pl.col(cdr3_col).str.len_chars() >= 3)
            & ~pl.col(cdr3_col).str.contains(r"[\*\_\#\s]")
        )
        .group_by(pl.col(cdr3_col).alias("cdr3_aa"))
        .agg(agg_expr)
    )

    df = query.collect(engine="streaming")

    if df.height == 0:
        raise ValueError("No valid productive CDR3 sequences found after filtering.")

    return df


def cluster_cdr3(
    df: pl.DataFrame, min_cluster_size: int
) -> tuple[pl.DataFrame, csr_matrix]:
    vectorizer = TfidfVectorizer(analyzer="char_wb", ngram_range=(3, 3), norm="l2")
    X: csr_matrix = vectorizer.fit_transform(df["cdr3_aa"].to_numpy())

    if df.height < min_cluster_size:
        print(
            f"Warning: dataset size ({df.height}) is smaller than the min_cluster_size",
            "Assigning all to noise",
        )
        cluster_ids = np.full(df.height, -1, dtype=int)
    else:
        clusterer = HDBSCAN(
            min_cluster_size=min_cluster_size,
            metric="euclidean",
            cluster_selection_method="leaf",
        )
        cluster_ids = clusterer.fit_predict(X)

    df = df.with_columns(
        pl.Series("cluster_id", cluster_ids),
        pl.int_range(0, df.height).alias("row_idx"),
    )

    clustered_groups = (
        df.filter((pl.col("cluster_id") != -1) & pl.col("cluster_id").is_not_null())
        .group_by("cluster_id")
        .agg(pl.col("row_idx"))
    )

    medoid_row_indices = []
    for row in clustered_groups.iter_rows(named=True):
        indices = np.asarray(row["row_idx"])
        if len(indices) == 1:
            medoid_row_indices.append(indices[0])
            continue
        sub = X[indices]
        centroid_sum = np.asarray(sub.sum(axis=0)).ravel()
        medoid_local_idx = int(np.argmax(sub.dot(centroid_sum)))
        medoid_row_indices.append(indices[medoid_local_idx])

    is_medoid = np.zeros(df.height, dtype=bool)
    if medoid_row_indices:
        is_medoid[medoid_row_indices] = True

    df = df.with_columns(pl.Series("is_medoid", is_medoid)).drop("row_idx")
    return df, X


def main() -> None:
    args = get_args()
    input_path: Path = Path(args.input_file).resolve()

    # filtering report
    df = load_data(input_path)
    print(df.head(20))

    df_clustered, X = cluster_cdr3(df, args.min_cluster_size)

    df_clustered.write_csv(args.output_file, separator="\t")
    print("Clusters generated")


if __name__ == "__main__":
    main()
