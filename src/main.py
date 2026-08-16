#!/usr/bin/env python3

import argparse
from pathlib import Path
import polars as pl


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
        "--min-cluster-size",
        required=True,
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
    return (
        pl.scan_csv(input_path, separator="\t")
        .select(["productive", "cdr3_aa"])
        .filter(
            pl.col("productive").cast(pl.String).is_in(["T", "True"])
            & pl.col("cdr3_aa").is_not_null()
            & ~pl.col("cdr3_aa").str.contains("r\*")
        )
        .group_by("cdr3_aa")
        .agg(pl.len().alias("hcdr3_count"))
        .collect(engine="streaming")
    )


def main() -> None:
    args = get_args()
    input_path: Path = Path(args.input_file).resolve()

    # filtering report
    df = load_data(input_path)
    print(df.head(20))


if __name__ == "__main__":
    main()
