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
    # peak at schema to handle columns
    schema = pl.scan_csv(input_path, separator="\t").collect_schema().names()

    cdr3_col = "cdr3_aa" if "cdr3_aa" in schema else "junction_aa"
    if cdr3_col not in schema:
        raise ValueError(
            f"Input file must contain 'cdr3_aa' or 'junction_aa'. Found: {schema}"
        )

    has_productive = "productive" in schema
    cols_to_select = [cdr3_col, "productive"] if has_productive else [cdr3_col]
    query = pl.scan_csv(input_path, separator="\t").select(cols_to_select)

    if has_productive:
        query = query.filter(
            pl.col("productive")
            .cast(pl.String)
            .str.to_uppercase()
            .is_in(["T", "TRUE", "1"])
        )

    # filtering out stops, ambiguous, missing entries, and short strings
    query = (
        query.filter(
            pl.col(cdr3_col).is_not_null()
            & (pl.col(cdr3_col).str.len_chars() >= 3)
            & ~pl.col(cdr3_col).str.contains(r"[\*\_\#\s]")
        )
        .group_by(pl.col(cdr3_col).alias("cdr3_aa"))
        .agg(pl.len().alias("read_count"))
    )

    df = query.collect(engine="streaming")

    if df.height == 0:
        raise ValueError(
                "No valid productive CDR3 sequences found after filtering"
        )

    return df


def main() -> None:
    args = get_args()
    input_path: Path = Path(args.input_file).resolve()

    # filtering report
    df = load_data(input_path)
    print(df.head(20))


if __name__ == "__main__":
    main()
