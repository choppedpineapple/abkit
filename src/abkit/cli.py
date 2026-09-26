from __future__ import annotations

import argparse
from pathlib import Path

from abkit.cluster import cluster_cdr3
from abkit.io import load_data


def check_file_exist(filepath: str | Path) -> str:
    path = Path(filepath)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"File {str(filepath)!r} not found!")
    if path.suffix.lower() not in (".tsv", ".csv"):
        raise argparse.ArgumentTypeError(
            f"Unsupported extension {path.suffix!r} (expected .tsv or .csv)"
        )
    return str(path)


def min_size(value: str) -> int:
    n = int(value)
    if n < 2:
        raise argparse.ArgumentTypeError("min_cluster_size must be at least 2")
    return n


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abkit",
        description="Fast antibody repertoire analysis from AIRR reports.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    cluster_parser = sub.add_parser(
        "cluster", help="Cluster CDR3 sequences and identify medoids."
    )
    cluster_parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=check_file_exist,
        help="Path to AIRR-compliant TSV/CSV report.",
    )
    cluster_parser.add_argument(
        "-o",
        "--output_file",
        default="clustered_output.tsv",
        help="Path for output TSV/CSV (default: clustered_output.tsv).",
    )
    cluster_parser.add_argument(
        "-m",
        "--min_cluster_size",
        default=5,
        type=min_size,
        help="HDBSCAN min_cluster_size (default: 5).",
    )

    return parser


def run_cluster(args: argparse.Namespace) -> None:
    df = load_data(args.input_file)
    df_clustered, _ = cluster_cdr3(df, min_cluster_size=args.min_cluster_size)

    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_delim = "\t" if output_path.suffix.lower() == ".tsv" else ","
    df_clustered.write_csv(output_path, separator=out_delim)
    print(f"Wrote {output_path}")


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()

    if args.command == "cluster":
        run_cluster(args)


if __name__ == "__main__":
    main()
