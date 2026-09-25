from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl


@dataclass(slots=True, frozen=True)
class AIRRDataLoader:
    input_path: Path
    delimiter: str = "\t"

    @classmethod
    def from_path(
        cls, input_path: str | Path, delim: str | None = None
    ) -> AIRRDataLoader:
        path = Path(input_path).resolve()
        if delim is None:
            delim = "\t" if path.suffix.lower() == ".tsv" else ","
        return cls(input_path=path, delimiter=delim)

    def load(self) -> pl.DataFrame:
        schema = (
            pl.scan_csv(self.input_path, separator=self.delimiter)
            .collect_schema()
            .names()
        )

        cdr3_col = "cdr3_aa" if "cdr3_aa" in schema else "junction_aa"
        if cdr3_col not in schema:
            raise ValueError(
                f"Input file must contain 'cdr3_aa' or 'junction_aa'. Found: {schema}"
            )

        has_productive = "productive" in schema
        count_cols = [
            "duplicate_count",
            "consensus_count",
            "read_count",
            "count",
        ]
        detected_count_col = next((c for c in count_cols if c in schema), None)

        cols_to_select = [cdr3_col]
        if has_productive:
            cols_to_select.append("productive")
        if detected_count_col:
            cols_to_select.append(detected_count_col)

        query = pl.scan_csv(self.input_path, separator=self.delimiter).select(
            cols_to_select
        )

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

        df = query.collect(engine="streaming").sort(
            ["read_count", "cdr3_aa"], descending=[True, False]
        )

        if df.height == 0:
            raise ValueError(
                "No valid productive CDR3 sequences found after filtering."
            )

        return df


def load_data(input_path: str | Path, delim: str | None = None) -> pl.DataFrame:
    return AIRRDataLoader.from_path(input_path, delim).load()
