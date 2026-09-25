# abkit

> Fast antibody repertoire analysis and CDR3 sequence clustering.

`abkit` is a toolkit for AIRR-seq repertoire clustering and representative
clone (medoid) discovery. Built with **Polars**, **sparse linear algebra**,
and **HDBSCAN**.

---

## Features

- **Fast ingestion** — streaming, validation and read-count aggregation of
  AIRR-compliant TSVs or CSVs via Polars.
- **Sparse feature extraction** — amino acid k-mer TF-IDF representation
  without dense matrix memory overhead.
- **Density-based clustering** — HDBSCAN, which handles varying clonal
  expansion densities and separates out noise rather than forcing every
  sequence into a cluster.
- **Linear-time medoid selection** — finds each cluster's representative
  sequence in O(N) instead of the usual O(N²) pairwise comparison. This is
  exact, not an approximation: see [validation](validation/VALIDATION.md).

---

## Install

Uses [uv](https://github.com/astral-sh/uv) for environment and dependency
management.

```bash
git clone https://github.com/choppedpineapple/abkit.git
cd abkit
uv sync
```

## Usage

```bash
uv run abkit cluster -i data/sample_airr.tsv
```

| Flag | Default | Description |
|---|---|---|
| `-i`, `--input_file` | required | IgBLAST AIRR report (TSV or CSV) |
| `-o`, `--output_file` | `clustered_output.tsv` | Output path |
| `-m`, `--min_cluster_size` | `5` | Minimum cluster size for HDBSCAN |

Input needs a `cdr3_aa` or `junction_aa` column. Productive sequences are kept
where a `productive` column exists, and read counts are aggregated from
whichever of `duplicate_count`, `consensus_count`, `read_count` or `count` is
present.

Output is the deduplicated sequence table plus `cluster_id`, `read_count` and
`is_medoid`. Sequences not assigned to a cluster get `cluster_id = -1`.

## Validation

The linear-time medoid selection is provably exact under two conditions:
length-normalised vectors and cosine distance. Both hold in `abkit`.

Tested across cluster shapes and sizes, and deliberately broken by removing
each condition in turn to confirm they matter. Full results and limitations in
[validation/VALIDATION.md](validation/VALIDATION.md).

```bash
uv run python validation/medoid_check.py
```

## Tests

```bash
uv run pytest
```

## Licence

MIT.
