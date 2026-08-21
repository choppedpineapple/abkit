# abkit

> Fast antibody repertoire analysis and CDR3 sequence clustering.

`abkit` is a high-performance toolkit for AIRR-seq repertoire clustering and representative clone (medoid) discovery. Built with **Polars**, **sparse linear algebra**, and **HDBSCAN**.

---

## Features

- **Fast Ingestion:** High-throughput streaming, validation, and read-count aggregation of AIRR-compliant TSVs via Polars.
- **Sparse Feature Extraction:** Amino acid k-mer TF-IDF representation without dense matrix memory overhead.
- **Density-Based Clustering:** Unsupervised cluster identification using HDBSCAN to model varying clonal expansion densities and isolate noise.
- **Fast Medoid Discovery:** $O(N)$ centroid-projection medoid extraction via sparse matrix operations instead of naive $O(N^2)$ pairwise distance matrices.

---

## Quickstart

This project uses [uv](https://github.com/astral-sh/uv) for fast, reproducible environment and dependency management.

### 1. Installation

```bash
git clone [https://github.com/](https://github.com/)<username>/abkit.git
cd abkit
uv sync
