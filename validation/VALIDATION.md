# Validation

`abkit` picks a representative sequence (the medoid) for each cluster.

The obvious way to do this is to measure the distance between every pair of
sequences and pick whichever one sits closest to all the others. That works,
but the cost grows with the square of the cluster size, so it gets slow fast.

`abkit` instead compares each sequence against the cluster's average. That
cost grows linearly. The question this document answers is whether the fast
way gives the same answer as the slow way.

Yes, exactly, provided two conditions hold:

1. Sequence vectors are length-normalised before use.
2. Distance is measured with cosine distance.

Both are true in `abkit`. This is not an approximation that happens to work
well. Under those two conditions the two methods are mathematically the same
calculation, so they cannot disagree.

The normalisation comes from `norm="l2"` in the TF-IDF vectoriser in
`src/abkit/cluster.py`. If that is ever changed, the guarantee below no longer holds
and nothing will warn you.

## How it was tested

Clusters were generated with a known shape: one large group of sequences plus
a smaller group sitting further away. The size of the smaller group was varied
from nothing up to half the cluster, since lopsided clusters are the hardest
case and the most realistic one — an expanded clone with a few scattered
relatives looks exactly like this.

For each setting, 20 clusters were generated and both methods run on each.

Two things were recorded. **Agreement** is how often the two methods picked the
same sequence. **Cost ratio** compares how good the two picks were: 1.000 means
identical, 1.500 means the fast method's pick was 50% further from the rest of
the cluster.

Reproduce with:

```bash
uv run python validation/medoid_check.py
```

## Results

### As `abkit` uses it - cosine distance, normalised

| Minority group | Agreement | Median ratio | Worst ratio |
|---|---|---|---|
| 0%  | 100% | 1.000 | 1.000 |
| 10% | 100% | 1.000 | 1.000 |
| 20% | 100% | 1.000 | 1.000 |
| 30% | 100% | 1.000 | 1.000 |
| 40% | 100% | 1.000 | 1.000 |
| 50% | 100% | 1.000 | 1.000 |

Identical answers everywhere. This also held when cluster size, number of
features and lobe separation were varied, which is expected - none of those
appear in the underlying maths.

### Without normalisation

| Minority group | Agreement | Median ratio | Worst ratio |
|---|---|---|---|
| 0%  | 95% | 1.000 | 1.001 |
| 10% | 80% | 1.000 | 1.017 |
| 20% | 70% | 1.000 | 1.035 |
| 30% |  0% | 1.672 | 1.734 |
| 40% |  0% | 1.293 | 1.318 |
| 50% | 30% | 1.007 | 1.025 |

Fails sharply. Most clusters are fine until the minority group reaches about
30%, at which point every cluster gets the wrong answer and the chosen sequence
can be nearly twice as far from the rest as it should be.

The reason: without normalisation, longer vectors pull the average towards
themselves more strongly. The average ends up reflecting vector size while the
distance measurement reflects sequence counts, so the two stop agreeing.

### With Euclidean distance instead of cosine

| Minority group | Agreement | Median ratio | Worst ratio |
|---|---|---|---|
| 0%  | 95% | 1.000 | 1.001 |
| 10% | 85% | 1.000 | 1.003 |
| 20% | 95% | 1.000 | 1.001 |
| 30% | 85% | 1.000 | 1.007 |
| 40% | 75% | 1.000 | 1.012 |
| 50% | 60% | 1.000 | 1.007 |

Close, but not exact. The two methods usually agree, and when they don't, the
fast method's pick is within about 1% of the best possible one. Agreement drops
as clusters become more lopsided.

With Euclidean distance, the sequence closest to the cluster's average is not
always the one closest to every other sequence, so the shortcut becomes an
approximation. It is only guaranteed to be exact with cosine distance.

## What this means in practice

The fast method is safe to use in `abkit` because both required conditions are
built in.

If you reuse this code elsewhere, keep them. Dropping normalisation is the more
dangerous mistake of the two: it looks fine on most clusters and then fails
badly on the lopsided ones, which are often the ones you care about.

## Limitations

- The mathematical equivalence was confirmed on synthetic clusters. The
  implementation was additionally run on a public AIRR dataset (PRJEB26509,
  12,552 unique CDR3 sequences), where every cluster returned exactly one medoid,
  including one cluster of 502 members.
- Where several sequences are equally good, the two methods may return
  different ones. The choice is arbitrary and the cost is identical.
- Only cosine and Euclidean distances were tested.
- The test dataset was already collapsed to unique sequences with almost no clonal
  expansion, so the read-count aggregation path is largely untested on real data.
  Behaviour on non-collapsed repertoires has not been checked.
