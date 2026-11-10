# polyDB rank-3 simple matroids

483 matroids, 2898 runs (2 starts × 3 seeds each). How the runs were made and what readability means: see `../README.md`.

## Summary

| | Count |
|---|---|
| Matroids drawn by at least one run | 459 / 483 (95%) |
| Runs drawn | 2754 / 2898 (95%) |
| Readable runs (readability ≥ 1) | 2640 / 2754 (96%) |
| Readability of drawn runs: median / mean | 1.92 / 1.95 |
| Seconds per run: median / 90th percentile / max | 0.4 / 2.2 / 26.9 |

## By requested start

| Requested start | Runs | Drawn | Readable | Median s | Drawn by: maximal_degenerations / circle / exact |
|---|---|---|---|---|---|
| maximal_degenerations | 1449 | 1377 | 1318 | 0.5 | 1329 / 27 / 21 |
| circle | 1449 | 1377 | 1322 | 0.4 | 0 / 1276 / 101 |

## Not drawn

| Reason | Runs | Matroids |
|---|---|---|
| fail: not char-0 realizable | 96 | n7_15, n8_33, n8_56, n9_104, n9_110, n9_111, n9_154, n9_163, n9_164, n9_294, n9_295, n9_330 … |
| fail: no real realization | 30 | n8_24, n9_112, n9_153, n9_282, n9_339 |
| fail: start=:exact found no start: Oscar's realization space has equations left in positive dimension. | 18 | n9_108, n9_109, n9_239 |

Matroids no run drew (24): n7_15, n8_24, n8_33, n8_56, n9_104, n9_108, n9_109, n9_110, n9_111, n9_112, n9_153, n9_154, n9_163, n9_164, n9_239, n9_282, n9_294, n9_295, n9_330, n9_339, n9_340, n9_342, n9_370, n9_43.

## Files

- `jobs.tsv`: one row per run: `start` is the requested start, `start_used` the start that found the drawing (the starts are tried in the order maximal_degenerations, circle, exact), `readability` (≥ 1 is readable), `png` the drawing saved by that run, `reason` why a run failed, `info` the grid points.
- `matroids.tsv`: the definition of every matroid: its name, number of elements, and all its nonbases (3-element dependent sets) on the elements 1..n; an empty list is the uniform matroid U(3,n). Each is `matroid_from_nonbases(nonbases, n)` in Oscar.
- `drawings/`: the first successful seed of each (matroid, requested start).
- `best/`: the most readable of those drawings for each matroid.
- `jobs_*.tsv`: the same rows split by the parallel worker that ran them.
- `log_*.txt`, `stdout_*.txt`: the raw terminal output of each worker.
