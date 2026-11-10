# Integer-grid matroids

80 matroids, 480 runs (2 starts × 3 seeds each). How the runs were made and what readability means: see `../README.md`.

## Summary

| | Count |
|---|---|
| Matroids drawn by at least one run | 79 / 80 (99%) |
| Runs drawn | 472 / 480 (98%) |
| Readable runs (readability ≥ 1) | 275 / 472 (58%) |
| Readability of drawn runs: median / mean | 1.09 / 1.12 |
| Seconds per run: median / 90th percentile / max | 2.1 / 14.9 / 42.6 |

## By requested start

| Requested start | Runs | Drawn | Readable | Median s | Drawn by: maximal_degenerations / circle / exact |
|---|---|---|---|---|---|
| maximal_degenerations | 240 | 235 | 137 | 7.8 | 88 / 115 / 32 |
| circle | 240 | 237 | 138 | 1.1 | 0 / 208 / 29 |

## Not drawn

| Reason | Runs | Matroids |
|---|---|---|
| fail: start=:exact found no start: Oscar's realization space has equations left in positive dimension. | 8 | grid13_n13, grid78_n13 |

Matroids no run drew (1): grid78_n13.

## Files

- `jobs.tsv`: one row per run: `start` is the requested start, `start_used` the start that found the drawing (the starts are tried in the order maximal_degenerations, circle, exact), `readability` (≥ 1 is readable), `png` the drawing saved by that run, `reason` why a run failed, `info` the grid points.
- `drawings/`: the first successful seed of each (matroid, requested start).
- `best/`: the most readable of those drawings for each matroid.
- `jobs_*.tsv`: the same rows split by the parallel worker that ran them.
- `log_*.txt`, `stdout_*.txt`: the raw terminal output of each worker.
