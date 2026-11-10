# Drawing database

Drawings of 563 matroids made with `draw_matroid` from RealizationSpaces, with one row per run recording how each drawing was found, how readable it is, and why failed runs failed.

| Folder | Matroids | Drawn | Readable runs |
|---|---|---|---|
| [`polydb/`](polydb/README.md) | 483 simple rank-3 matroids on 6–9 elements from polyDB `Matroids.Small` | 459 (95%) | 2640 / 2754 drawn runs (96%) |
| [`integer_grid/`](integer_grid/README.md) | 80 matroids of 10–13 random points on the 4×4 integer grid | 79 (99%) | 275 / 472 drawn runs (58%) |

Each folder's README has the summary, the results per start, every reason for failure with the matroids it affects, and a description of its files: `jobs.tsv` holds every run, `drawings/` the first successful drawing of each (matroid, start), and `best/` the most readable drawing of each matroid.

## How the runs were made

- Software: RealizationSpaces with the drawing code (`draw_matroid`), Julia 1.12.6, Oscar 1.7.3, HomotopyContinuation 2.20.0, 2 threads per run.
- Every matroid is drawn with `start=:maximal_degenerations` and with `start=:circle`, each with seeds 1, 2 and 3. When the requested start fails, `draw_matroid` falls through to the next start in the order maximal degenerations, circle, exact; the column `start_used` records which start produced the drawing. The numerical irreducible decomposition start (`:nid`) was not used.
- Each run was limited to 1800 s and 5 GB of memory.
- Reruns: a run can fail for numerical reasons alone (the optimizer slides into a degeneration, or a run uses too much memory); such a run is rerun with a different seed. Ten runs were rerun this way, and their rows show the new seed: n9_111 (circle, seed 2), grid2_n13 and grid74_n13 (maximal_degenerations, seed 2) and grid76_n13 (seed 3) with seed 4, and grid67_n13 and grid75_n12 (maximal_degenerations, seeds 1–3) with seeds 4–6. All nine grid reruns drew; n9_111 failed as it should, since it is not realizable over ℚ.
- Matroid definitions: polyDB matroids are named `n<elements>_<index>` and defined by their nonbases in `polydb/matroids.tsv`. Grid matroids are named `grid<index>_n<points>`, and the `info` column of `integer_grid/jobs.tsv` lists their points.

## Readability and validation

- A drawing is accepted only if every nonbasis determinant is at most 10⁻⁷, every basis determinant is above 10⁻⁷, and the largest nonbasis determinant is at most 10⁻⁵ times the smallest basis determinant.
- Readability is the smaller of (closest two points) / (10% of the box width) and (closest point to a drawn line it is not on) / (6% of the box width), and 0 if a label leaves the box or overlaps a point or another label. A drawing is readable when its readability is at least 1.

## Reasons a matroid is not drawn

| Reason | polyDB | Grid |
|---|---|---|
| Not realizable over any field of characteristic 0 (proved by Oscar's exact realization space) | 16 | 0 |
| No real realization: the exact realization space is finite with no real point | 5 | 0 |
| No start found: the exact realization space has equations left in positive dimension, and the numerical starts failed | 3 (n9_108, n9_109, n9_239) | 1 (grid78_n13) |

The three polyDB matroids in the last row also have no real realization: in each, a restriction to 8 points has a finite realization space with no real point (n9_108 and n9_109: points 1–8; n9_239: points {1,2,3,5,6,7,8,9}). So every undrawn polyDB matroid provably has no drawing; grid78_n13 is undecided.

## Dependence on the Oscar version

An earlier run of the same scan with Oscar 1.8.0 drew the same polyDB matroids, but not grid23_n12 and grid53_n12. For the same matroid, `realization_space(M; char=0)` has the same dimension under both versions, but Oscar 1.8.0 eliminates fewer variables:

| Matroid | Oscar 1.7.3 | Oscar 1.8.0 |
|---|---|---|
| grid23_n12 | 3 variables, ideal 0, dim 3 | 4 variables, 1 equation, dim 3 |
| grid53_n12 | 4 variables, ideal 0, dim 4 | 5 variables, 1 equation, dim 4 |
| grid78_n13 | 4 variables, 1 equation, dim 3 | 5 variables, 3 equations, dim 3 |
| grid13_n13 | 4 variables, 1 equation, dim 3 | 6 variables, 14 equations, dim 3 |

With ideal 0 the exact start picks random rational parameters and draws the matroid; with an equation left in positive dimension it gives up. The other differences between the two runs were numerical failures that a different seed removed.
