# spopt (development version)

* `spider_lines()` gains a `facilities` argument to draw links for selected facilities only. Filtering happens before any geometry is built, so one store's trade area from a large `huff()` result takes milliseconds instead of seconds.

## Bug fixes

* `bridge_islands = TRUE` could return a graph that was still disconnected:
  each component was linked once to its closest unit in any other component,
  so two islands closer to each other than to the mainland linked to each
  other only. Bridging now repeats until one component remains. `skater()` on
  such a forest previously returned one region more than requested.

## Performance (regionalization)

Measured on synthetic hexagonal grids with four smooth attributes
(`benchmarks/stage-timing.R`); partitions are unchanged on every captured
fixture (`benchmarks/equivalence.R`).

* `skater()`: cuts are scored from per-subtree sums in one pass per region and
  each region's best cut is cached, O(k n p) instead of O(k n^2 p).
  2,500 units: 2.2 s -> 0.09 s; 10,000 units: 31 s -> 0.8 s.
* `ward_spatial()`: candidate merges are kept in a priority queue with lazy
  invalidation instead of rescanning every adjacent pair per merge. Results are
  now deterministic (tie-breaking previously followed hash-map order).
  2,500 units: 0.39 s -> 0.09 s; 10,000 units: 7.3 s -> 0.8 s.
* `azp()` and `max_p_regions()` improvement phase: the contiguity check for a
  move stops as soon as the moving area's in-region neighbours are mutually
  reachable (breadth-first, stamp array instead of hash sets), and member
  removal is O(1) through a position index. 2,500 units, SA from a SKATER
  start: 6.6 s -> 0.9 s; 10,000 units: 107 s -> 5.1 s.
* Island bridging finds each component's closest cross-component pair with one
  distance matrix instead of a nested R loop (2,235 units, three islands:
  2.3 s -> 0.1 s excluding `poly2nb`).

## Regionalization: AZP and max-p search

* `azp()` gains `initial_regions` for warm starts from any partition, for
  example a `skater()` or `ward_spatial()` result. `n_regions` can be omitted
  when a warm start is given.
* The AZP search was rewritten on a shared local-search engine
  (`src/rust/src/region/local_search.rs`):
  - `method = "basic"` now runs Openshaw's procedure to convergence. One
    iteration is a full pass over all regions; previously one iteration was a
    single move, so the default budget of 100 iterations stopped far short of
    a local optimum.
  - `method = "sa"` cools once per sweep of the border and calibrates its start
    temperature from a sample of uphill moves. Previously the temperature was
    set to the objective divided by the number of regions, which accepted
    almost every uphill move for the whole run and returned a partition close
    to the random start.
  - `method = "tabu"` stops after `max(100, 10 * tabu_length)` consecutive
    non-improving moves instead of a fixed 100 moves in total.
  - Tabu and SA results are finished with a greedy pass, so every result is a
    local optimum under single-area moves.
  - Move evaluation uses the exact O(p) update of a within-region sum of
    squares, and the contiguity check runs only for moves that would be
    accepted. Searches that previously made 100 moves now run thousands of
    sweeps in comparable time.
  - `max_iterations` now defaults to 1000 and its meaning depends on `method`;
    see the documentation.
  - Metadata gains `initial_objective` and `warm_start`.
* `max_p_regions()` construction phase: when one frontier area is enough to
  reach the threshold, the one that overshoots it least is taken, so the
  remaining areas can form more regions; enclaves join the neighbouring region
  whose heterogeneity they raise least rather than a random one; and among
  constructions with the same number of regions the most homogeneous is kept.
* `max_p_regions()` improvement phase for the default objective
  (`homogeneous = TRUE`, `compact = FALSE`) is now a floor-constrained AZP:
  greedy pass, simulated annealing over the border for `n_sa_iterations`
  sweeps, greedy polish. `n_sa_iterations` now counts sweeps and defaults to
  1000. The previous annealing phase rarely improved on the construction. The
  compactness and heterogeneous objectives keep the previous routine.
* New tests for AZP and max-p search behaviour (`tests/testthat/test-azp.R`).

# spopt 0.1.3

* New `spider_lines()` turns location-allocation results from `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` into an sf layer of demand-to-facility lines. It supports split `cflp()` allocations, `huff()` probabilities, and road-following lines through a user-supplied `route_fun`. `r5r_route_fun()` builds a `route_fun` from an r5r network.
* `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` now add a `.cost` column to the demand output: the cost to the assigned facility (or primary store), taken from the original cost matrix before any NA/Inf replacement.
* `p_median()`, `cflp()`, and `mclp()` now record `weight_col` in their result metadata.
* `distance_matrix()` now reduces any non-POINT geometry (including MULTIPOINT) to its centroid, so each feature gets one row or column. Previously, MULTIPOINT inputs produced a mis-sized matrix in projected CRSs.

# spopt 0.1.2

* Initial CRAN release.
