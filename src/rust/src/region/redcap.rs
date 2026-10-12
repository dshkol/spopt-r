//! REDCAP (REgionalization with Dynamically Constrained Agglomerative
//! clustering and Partitioning)
//!
//! Guo, D. (2008). Regionalization with dynamically constrained agglomerative
//! clustering and partitioning. International Journal of Geographical
//! Information Science, 22(7), 801-823.
//!
//! REDCAP runs a contiguity-constrained hierarchical clustering, records for
//! each merge the cheapest contiguity edge between the two clusters joined,
//! and treats those n-1 edges as a spanning tree. The tree is then partitioned
//! exactly as in SKATER, by repeatedly removing the edge whose removal lowers
//! the total within-region sum of squared deviations the most. Compared with
//! SKATER's minimum spanning tree, the clustering-derived tree reflects
//! cluster-level rather than edge-level similarity, which in practice gives
//! more homogeneous regions at the same k.
//!
//! Implemented variant: full-order Ward linkage ("fullorder-ward"), the one
//! that performed best in the comparisons this package was benchmarked
//! against. Other linkages can reuse `cut_tree` with a different merge rule.

use extendr_api::prelude::*;

use super::skater::{compute_edge_dissimilarity, cut_tree};
use super::ward::{compute_objective, merge_sequence};

/// Solve REDCAP regionalization
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
    method: &str,
    floor_var: Option<Vec<f64>>,
    floor_value: f64,
) -> List {
    let n = attrs.nrows();
    let n_attrs = attrs.ncols();

    if n == 0 || n_regions == 0 || n_regions > n {
        return list!(
            labels = Vec::<i32>::new(),
            n_regions = 0i32,
            objective = 0.0f64
        );
    }

    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    let tree: Vec<(usize, usize)> = match method {
        "fullorder-ward" | "fullorder-wardlinkage" => {
            let costs = compute_edge_dissimilarity(&attr_vecs, adj_i, adj_j);
            merge_sequence(&attr_vecs, adj_i, adj_j, &costs)
                .into_iter()
                .map(|m| m.edge)
                .collect()
        }
        other => panic!(
            "Unknown REDCAP method '{}'; implemented: \"fullorder-ward\"",
            other
        ),
    };

    let labels0 = cut_tree(&attr_vecs, &tree, n_regions, floor_var.as_deref(), floor_value);
    let labels: Vec<i32> = labels0.iter().map(|&l| l + 1).collect();
    let actual_n_regions = *labels.iter().max().unwrap_or(&0) as usize;
    let objective = compute_objective(&attr_vecs, &labels, actual_n_regions);

    list!(
        labels = labels,
        n_regions = actual_n_regions as i32,
        objective = objective
    )
}
