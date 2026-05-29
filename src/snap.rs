use crate::eigenvalue::EigenvalueQuantization;

/// Practical tools for working with quantized eigenvalue spectra
pub struct SnapTools;

impl SnapTools {
    /// Snap eigenvalues to nearest algebraic numbers (rationals and sqrt rationals)
    pub fn snap_eigenvalues(eigenvalues: &[f64], tolerance: f64) -> Vec<f64> {
        eigenvalues
            .iter()
            .map(|&lambda| {
                // Try rationals p/q for small q
                for q in 1..=20u64 {
                    let p = (lambda * q as f64).round() as i64;
                    let rational = p as f64 / q as f64;
                    if (lambda - rational).abs() <= tolerance {
                        return rational;
                    }
                }
                // Try sqrt(p/q)
                for q in 1..=10u64 {
                    for p in 1..=50u64 {
                        let val = (p as f64 / q as f64).sqrt();
                        if (lambda - val).abs() <= tolerance {
                            return val;
                        }
                        if (lambda + val).abs() <= tolerance {
                            return -val;
                        }
                    }
                }
                lambda
            })
            .collect()
    }

    /// Find the minimal polynomial approximately satisfied by each eigenvalue
    pub fn spectral_fingerprint(eigenvalues: &[f64]) -> Vec<String> {
        eigenvalues
            .iter()
            .map(|&lambda| {
                // Check if it's close to p/q
                for q in 1..=20u64 {
                    let p = (lambda * q as f64).round() as i64;
                    if (lambda - p as f64 / q as f64).abs() < 0.01 {
                        if q == 1 {
                            return format!("x - {} = 0", p);
                        }
                        return format!("{}x - {} = 0", q, p);
                    }
                }
                // Check quadratic: x² = p/q
                for q in 1..=10u64 {
                    for p in 1..=50u64 {
                        if (lambda * lambda - p as f64 / q as f64).abs() < 0.01 {
                            if q == 1 {
                                return format!("x² - {} = 0", p);
                            }
                            return format!("{}x² - {} = 0", q, p);
                        }
                    }
                }
                format!("≈ {:.4}", lambda)
            })
            .collect()
    }

    /// Check if a graph's spectrum is "snapped" (quantized)
    pub fn is_quantized(eigenvalues: &[f64], tolerance: f64) -> bool {
        let snap_dists = EigenvalueQuantization::rational_snap_distance(eigenvalues, 20);
        let avg_snap = snap_dists.iter().sum::<f64>() / snap_dists.len().max(1) as f64;
        avg_snap < tolerance
    }

    /// Conservation ratio after snapping eigenvalues
    pub fn snapped_cr(adj: &[Vec<f64>], tolerance: f64) -> f64 {
        let eigs = EigenvalueQuantization::eigenvalues(adj);
        let snapped = Self::snap_eigenvalues(&eigs, tolerance);

        if snapped.is_empty() {
            return 0.0;
        }

        let sum: f64 = snapped.iter().sum();
        let sum_sq: f64 = snapped.iter().map(|x| x * x).sum();

        if sum_sq.abs() < 1e-12 {
            return 0.0;
        }

        (sum * sum) / (snapped.len() as f64 * sum_sq)
    }

    /// Build a graph whose eigenvalues approximate a target spectrum
    /// Uses a simple gradient-free search over integer-weighted graphs
    pub fn approximate_spectrum(target: &[f64], n_nodes: usize) -> Vec<Vec<f64>> {
        let mut best_adj = vec![vec![0.0f64; n_nodes]; n_nodes];
        let mut best_error = f64::MAX;

        // Try several random integer graphs and pick the best
        let mut state: u64 = 42;
        for _ in 0..200 {
            // Generate candidate
            let mut adj = vec![vec![0.0f64; n_nodes]; n_nodes];
            for i in 0..n_nodes {
                for j in (i + 1)..n_nodes {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    let w = (state % 4) as f64;
                    adj[i][j] = w;
                    adj[j][i] = w;
                }
            }

            let eigs = EigenvalueQuantization::eigenvalues(&adj);
            let mut sorted_eigs: Vec<f64> = eigs.into_iter().collect();
            sorted_eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());

            // Compare with target (sorted)
            let mut sorted_target: Vec<f64> = target.to_vec();
            sorted_target.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let len = sorted_eigs.len().min(sorted_target.len());
            let error: f64 = (0..len)
                .map(|i| (sorted_eigs[i] - sorted_target[i]).powi(2))
                .sum();

            if error < best_error {
                best_error = error;
                best_adj = adj;
            }
        }

        best_adj
    }
}
