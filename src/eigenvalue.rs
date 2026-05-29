use crate::power_iteration_eigenvalues;

/// Eigenvalues of integer-weighted Laplacians snap to algebraic numbers.
pub struct EigenvalueQuantization;

impl EigenvalueQuantization {
    /// Generate a random integer-weighted graph (deterministic PRNG)
    pub fn random_integer_graph(n: usize, max_weight: u64) -> Vec<Vec<f64>> {
        let mut adj = vec![vec![0.0f64; n]; n];
        // Simple LCG PRNG for determinism
        let mut state: u64 = 12345;
        for i in 0..n {
            for j in (i + 1)..n {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let w = (state % (max_weight + 1)) as f64;
                adj[i][j] = w;
                adj[j][i] = w;
            }
        }
        adj
    }

    /// Compute eigenvalues using power iteration with deflation
    pub fn eigenvalues(adj: &[Vec<f64>]) -> Vec<f64> {
        let n = adj.len();
        if n == 0 {
            return vec![];
        }

        // Compute Laplacian
        let degrees: Vec<f64> = adj.iter().map(|row| row.iter().sum()).collect();
        let mut lap = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    lap[i][j] = degrees[i];
                } else {
                    lap[i][j] = -adj[i][j];
                }
            }
        }

        power_iteration_eigenvalues(&lap, n)
    }

    /// How close are eigenvalues to rational numbers p/q with q <= max_denominator?
    pub fn rational_snap_distance(eigenvalues: &[f64], max_denominator: u64) -> Vec<f64> {
        eigenvalues
            .iter()
            .map(|&lambda| {
                let mut best = f64::MAX;
                for q in 1..=max_denominator {
                    let p = (lambda * q as f64).round() as i64;
                    let rational = p as f64 / q as f64;
                    let dist = (lambda - rational).abs();
                    if dist < best {
                        best = dist;
                    }
                }
                best
            })
            .collect()
    }

    /// How close are eigenvalues to roots of low-degree polynomials?
    /// Tests polynomials x^d + c_{d-1}*x^{d-1} + ... + c_0 for small integer coefficients
    pub fn algebraic_snap_distance(eigenvalues: &[f64], max_degree: usize) -> Vec<f64> {
        eigenvalues
            .iter()
            .map(|&lambda| {
                let mut best = f64::MAX;
                for deg in 1..=max_degree {
                    // Test integer polynomials with coefficients in [-3, 3]
                    let range = 3i64;
                    match deg {
                        1 => {
                            // x + c = 0 => x = -c
                            for c in -range..=range {
                                let root = -c as f64;
                                let dist = (lambda - root).abs();
                                if dist < best {
                                    best = dist;
                                }
                            }
                        }
                        2 => {
                            // x² + bx + c = 0
                            for b in -range..=range {
                                for c in -range..=range {
                                    let disc = (b * b - 4 * c) as f64;
                                    if disc >= 0.0 {
                                        let sq = disc.sqrt();
                                        let r1 = (-b as f64 + sq) / 2.0;
                                        let r2 = (-b as f64 - sq) / 2.0;
                                        let d1 = (lambda - r1).abs();
                                        let d2 = (lambda - r2).abs();
                                        if d1 < best {
                                            best = d1;
                                        }
                                        if d2 < best {
                                            best = d2;
                                        }
                                    }
                                }
                            }
                        }
                        3 => {
                            // x³ + ax² + bx + c = 0 — try integer roots
                            for a in -range..=range {
                                for b in -range..=range {
                                    for c in -range..=range {
                                        for trial_root in -5i64..=5 {
                                            let x = trial_root as f64;
                                            let val = x * x * x
                                                + a as f64 * x * x
                                                + b as f64 * x
                                                + c as f64;
                                            if val.abs() < 0.5 {
                                                // Near a root, check actual distance
                                                // Newton refine
                                                let mut rx = x;
                                                for _ in 0..10 {
                                                    let f = rx * rx * rx
                                                        + a as f64 * rx * rx
                                                        + b as f64 * rx
                                                        + c as f64;
                                                    let fp = 3.0 * rx * rx
                                                        + 2.0 * a as f64 * rx
                                                        + b as f64;
                                                    if fp.abs() < 1e-12 {
                                                        break;
                                                    }
                                                    rx -= f / fp;
                                                }
                                                let dist = (lambda - rx).abs();
                                                if dist < best {
                                                    best = dist;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        _ => {
                            // For higher degrees, just check integer roots
                            for trial_root in -5i64..=5 {
                                let dist = (lambda - trial_root as f64).abs();
                                if dist < best {
                                    best = dist;
                                }
                            }
                        }
                    }
                }
                best
            })
            .collect()
    }

    /// Quantization quality: minimum gap / maximum gap ratio among sorted eigenvalues
    pub fn quantization_quality(eigenvalues: &[f64]) -> f64 {
        if eigenvalues.len() < 2 {
            return 0.0;
        }
        let mut sorted: Vec<f64> = eigenvalues.iter().cloned().filter(|x| x.abs() > 1e-10).collect();
        if sorted.len() < 2 {
            return 0.0;
        }
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let gaps: Vec<f64> = sorted.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
        let min_gap = gaps.iter().cloned().fold(f64::MAX, f64::min);
        let max_gap = gaps.iter().cloned().fold(0.0f64, f64::max);

        if max_gap < 1e-12 {
            return 0.0;
        }
        min_gap / max_gap
    }

    /// Compare quantization for integer vs random-real weights
    pub fn integer_vs_real_quantization(n: usize, trials: usize) -> (f64, f64) {
        let mut int_total = 0.0f64;
        let mut real_total = 0.0f64;

        for t in 0..trials {
            // Integer graph
            let int_adj = Self::random_integer_graph(n, 5);
            let int_eigs = Self::eigenvalues(&int_adj);
            let int_snap = Self::rational_snap_distance(&int_eigs, 20);
            let int_avg = int_snap.iter().sum::<f64>() / int_snap.len().max(1) as f64;
            int_total += int_avg;

            // Real-weighted graph (using same PRNG but producing floats)
            let mut adj = vec![vec![0.0f64; n]; n];
            let mut state: u64 = 12345u64.wrapping_add(t as u64 * 997);
            for i in 0..n {
                for j in (i + 1)..n {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    let frac = (state as f64) / (u64::MAX as f64) * 5.0;
                    adj[i][j] = frac;
                    adj[j][i] = frac;
                }
            }
            let real_eigs = Self::eigenvalues(&adj);
            let real_snap = Self::rational_snap_distance(&real_eigs, 20);
            let real_avg = real_snap.iter().sum::<f64>() / real_snap.len().max(1) as f64;
            real_total += real_avg;
        }

        // Return as quality (lower snap distance = better quantization, so invert)
        let int_avg = int_total / trials as f64;
        let real_avg = real_total / trials as f64;
        // Return snap distances — integer should be smaller (closer to rationals)
        (int_avg, real_avg)
    }
}
