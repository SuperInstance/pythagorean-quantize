use crate::compute_cr;

/// Eisenstein integers: a + bω where ω = e^(2πi/6) = (-1+i√3)/2
/// The hexagonal lattice.
pub struct EisensteinInt;

impl EisensteinInt {
    /// The fundamental sixth root of unity: ω = (-1/2, √3/2)
    pub fn omega() -> (f64, f64) {
        (-0.5, (3f64).sqrt() / 2.0)
    }

    /// Norm: |a + bω|² = a² - ab + b²
    pub fn norm(a: i64, b: i64) -> u64 {
        let a = a as i64;
        let b = b as i64;
        (a * a - a * b + b * b) as u64
    }

    /// Generate Eisenstein integers with norm <= max_norm
    pub fn with_norm(max_norm: u64) -> Vec<(i64, i64, u64)> {
        let mut result = Vec::new();
        let bound = (max_norm as f64).sqrt() as i64 + 1;
        for a in -bound..=bound {
            for b in -bound..=bound {
                let n = Self::norm(a, b);
                if n <= max_norm && n > 0 {
                    result.push((a, b, n));
                }
            }
        }
        result.sort_by_key(|x| x.2);
        result.dedup_by_key(|x| (x.0, x.1));
        result
    }

    /// Build hexagonal lattice adjacency (6 neighbors in Eisenstein lattice)
    pub fn lattice_graph(max_norm: u64) -> Vec<Vec<f64>> {
        let points = Self::with_norm(max_norm);
        let n = points.len();
        let mut adj = vec![vec![0.0f64; n]; n];

        // The 6 unit neighbors in Eisenstein lattice: ±1, ±ω, ±(1+ω)
        let neighbors = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)];

        let mut idx_map = std::collections::HashMap::new();
        for (i, &(a, b, _)) in points.iter().enumerate() {
            idx_map.insert((a, b), i);
        }

        for (i, &(a, b, _)) in points.iter().enumerate() {
            for &(da, db) in &neighbors {
                if let Some(&j) = idx_map.get(&(a + da, b + db)) {
                    adj[i][j] = 1.0;
                }
            }
        }

        adj
    }

    /// Conservation ratio of the hexagonal lattice
    pub fn hexagonal_cr(n_rings: usize) -> f64 {
        // n_rings of hexagonal lattice: norm increases
        let max_norm = (n_rings as u64).pow(2);
        let adj = Self::lattice_graph(max_norm);
        compute_cr(&adj)
    }

    /// Compare hexagonal CR vs square CR
    pub fn hex_vs_square_cr(n: usize) -> (f64, f64) {
        let max_norm = (n as u64).pow(2);
        let hex_adj = Self::lattice_graph(max_norm);
        let hex_cr = compute_cr(&hex_adj);

        // Build square lattice (Z²) with similar number of points
        let bound = n as i64;
        let mut sq_points = Vec::new();
        for x in -bound..=bound {
            for y in -bound..=bound {
                if x * x + y * y <= bound * bound && (x != 0 || y != 0) {
                    sq_points.push((x, y));
                }
            }
        }
        let n_pts = sq_points.len();
        let mut sq_adj = vec![vec![0.0f64; n_pts]; n_pts];
        let mut sq_idx = std::collections::HashMap::new();
        for (i, &(x, y)) in sq_points.iter().enumerate() {
            sq_idx.insert((x, y), i);
        }
        let sq_neighbors = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)];
        for (i, &(x, y)) in sq_points.iter().enumerate() {
            for &(dx, dy) in &sq_neighbors {
                if let Some(&j) = sq_idx.get(&(x + dx, y + dy)) {
                    sq_adj[i][j] = 1.0;
                }
            }
        }
        let sq_cr = compute_cr(&sq_adj);

        (hex_cr, sq_cr)
    }
}
