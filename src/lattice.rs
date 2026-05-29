use crate::compute_cr;
use crate::eigenvalue::EigenvalueQuantization;

/// Compare different lattices spectrally
pub struct LatticeComparison {
    pub name: String,
    pub adj: Vec<Vec<f64>>,
    pub cr: f64,
    pub spectral_gap: f64,
    pub quantization: f64,
}

impl LatticeComparison {
    fn from_adj(name: &str, adj: Vec<Vec<f64>>) -> Self {
        let cr = compute_cr(&adj);
        let eigs = EigenvalueQuantization::eigenvalues(&adj);
        let mut sorted_eigs: Vec<f64> = eigs.iter().cloned().filter(|x| x.abs() > 1e-10).collect();
        sorted_eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let spectral_gap = if sorted_eigs.len() >= 2 {
            sorted_eigs[0]
        } else {
            0.0
        };
        let quantization = EigenvalueQuantization::quantization_quality(&eigs);
        Self {
            name: name.to_string(),
            adj,
            cr,
            spectral_gap,
            quantization,
        }
    }

    /// Square lattice (Z²) side x side with 4-connectivity
    pub fn square(side: usize) -> LatticeComparison {
        let n = side * side;
        let mut adj = vec![vec![0.0f64; n]; n];
        for r in 0..side {
            for c in 0..side {
                let idx = r * side + c;
                if c + 1 < side {
                    adj[idx][idx + 1] = 1.0;
                    adj[idx + 1][idx] = 1.0;
                }
                if r + 1 < side {
                    adj[idx][idx + side] = 1.0;
                    adj[idx + side][idx] = 1.0;
                }
            }
        }
        Self::from_adj("square", adj)
    }

    /// Hexagonal lattice (Eisenstein) with n_rings
    pub fn hexagonal(n_rings: usize) -> LatticeComparison {
        let max_norm = (n_rings as u64).pow(2);
        let adj = crate::eisenstein::EisensteinInt::lattice_graph(max_norm);
        Self::from_adj("hexagonal", adj)
    }

    /// Triangular lattice — 6-connectivity on square grid
    pub fn triangular(side: usize) -> LatticeComparison {
        let n = side * side;
        let mut adj = vec![vec![0.0f64; n]; n];
        for r in 0..side {
            for c in 0..side {
                let idx = r * side + c;
                // Right
                if c + 1 < side {
                    adj[idx][idx + 1] = 1.0;
                    adj[idx + 1][idx] = 1.0;
                }
                // Down
                if r + 1 < side {
                    adj[idx][idx + side] = 1.0;
                    adj[idx + side][idx] = 1.0;
                }
                // Diagonal (upper-right to lower-left)
                if r + 1 < side && c > 0 {
                    adj[idx][(r + 1) * side + c - 1] = 1.0;
                    adj[(r + 1) * side + c - 1][idx] = 1.0;
                }
            }
        }
        Self::from_adj("triangular", adj)
    }

    /// Fibonacci substitution lattice (1D Fibonacci chain — Penrose-like quasicrystal)
    pub fn penrose(inflation_levels: usize) -> LatticeComparison {
        // Fibonacci substitution: L → LS, S → L
        let mut chain = vec!['L'];
        for _ in 0..inflation_levels {
            let mut next = String::new();
            for ch in &chain {
                match ch {
                    'L' => next.push_str("LS"),
                    'S' => next.push('L'),
                    _ => {}
                }
            }
            chain = next.chars().collect();
        }

        let n = chain.len();
        if n < 2 {
            return Self::from_adj("penrose", vec![vec![0.0; 1]; 1]);
        }

        // Assign positions: L = golden ratio spacing, S = 1
        let phi = (1.0 + 5f64.sqrt()) / 2.0;
        let mut positions = vec![0.0f64];
        for &ch in &chain {
            let spacing = match ch {
                'L' => phi,
                'S' => 1.0,
                _ => 1.0,
            };
            positions.push(positions.last().unwrap() + spacing);
        }

        // Connect nearest neighbors
        let mut adj = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            if i + 1 < n {
                adj[i][i + 1] = 1.0;
                adj[i + 1][i] = 1.0;
            }
        }

        Self::from_adj("penrose", adj)
    }

    /// Fibonacci lattice (nodes at Fibonacci positions on a line)
    pub fn fibonacci(generations: usize) -> LatticeComparison {
        let mut fibs = vec![1u64, 1u64];
        for _ in 2..generations {
            fibs.push(fibs[fibs.len() - 1] + fibs[fibs.len() - 2]);
        }
        let n = fibs.len();
        let mut adj = vec![vec![0.0f64; n]; n];

        // Connect Fibonacci neighbors (consecutive Fibonacci numbers)
        for i in 0..n {
            if i + 1 < n {
                adj[i][i + 1] = 1.0;
                adj[i + 1][i] = 1.0;
            }
            // Also connect i to i+2 with weaker weight
            if i + 2 < n {
                adj[i][i + 2] = 0.5;
                adj[i + 2][i] = 0.5;
            }
        }

        Self::from_adj("fibonacci", adj)
    }

    /// Compare all lattices at given size parameter
    pub fn compare_all(size: usize) -> Vec<LatticeComparison> {
        vec![
            Self::square(size),
            Self::hexagonal(size),
            Self::triangular(size),
            Self::penrose(size),
            Self::fibonacci(size),
        ]
    }
}
