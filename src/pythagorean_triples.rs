/// Pythagorean triples and their graph structure
#[derive(Debug, Clone)]
pub struct PythagoreanTriple {
    pub a: u64,
    pub b: u64,
    pub c: u64,
}

impl PythagoreanTriple {
    fn new(a: u64, b: u64, c: u64) -> Self {
        // Ensure a <= b for canonical form
        if a <= b {
            Self { a, b, c }
        } else {
            Self { a: b, b: a, c }
        }
    }

    pub fn gcd(a: u64, b: u64) -> u64 {
        let mut a = a;
        let mut b = b;
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    /// Generate primitive triples using Euclid's formula:
    /// a = m²-n², b = 2mn, c = m²+n² for coprime m>n, not both odd
    pub fn primitive(max_c: u64) -> Vec<PythagoreanTriple> {
        let mut triples = Vec::new();
        let m_max = ((max_c as f64).sqrt() as u64) + 1;
        for m in 2..=m_max {
            for n in 1..m {
                if (m - n) % 2 == 0 {
                    continue; // both odd or both even
                }
                if Self::gcd(m, n) != 1 {
                    continue;
                }
                let a = m * m - n * n;
                let b = 2 * m * n;
                let c = m * m + n * n;
                if c > max_c {
                    break;
                }
                triples.push(Self::new(a, b, c));
            }
        }
        triples.sort_by_key(|t| t.c);
        triples
    }

    /// All triples up to max_c (primitive * k)
    pub fn all(max_c: u64) -> Vec<PythagoreanTriple> {
        let primitives = Self::primitive(max_c);
        let mut all = Vec::new();
        for prim in &primitives {
            let mut k = 1u64;
            loop {
                let c = prim.c * k;
                if c > max_c {
                    break;
                }
                all.push(Self::new(prim.a * k, prim.b * k, c));
                k += 1;
            }
        }
        all.sort_by_key(|t| t.c);
        all
    }

    /// The "snap distance" — how close is a real triple (x,y,z) to an integer one?
    pub fn snap_distance(x: f64, y: f64, z: f64) -> f64 {
        // Find nearest integer triple that satisfies a² + b² = c²
        let ix = x.round() as u64;
        let iy = y.round() as u64;
        let iz = z.round() as u64;

        // Check if rounded values form a valid triple
        if ix * ix + iy * iy == iz * iz {
            return ((x - ix as f64).powi(2)
                + (y - iy as f64).powi(2)
                + (z - iz as f64).powi(2))
            .sqrt();
        }

        // Search nearby triples
        let max_c = (z.ceil() as u64) + 2;
        let triples = Self::all(max_c.max(5));
        let mut min_dist = f64::MAX;
        for t in &triples {
            let d = ((x - t.a as f64).powi(2)
                + (y - t.b as f64).powi(2)
                + (z - t.c as f64).powi(2))
            .sqrt();
            if d < min_dist {
                min_dist = d;
            }
            // Also check swapped a,b
            let d2 = ((x - t.b as f64).powi(2)
                + (y - t.a as f64).powi(2)
                + (z - t.c as f64).powi(2))
            .sqrt();
            if d2 < min_dist {
                min_dist = d2;
            }
        }
        if min_dist == f64::MAX { 0.0 } else { min_dist }
    }

    /// Build adjacency matrix of triples where edges connect triples that share a number
    pub fn triple_graph(triples: &[PythagoreanTriple]) -> Vec<Vec<f64>> {
        let n = triples.len();
        let mut adj = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let ti = [&triples[i].a, &triples[i].b, &triples[i].c];
                let tj = [&triples[j].a, &triples[j].b, &triples[j].c];
                let shares = ti.iter().any(|x| tj.iter().any(|y| x == y));
                if shares {
                    adj[i][j] = 1.0;
                    adj[j][i] = 1.0;
                }
            }
        }
        adj
    }
}

/// Conservation ratio of the Pythagorean triple graph
pub fn pythagorean_cr(max_c: u64) -> f64 {
    let triples = PythagoreanTriple::all(max_c);
    if triples.len() < 2 {
        return 0.0;
    }
    let adj = PythagoreanTriple::triple_graph(&triples);
    compute_cr(&adj)
}

/// Compute conservation ratio: ratio of eigenvalue sum² to sum of squares
/// Measures how "coherent" the graph spectrum is
pub fn compute_cr(adj: &[Vec<f64>]) -> f64 {
    let n = adj.len();
    if n == 0 {
        return 0.0;
    }

    // Compute degree matrix
    let degrees: Vec<f64> = adj
        .iter()
        .map(|row| row.iter().sum())
        .collect();

    // Compute Laplacian: L = D - A
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

    let eigenvalues = power_iteration_eigenvalues(&lap, n.min(20));

    if eigenvalues.is_empty() {
        return 0.0;
    }

    let sum: f64 = eigenvalues.iter().sum();
    let sum_sq: f64 = eigenvalues.iter().map(|x| x * x).sum();

    if sum_sq.abs() < 1e-12 {
        return 0.0;
    }

    // CR = (sum²) / (n * sum_sq) — normalized coherence
    // For perfectly coherent spectrum, CR approaches 1
    (sum * sum) / (eigenvalues.len() as f64 * sum_sq)
}

/// Power iteration method for eigenvalue computation
pub fn power_iteration_eigenvalues(mat: &[Vec<f64>], k: usize) -> Vec<f64> {
    let n = mat.len();
    if n == 0 {
        return vec![];
    }
    let k = k.min(n);

    let mut result = Vec::with_capacity(k);
    let mut deflated = mat.to_vec();

    for _ in 0..k {
        let mut v: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) / n as f64).collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in v.iter_mut() {
            *x /= norm;
        }

        let mut eigenvalue = 0.0f64;
        for _ in 0..200 {
            // Multiply
            let mut new_v = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    new_v[i] += deflated[i][j] * v[j];
                }
            }

            let norm = new_v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                break;
            }
            eigenvalue = 0.0;
            for i in 0..n {
                eigenvalue += v[i] * new_v[i];
            }
            for x in new_v.iter_mut() {
                *x /= norm;
            }
            v = new_v;
        }

        if eigenvalue.abs() < 1e-12 {
            break;
        }
        result.push(eigenvalue);

        // Deflate
        for i in 0..n {
            for j in 0..n {
                deflated[i][j] -= eigenvalue * v[i] * v[j];
            }
        }
    }

    result
}
