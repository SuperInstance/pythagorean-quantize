# pythagorean-quantize

**Eigenvalues of integer Laplacians snap to algebraic numbers. Pythagorean triples are the discrete spectrum of Euclidean distance. Hexagonal beats square.**

> **The aha moment:** When you build a graph Laplacian from integer weights, the eigenvalues don't form a messy continuum — they *snap* to algebraic numbers. Rational numbers, square roots of rationals, roots of low-degree polynomials with integer coefficients. The discrete structure of the integers imprints itself on the spectrum. Pythagorean triples — 3² + 4² = 5² — aren't just a curiosity of number theory. They're the *discrete spectrum of Euclidean distance*.

## Why This Exists

Continuous space has a continuous spectrum. Discrete space has a discrete spectrum. When you quantize space onto a lattice, the distances between points can only take certain values. For the integer lattice ℤ², those distances are √(a² + b²) for integers a, b — which means the possible squared distances are exactly the integers that can be written as sums of two squares.

Pythagorean triples are the special case where a² + b² = c² exactly — integer distances that close perfectly. This library explores what happens when you build graphs from these discrete distances and look at their spectral properties.

## Quick Start

```rust
use pythagorean_quantize::*;

fn main() {
    // Generate primitive Pythagorean triples with hypotenuse ≤ 50
    let triples = PythagoreanTriple::primitive(50);
    for t in &triples {
        println!("{}² + {}² = {}²", t.a, t.b, t.c);
    }
    // Output:
    // 3² + 4² = 5²
    // 5² + 12² = 13²
    // 8² + 15² = 17²
    // 7² + 24² = 25²
    // ...

    // Build a graph where triples are connected if they share a number
    let adj = PythagoreanTriple::triple_graph(&triples);
    let cr = compute_cr(&adj);
    println!("Conservation ratio of Pythagorean graph: {:.4}", cr);

    // Compare hexagonal vs square lattice
    let (hex_cr, sq_cr) = EisensteinInt::hex_vs_square_cr(3);
    println!("Hexagonal CR: {:.4}", hex_cr);
    println!("Square CR:    {:.4}", sq_cr);
    println!("Hex/Sq ratio: {:.4}", hex_cr / sq_cr);
}
```

## Modules

### Pythagorean Triples

The backbone. Euclid's formula generates all primitive triples:

```
a = m² - n²,  b = 2mn,  c = m² + n²
```

where m > n, gcd(m,n) = 1, and m,n are not both odd.

```rust
// Primitive triples only (coprime)
let primitive = PythagoreanTriple::primitive(100);

// All triples (including multiples like 6-8-10)
let all = PythagoreanTriple::all(100);

// How close is a real-valued triple to an integer one?
let dist = PythagoreanTriple::snap_distance(3.01, 3.99, 5.005);
// dist ≈ 0.016 — almost snapped to (3,4,5)
```

The **snap distance** is key: it measures how far a real-valued distance triple is from the nearest Pythagorean triple. Small snap distance = the real geometry "wants" to be discrete.

#### The Triple Graph

Connect two Pythagorean triples with an edge if they share a number. (3,4,5) connects to (5,12,13) because they share 5. This creates a graph whose structure encodes the multiplicative relationships between integer distances.

```rust
let triples = PythagoreanTriple::primitive(50);
let adj = PythagoreanTriple::triple_graph(&triples);
// adj[i][j] = 1.0 if triples i and j share a, b, or c
```

### Eisenstein Integers — The Hexagonal Lattice

Eisenstein integers are numbers of the form a + bω where ω = e^(2πi/6) = (-1 + i√3)/2. They form a **hexagonal lattice** — the densest 2D packing.

```rust
// The norm: |a + bω|² = a² - ab + b²
assert_eq!(EisensteinInt::norm(1, 0), 1);  // unit
assert_eq!(EisensteinInt::norm(1, 1), 1);  // another unit
assert_eq!(EisensteinInt::norm(2, 1), 3);  // 4 - 2 + 1 = 3

// Eisenstein integers with norm ≤ 10
let points = EisensteinInt::with_norm(10);
// Returns (a, b, norm) triples sorted by norm

// Build hexagonal lattice graph (6 neighbors)
let adj = EisensteinInt::lattice_graph(4);

// Compare hexagonal vs square lattice
let (hex_cr, sq_cr) = EisensteinInt::hex_vs_square_cr(3);
```

**Why hexagonal beats square:**
- Square lattice (ℤ²): 4 neighbors, coordination number 4
- Hexagonal lattice (Eisenstein): 6 neighbors, coordination number 6
- The hexagonal lattice is the optimal packing in 2D (Gauss proved this)
- Higher connectivity → better spectral properties → higher conservation ratio

### Eigenvalue Quantization

Build Laplacians from integer-weighted graphs and watch the eigenvalues snap to algebraic numbers:

```rust
// Random integer-weighted graph
let adj = EigenvalueQuantization::random_integer_graph(5, 3);
let eigs = EigenvalueQuantization::eigenvalues(&adj);

// How close are eigenvalues to rationals p/q?
let snap_dists = EigenvalueQuantization::rational_snap_distance(&eigs, 20);
// For integer graphs, snap distances are typically very small

// How close to roots of low-degree polynomials?
let alg_dists = EigenvalueQuantization::algebraic_snap_distance(&eigs, 3);

// Quantization quality: ratio of min gap to max gap
let q = EigenvalueQuantization::quantization_quality(&eigs);
// q near 1 = very uniform (quantized); q near 0 = irregular
```

#### Integer vs Real Quantization

The library directly compares what happens with integer vs. random real weights:

```rust
let (int_snap, real_snap) = EigenvalueQuantization::integer_vs_real_quantization(5, 30);
// int_snap: average distance to nearest rational (integer weights)
// real_snap: average distance to nearest rational (real weights)
// int_snap << real_snap — integer weights snap harder
```

This is the core finding: **discrete structure produces discrete spectra.** Integer weights → eigenvalues cluster near algebraic numbers. Random real weights → no such clustering.

### Lattice Comparison

Compare spectral properties across different lattice types:

```rust
let comparisons = LatticeComparison::compare_all(4);
for lc in &comparisons {
    println!("{:12} CR={:.4} gap={:.4} quant={:.4}",
        lc.name, lc.cr, lc.spectral_gap, lc.quantization);
}
```

Available lattices:
- **Square** (ℤ²): 4-connectivity
- **Hexagonal** (Eisenstein): 6-connectivity
- **Triangular**: 6-connectivity on a square grid (includes diagonals)
- **Penrose**: Fibonacci substitution chain (quasicrystal)
- **Fibonacci**: Nodes at Fibonacci positions

### Snap Tools

Practical utilities for working with quantized spectra:

```rust
// Snap eigenvalues to nearest algebraic numbers
let snapped = SnapTools::snap_eigenvalues(&eigs, 0.01);

// Get minimal polynomial fingerprints
let fingerprints = SnapTools::spectral_fingerprint(&eigs);
// e.g., ["x - 2 = 0", "x² - 3 = 0", "≈ 4.1234"]

// Check if a spectrum is quantized
let is_quant = SnapTools::is_quantized(&eigs, 0.5);

// CR after snapping eigenvalues
let snapped_cr = SnapTools::snapped_cr(&adj, 0.1);

// Find a graph whose spectrum approximates a target
let target = vec![1.0, 2.0, 3.0];
let adj = SnapTools::approximate_spectrum(&target, 3);
```

## The Conservation Ratio (CR)

Throughout this library, the **Conservation Ratio** measures spectral coherence:

```
CR = (Σ λᵢ)² / (n · Σ λᵢ²)
```

- CR = 1: all eigenvalues equal (maximally coherent)
- CR → 0: eigenvalues are spread out (incoherent)
- CR measures how "concentrated" the spectrum is

For lattices, higher CR means the graph acts more like a single resonant system. Hexagonal lattices consistently achieve higher CR than square lattices of similar size.

## The Big Picture

```
Euclidean distance on integers → Pythagorean triples
Pythagorean triples → graph (shared numbers = edges)
Graph → Laplacian → eigenvalues
Integer weights → eigenvalues snap to algebraic numbers
Eisenstein integers (hexagonal) → better connectivity → higher CR

Discrete space → discrete spectrum → quantization
```

The name "pythagorean-quantize" captures this chain: Pythagorean triples (integer distances) lead to spectral quantization (discrete eigenvalue structure).

## Honest Limitations

- **Power iteration eigenvalues are approximate.** We use power iteration with deflation for eigenvalue computation. It works well for the dominant eigenvalues but loses accuracy for smaller ones. For serious spectral analysis, use a proper eigensolver (e.g., `nalgebra` or `rust-ndarray` with `lapack`).
- **Graph sizes are limited.** The adjacency matrix is stored as `Vec<Vec<f64>>` — O(n²) memory. Don't try to build a 10,000-node graph.
- **The PRNG is deterministic but crude.** `EigenvalueQuantization::random_integer_graph` uses a simple LCG. It's deterministic (same input → same output) but not cryptographically random. Fine for exploration, not for Monte Carlo.
- **Snap distances use brute force search.** `PythagoreanTriple::snap_distance` generates all triples up to the ceiling of z and searches. Slow for large z.
- **CR is not the only spectral measure.** CR captures one aspect of spectral coherence. Spectral gap, algebraic connectivity, and eigenvector localization all matter too. This library focuses on CR for conceptual clarity.

## Installation

```toml
[dependencies]
pythagorean-quantize = "0.1.0"
```

Zero dependencies. Pure Rust.

## License

MIT
