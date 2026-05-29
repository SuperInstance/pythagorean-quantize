# pythagorean-quantize — Discrete Spectra from Integer Lattices

Eigenvalues of integer Laplacians snap to algebraic numbers. Pythagorean triples are the discrete spectrum of Euclidean distance. Hexagonal beats square. Explore the spectral properties of quantized space.

## What This Gives You

- **Pythagorean triples** — generate primitives via Euclid's formula, build connection graphs
- **Conservation ratios** — spectral CR for Pythagorean graphs and integer lattices
- **Eisenstein integers** — hexagonal lattice vs. square lattice comparison
- **Eigenvalue snapping** — watch continuous spectra collapse to discrete algebraic numbers
- **Lattice quantization** — snap continuous-space graphs to integer coordinates

## Quick Start

```rust
use pythagorean_quantize::*;

// Generate primitive Pythagorean triples
let triples = PythagoreanTriple::primitive(50);
for t in &triples {
    println!("{}² + {}² = {}²", t.a, t.b, t.c);
}

// Build a connection graph and compute CR
let adj = PythagoreanTriple::triple_graph(&triples);
let cr = compute_cr(&adj);
println!("Conservation ratio: {:.4}", cr);

// Hexagonal vs square lattice
let (hex_cr, sq_cr) = EisensteinInt::hex_vs_square_cr(3);
println!("Hex CR: {:.4} | Square CR: {:.4} | Ratio: {:.4}", hex_cr, sq_cr, hex_cr / sq_cr);
```

## Modules

| Module | Description |
|--------|-------------|
| `pythagorean_triples` | Euclid's formula, primitive/all triple generation, connection graphs |
| `eisenstein` | Hexagonal lattice via Eisenstein integers, hex vs square CR comparison |
| `lattice` | Integer lattice quantization, coordinate snapping |
| `eigenvalue` | Laplacian eigenvalue computation for integer-weighted graphs |
| `snap` | Spectral snapping: continuous → discrete eigenvalue mapping |

## How It Fits

Part of the spectral graph theory stack. The conservation ratio computed here is the same CR used in [plato-room](https://github.com/SuperInstance/plato-room) for knowledge graphs and [plato-fleet](https://github.com/SuperInstance/plato-fleet) for fleet topology. The hexagonal-beats-square result motivates the Eisenstein lattice used throughout the ecosystem.

## Installation

```toml
[dependencies]
pythagorean-quantize = "0.1"
```

## License

MIT
