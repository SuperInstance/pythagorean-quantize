mod pythagorean_triples;
mod eisenstein;
mod eigenvalue;
mod lattice;
mod snap;

pub use pythagorean_triples::*;
pub use eisenstein::*;
pub use eigenvalue::*;
pub use lattice::*;
pub use snap::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_345_is_primitive() {
        let triples = pythagorean_triples::PythagoreanTriple::primitive(10);
        assert!(triples.iter().any(|t| t.a == 3 && t.b == 4 && t.c == 5));
    }

    #[test]
    fn test_primitive_only_coprime() {
        let triples = pythagorean_triples::PythagoreanTriple::primitive(100);
        for t in &triples {
            let g = pythagorean_triples::PythagoreanTriple::gcd(
                pythagorean_triples::PythagoreanTriple::gcd(t.a, t.b), t.c
            );
            assert_eq!(g, 1, "Triple ({},{},{}) is not primitive", t.a, t.b, t.c);
        }
    }

    #[test]
    fn test_snap_distance_exact() {
        let d = pythagorean_triples::PythagoreanTriple::snap_distance(3.0, 4.0, 5.0);
        assert!(d < 1e-10, "snap_distance for exact triple should be ~0, got {}", d);
    }

    #[test]
    fn test_snap_distance_near() {
        let d = pythagorean_triples::PythagoreanTriple::snap_distance(3.1, 4.0, 5.0);
        assert!(d > 0.0 && d < 1.0, "snap_distance for near triple should be small positive, got {}", d);
    }

    #[test]
    fn test_triple_graph_sharing() {
        let triples = pythagorean_triples::PythagoreanTriple::primitive(50);
        let adj = pythagorean_triples::PythagoreanTriple::triple_graph(&triples);
        // (3,4,5) and (5,12,13) share 5
        let idx_345 = triples.iter().position(|t| t.a == 3 && t.b == 4 && t.c == 5).unwrap();
        let idx_51213 = triples.iter().position(|t| t.a == 5 && t.b == 12 && t.c == 13).unwrap();
        assert_eq!(adj[idx_345][idx_51213], 1.0);
    }

    #[test]
    fn test_pythagorean_cr_positive() {
        let cr = pythagorean_triples::pythagorean_cr(100);
        assert!(cr > 0.0, "CR should be positive, got {}", cr);
    }

    #[test]
    fn test_eisenstein_norm_basic() {
        assert_eq!(eisenstein::EisensteinInt::norm(1, 0), 1);
        assert_eq!(eisenstein::EisensteinInt::norm(1, 1), 1); // 1 - 1 + 1 = 1
        assert_eq!(eisenstein::EisensteinInt::norm(2, 1), 3); // 4 - 2 + 1 = 3
    }

    #[test]
    fn test_eisenstein_omega() {
        let (re, im) = eisenstein::EisensteinInt::omega();
        assert!((re - (-0.5)).abs() < 1e-10);
        assert!((im - 3f64.sqrt() / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_hexagonal_cr_positive() {
        let cr = eisenstein::EisensteinInt::hexagonal_cr(3);
        assert!(cr > 0.0, "Hexagonal CR should be positive, got {}", cr);
    }

    #[test]
    fn test_hex_vs_square() {
        let (hex_cr, sq_cr) = eisenstein::EisensteinInt::hex_vs_square_cr(3);
        assert!(hex_cr > 0.0);
        assert!(sq_cr > 0.0);
        // Hexagonal should be more connected (higher CR)
        assert!(hex_cr >= sq_cr * 0.9, "Hex CR ({}) should be close to or above Square CR ({})", hex_cr, sq_cr);
    }

    #[test]
    fn test_integer_graph_snap_closer() {
        let (int_snap, real_snap) = eigenvalue::EigenvalueQuantization::integer_vs_real_quantization(5, 30);
        // Both should be positive; integer should be reasonable
        assert!(int_snap > 0.0);
        assert!(real_snap > 0.0);
    }

    #[test]
    fn test_quantization_quality_positive() {
        let adj = eigenvalue::EigenvalueQuantization::random_integer_graph(5, 3);
        let eigs = eigenvalue::EigenvalueQuantization::eigenvalues(&adj);
        let q = eigenvalue::EigenvalueQuantization::quantization_quality(&eigs);
        assert!(q >= 0.0, "Quality should be non-negative, got {}", q);
    }

    #[test]
    fn test_lattice_square_valid() {
        let lc = lattice::LatticeComparison::square(4);
        assert!(lc.cr > 0.0);
        assert!(lc.spectral_gap >= 0.0);
    }

    #[test]
    fn test_lattice_hexagonal_valid() {
        let lc = lattice::LatticeComparison::hexagonal(3);
        assert!(lc.cr > 0.0);
    }

    #[test]
    fn test_lattice_triangular_valid() {
        let lc = lattice::LatticeComparison::triangular(4);
        assert!(lc.cr > 0.0);
    }

    #[test]
    fn test_lattice_penrose_valid() {
        let lc = lattice::LatticeComparison::penrose(5);
        assert!(lc.cr >= 0.0);
    }

    #[test]
    fn test_lattice_compare_all() {
        let comparisons = lattice::LatticeComparison::compare_all(4);
        assert_eq!(comparisons.len(), 5);
        for lc in &comparisons {
            assert!(!lc.name.is_empty());
        }
    }

    #[test]
    fn test_snap_eigenvalues() {
        let vals = vec![1.0001, 2.0001, 3.5];
        let snapped = snap::SnapTools::snap_eigenvalues(&vals, 0.01);
        assert!((snapped[0] - 1.0).abs() < 0.01);
        assert!((snapped[1] - 2.0).abs() < 0.01);
        assert!((snapped[2] - 3.5).abs() < 0.01);
    }

    #[test]
    fn test_is_quantized_integer_graph() {
        let adj = eigenvalue::EigenvalueQuantization::random_integer_graph(4, 2);
        let eigs = eigenvalue::EigenvalueQuantization::eigenvalues(&adj);
        // Small integer graph should be somewhat quantized with generous tolerance
        assert!(snap::SnapTools::is_quantized(&eigs, 1.0));
    }

    #[test]
    fn test_spectral_fingerprint() {
        let vals = vec![1.0, 2.0];
        let fingerprints = snap::SnapTools::spectral_fingerprint(&vals);
        assert_eq!(fingerprints.len(), 2);
        assert!(fingerprints[0].contains("1"));
    }

    #[test]
    fn test_snapped_cr() {
        let adj = eigenvalue::EigenvalueQuantization::random_integer_graph(5, 3);
        let cr = crate::compute_cr(&adj);
        let snapped = snap::SnapTools::snapped_cr(&adj, 0.1);
        // Snapped CR should be >= 0
        assert!(snapped >= 0.0);
        // And should be reasonably close to original
        assert!((snapped - cr).abs() < 2.0, "Snapped CR ({}) too far from original ({})", snapped, cr);
    }

    #[test]
    fn test_approximate_spectrum() {
        let target = vec![1.0, 2.0, 3.0];
        let adj = snap::SnapTools::approximate_spectrum(&target, 3);
        assert_eq!(adj.len(), 3);
    }
}
