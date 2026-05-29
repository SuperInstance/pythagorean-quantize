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

fn main() {
    println!("=== Pythagorean Quantize ===\n");

    let triples = pythagorean_triples::PythagoreanTriple::primitive(100);
    println!("Primitive triples with c <= 100: {}", triples.len());
    for t in triples.iter().take(5) {
        println!("  {}² + {}² = {}²", t.a, t.b, t.c);
    }

    let cr = pythagorean_triples::pythagorean_cr(200);
    println!("\nPythagorean CR (c<=200): {:.6}", cr);

    let (hex_cr, sq_cr) = eisenstein::EisensteinInt::hex_vs_square_cr(5);
    println!("\nHexagonal CR: {:.6}", hex_cr);
    println!("Square CR:    {:.6}", sq_cr);
    println!("Hex > Square: {}", hex_cr > sq_cr);

    let comparisons = lattice::LatticeComparison::compare_all(6);
    println!("\n=== Lattice Comparison ===");
    for lc in &comparisons {
        println!("{:>12}: CR={:.4} gap={:.4} quant={:.4}",
            lc.name, lc.cr, lc.spectral_gap, lc.quantization);
    }

    let (int_q, real_q) = eigenvalue::EigenvalueQuantization::integer_vs_real_quantization(6, 50);
    println!("\nInteger graph quantization: {:.6}", int_q);
    println!("Real graph quantization:    {:.6}", real_q);

    println!("\nDone.");
}
