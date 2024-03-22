use super::bitvector::BitVector;

/// print `bv` to stdout in 16 groups of 4 per row
pub fn print_bit_vec(bv: &[usize]) {
    for line in bv.chunks(16 * 4) {
        for chunk in line.chunks(4) {
            print!(" ");
            for elt in chunk {
                print!("{elt}");
            }
        }
        println!();
    }
}

/// return the number of bits in the intersection of a and b
pub fn intersect(a: &BitVector, b: &BitVector) -> usize {
    a.iter()
        .zip(b.iter())
        .fold(0, |acc, (a, b)| acc + (a & b).count_ones()) as usize
}

/// Computes the Tanimoto distance between bit vectors a and b
///
/// T(a, b) = (a ∩ b) / (a + b - a ∩ b), at least according to
/// featurebase.com/blog/tanimoto-and-chemical-similarity-in-featurebase
pub fn tanimoto(a: &BitVector, b: &BitVector) -> f64 {
    let num = intersect(a, b);
    let den: usize = a.count() + b.count() - num;
    num as f64 / den as f64
}
