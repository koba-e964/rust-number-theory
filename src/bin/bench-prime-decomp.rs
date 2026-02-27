use std::hint::black_box;
use std::time::Instant;

use num::BigInt;
use rust_number_theory::{
    algebraic::Algebraic, integral_basis, polynomial::Polynomial, prime_decomp::decompose,
};

fn run_case(name: &str, poly: Polynomial<BigInt>, p: BigInt, iterations: usize) {
    let theta = Algebraic::new(poly);
    let int_basis = integral_basis::find_integral_basis(&theta);
    let mult_table = int_basis.get_mult_table(&theta);
    let mut samples_ns = Vec::with_capacity(iterations);

    for _ in 0..iterations {
        let start = Instant::now();
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        black_box(result);
        samples_ns.push(start.elapsed().as_nanos() as f64);
    }

    let count = samples_ns.len() as f64;
    let mean = samples_ns.iter().sum::<f64>() / count;
    let variance = if samples_ns.len() <= 1 {
        0.0
    } else {
        let sum_sq = samples_ns
            .iter()
            .map(|x| {
                let d = x - mean;
                d * d
            })
            .sum::<f64>();
        sum_sq / (count - 1.0)
    };
    let stddev = variance.sqrt();
    let total_ns = samples_ns.iter().sum::<f64>();
    println!("case={name}");
    println!("  iterations={iterations}");
    println!("  total_ms={:.3}", total_ns / 1_000_000.0);
    println!("  mean_us={:.3}", mean / 1_000.0);
    println!("  stddev_us={:.3}", stddev / 1_000.0);
}

fn main() {
    let iterations = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(200);

    // Q(sqrt(5)): index (Z_K : Z[sqrt(5)]) is 2.
    let sqrt5 = Polynomial::from_raw(vec![(-5).into(), 0.into(), 1.into()]);
    run_case("simple_path_sqrt5_p3", sqrt5.clone(), 3.into(), iterations);
    run_case("bl_path_sqrt5_p2", sqrt5, 2.into(), iterations);

    // A cubic field example used elsewhere in the repository.
    let cubic = Polynomial::from_raw(vec![(-2).into(), 0.into(), 0.into(), 1.into()]);
    run_case("cubic_path_p3", cubic, 3.into(), iterations);
}
