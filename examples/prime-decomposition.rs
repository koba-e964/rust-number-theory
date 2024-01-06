extern crate num;
extern crate rust_number_theory;

use num::BigInt;
use rust_number_theory::{
    algebraic::Algebraic, integral_basis, polynomial::Polynomial, prime_decomp::decompose,
};

fn main() {
    // x^3 + 5x^2 + 6x - 3
    let poly: Polynomial<BigInt> =
        Polynomial::from_raw(vec![(-3).into(), 6.into(), 5.into(), 1.into()]);

    let p: BigInt = 3.into();
    let theta = Algebraic::new(poly);
    let int_basis = integral_basis::find_integral_basis(&theta);
    let mult_table = int_basis.get_mult_table(&theta);
    let result = decompose(&theta, &int_basis, &mult_table, &p);
    println!("p = {}", p);
    for (id, e) in result {
        println!("N(ideal) = {}, e = {}", id.norm(), e);
    }
}
