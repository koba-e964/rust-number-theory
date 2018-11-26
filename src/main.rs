extern crate num;
mod polynomial;
mod resultant;
use polynomial::{Polynomial, div_rem};

fn main() {
    // x^4 + x^2 + 1
    let p = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into(), 0.into(), 1.into()]);
    // x^2 + 2x + 3
    let q = Polynomial::from_raw(vec![3.into(), 2.into(), 1.into()]);
    eprintln!("p % q = {:?}", div_rem(&p, &q));
}
