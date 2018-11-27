extern crate num;
pub mod polynomial;
pub mod resultant;
use num::BigInt;
use polynomial::{Polynomial, div_rem_bigint, div_rem_bigrational};

fn main() {
    // 9x^5 + 6x^4 + 2x^2 + 5
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()]);
    // 7x^4 + x^3 + 6x^2 + 6x + 6
    let q: Polynomial<BigInt> = Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
    let pr = Polynomial::from_raw(
        p.dat.into_iter().map(|x| x.into()).collect());
    let qr = Polynomial::from_raw(
        q.dat.into_iter().map(|x| x.into()).collect());
    eprintln!("pr % qr = {:?}", div_rem_bigrational(&pr, &qr));
}
