extern crate num;

pub mod polynomial;
pub mod resultant;
pub mod discriminant;

use num::BigInt;
use polynomial::Polynomial;
use resultant::resultant;

fn main() {
    // 9x^5 + 6x^4 + 2x^2 + 5
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()]);
    // 7x^4 + x^3 + 6x^2 + 6x + 6
    let q: Polynomial<BigInt> = Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
    let res = resultant(&p, &q);
    eprintln!("resultant = {}", res);
    assert_eq!(res, 335159672.into());
    // 2x^3 + x^2 - 2x + 3
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![3.into(), (-2).into(), 1.into(), 2.into()]);
    eprintln!("det = {}", discriminant::discriminant(&p));
}
