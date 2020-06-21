extern crate num;
extern crate rust_number_theory;

use num::BigInt;

use rust_number_theory::polynomial::Polynomial;
use rust_number_theory::resultant::resultant;

fn main() {
    // 9x^5 + 6x^4 + 2x^2 + 5
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![
        5.into(),
        0.into(),
        2.into(),
        0.into(),
        6.into(),
        9.into(),
    ]);
    // 7x^4 + x^3 + 6x^2 + 6x + 6
    let q: Polynomial<BigInt> =
        Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
    let res = resultant(&p, &q);
    eprintln!("resultant = {}", res);
    assert_eq!(res, 335159672.into());

    // x^4 + x^2 + 2
    let p: Polynomial<BigInt> =
        Polynomial::from_raw(vec![2.into(), 0.into(), 1.into(), 0.into(), 1.into()]);
    // x^2 + 1
    let q: Polynomial<BigInt> = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into()]);
    let res = resultant(&p, &q);
    eprintln!("resultant = {}", res);
    assert_eq!(res, 4.into());
}
