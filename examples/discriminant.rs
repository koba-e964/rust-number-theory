extern crate num;
extern crate rust_number_theory;

use num::traits::Pow;
use num::BigInt;

use rust_number_theory::discriminant;
use rust_number_theory::polynomial::Polynomial;

fn main() {
    // 2x^3 + x^2 - 2x + 3
    let p: Polynomial<BigInt> =
        Polynomial::from_raw(vec![3.into(), (-2).into(), 1.into(), 2.into()]);
    let disc = discriminant::discriminant(&p);
    eprintln!("disc = {}", disc);
    assert_eq!(disc, (-1132).into());

    // x^3 + 9x + 1
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![1.into(), 9.into(), 0.into(), 1.into()]);
    let disc = discriminant::discriminant(&p);
    eprintln!("disc = {}", disc);
    assert_eq!(disc, (-2943).into());

    // 6x^5 - 7x^4 + 6x^3 - 7x^2 + 6x + 5
    let p: Polynomial<BigInt> = Polynomial::from_raw(
        vec![5, 6, -7, 6, -7, 6]
            .into_iter()
            .map(BigInt::from)
            .collect(),
    );
    let disc = discriminant::discriminant(&p);
    eprintln!("disc = {}", disc);
    assert_eq!(disc, BigInt::from(6).pow(4u64) * BigInt::from(7601837));
}
