extern crate num;
extern crate rust_number_theory;

use num::traits::Pow;
use num::BigInt;
use std::str::FromStr;

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

    // GNFS164
    let p: Polynomial<BigInt> = Polynomial::from_raw(
        vec![
            -12171622290476241497444980012311021i128,
            -13601173202899548432935219131949,
            176917216602508818430161036,
            557524556427309931902111,
            5627796025215486707,
            8293702863045600,
        ]
        .into_iter()
        .map(BigInt::from)
        .collect(),
    );
    let disc = discriminant::discriminant(&p);
    eprintln!("disc = {}", disc);
    let small = BigInt::from(16u128 * 27 * 5 * 7 * 13 * 289 * 31 * 71 * 83 * 2927);
    let large1 = BigInt::from(50920368630756370485097u128);
    let large2 = BigInt::from(1106020951959326672306393u128);
    let p143 = BigInt::from_str(
        &("92490304539205889793185046492402267591607848264472".to_owned()
            + "97142254169034426859184808668955053078096124820133"
            + "9437358977285268760202133522385835558594971"),
    )
    .unwrap();
    assert_eq!(disc, -small * large1 * large2 * p143);
}
