#![allow(clippy::needless_range_loop)]

extern crate num;
extern crate rust_number_theory;

use rust_number_theory::algebraic::Algebraic;
use rust_number_theory::integral_basis;
use rust_number_theory::polynomial::Polynomial;

fn main() {
    let p;
    let case_nr = 1;
    if case_nr == 1 {
        p = Polynomial::from_raw(vec![
            5.into(),
            6.into(),
            (-7).into(),
            6.into(),
            (-7).into(),
            6.into(),
        ]);
    } else if case_nr == 2 {
        p = Polynomial::from_raw(vec![1.into(), 9.into(), 0.into(), 1.into()]);
    } else if case_nr == 3 {
        p = Polynomial::from_raw(vec![37.into(), 2.into(), 1.into()]);
    } else {
        panic!();
    }
    let theta = Algebraic::new(p);
    find_integral_basis(&theta);
}

fn find_integral_basis(theta: &Algebraic) {
    let o = integral_basis::find_integral_basis(theta);
    eprintln!("Z_K = {}", o);

    eprintln!("D(Z_K) = {}", o.discriminant(theta));
}
