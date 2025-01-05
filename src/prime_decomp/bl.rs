use std::convert::TryInto;

use num::{BigInt, BigRational, Zero};
use number_theory_linear::{hnf::HNF, subspace::image_mod_p};

use crate::{
    algebraic::Algebraic,
    ideal::Ideal,
    mult_table::MultTable,
    order::{self, Order},
    poly_mod,
    polynomial::Polynomial,
};

// 6.2.2 of [Cohen]. Returns a list of pairs (P, e).
pub fn decompose<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
) -> Vec<(Ideal<'mul>, usize)> {
    panic!()
}

// 6.2.5 of [Cohen]. Multiplies two ideals I/pO and J/pO.
// TODO: A type for ideals over O/pO must be defined and used here.
#[allow(clippy::needless_range_loop)]
pub fn multiply<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
    i: &Ideal<'mul>,
    j: &Ideal<'mul>,
) -> Ideal<'mul> {
    let n = mult_table.deg();
    // 1. [Compute Matrix]
    let iv = i.as_hnf().as_vecs();
    let jv = j.as_hnf().as_vecs();
    assert_eq!(iv.len(), n);
    assert_eq!(jv.len(), n);
    let r = iv.len();
    let m = jv.len();
    let mut mat = vec![vec![BigInt::from(0); n]; r * m];
    for i in 0..r {
        for j in 0..m {
            let mut mul = mult_table.mul(&iv[i], &jv[j]);
            for k in 0..n {
                mat[i * m + j][k] = core::mem::take(&mut mul[k]) % p;
            }
        }
    }
    // 2. [Compute Image]
    let image = image_mod_p(&mat, p);
    let hnf = HNF::new(&image);
    Ideal::new(hnf, mult_table)
}
