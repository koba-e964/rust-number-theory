use std::convert::TryInto;

use num::{BigInt, BigRational, Zero};
use number_theory_linear::{hnf::HNF, subspace::image_mod_p};

use crate::{
    algebraic::Algebraic, ideal::Ideal, mult_table::MultTable, order::Order, poly_mod,
    polynomial::Polynomial,
};

// 6.2.2 of [Cohen]. Returns a list of pairs (P, e).
pub fn decompose<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
) -> Vec<(Ideal<'mul>, usize)> {
    let result = poly_mod::factorize_mod_p::<BigInt>(&theta.min_poly, p, p.try_into().unwrap_or(0));
    result
        .into_iter()
        .map(|(poly, mul)| {
            let poly = Polynomial::from_raw(
                poly.into_vec()
                    .into_iter()
                    .map(BigRational::from_integer)
                    .collect(),
            );
            let mut elem = if poly.deg() >= theta.min_poly.deg() {
                vec![BigInt::from(0); theta.deg()]
            } else {
                int_basis.to_z_basis_int(&Algebraic::with_expr(theta.min_poly.clone(), poly))
            };
            while elem.iter().any(|coef| !coef.is_zero())
                && elem.iter().all(|coef| coef % p == BigInt::zero())
            {
                for coef in &mut elem {
                    *coef /= p;
                }
            }
            let ancilla = Ideal::principal(&elem, mult_table);
            let mut pelem = vec![BigInt::from(0); theta.deg()];
            pelem[0] = p.clone();
            let pz = Ideal::principal(&pelem, mult_table);
            (&ancilla + &pz, mul)
        })
        .collect()
}

// 6.2.5 of [Cohen]. Multiplies two ideals I/pO and J/pO.
// TODO: A type for ideals over O/pO must be defined and used here.
#[allow(clippy::needless_range_loop)]
pub fn multiply<'mul>(
    _theta: &Algebraic,
    _int_basis: &Order,
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
