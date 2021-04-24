#![allow(clippy::many_single_char_names, clippy::needless_range_loop, unused)]
use num::traits::{Num, NumAssign, NumOps, Zero};
use num::Integer;
use rand::distributions::uniform::SampleUniform;
use rand::{thread_rng, Rng};
use std::ops::Neg;

use crate::poly_mod::prim::{
    differential, divide_by_x_a, modinv, modpow, poly_divrem, poly_gcd, poly_mod, poly_modpow,
    poly_of_mod,
};
use crate::polynomial::Polynomial;

/// If p is very large (so that p does not fit in usize), this field is ignored.
/// In that case, the caller can pass any value.
pub fn squarefree<
    Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + SampleUniform + From<i32>,
>(
    poly: &Polynomial<Int>,
    p: &Int,
    pusize: usize,
) -> Vec<(Polynomial<Int>, usize)>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    if poly.is_zero() {
        panic!();
    }
    let mut e = 1;
    let mut t0 = poly_mod(poly, p);
    let mut result = vec![];
    'outer: while t0.deg() != 0 {
        let der = differential(&t0, p);
        let mut t = poly_gcd::<Int>(&t0, &der, p);
        let mut v = poly_divrem(&t0, &t, p).0;
        let mut k = 0;
        loop {
            if v.deg() == 0 {
                let mut raw = vec![Int::zero(); t.deg() / pusize + 1];
                for i in 0..=t.deg() / pusize {
                    raw[i] = t.coef_at(pusize * i);
                }
                t0 = Polynomial::from_raw(raw);
                e *= pusize;
                continue 'outer;
            }
            k += 1;
            let w = poly_gcd::<Int>(&t, &v, p);
            let aek = poly_divrem(&v, &w, p).0;
            v = w;
            t = poly_divrem(&t, &v, p).0;
            if aek.deg() != 0 {
                result.push((aek, e * k));
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn squarefree_test_1() {
        let p = 5;
        let poly = Polynomial::from_raw(vec![1, 2, 1]);
        let ans = squarefree::<i64>(&poly, &p, 5);
        // squarefree factorization of x^2 + 2x + 1 is (x+1)^2.
        // Since we don't distinguish two polynomials equal upto scalar-multiplication,
        // we need to check fac ~ x + 1 rather than fac == x + 1.
        assert_eq!(ans.len(), 1);
        let &(ref fac, e) = &ans[0];
        assert_eq!(e, 2);
        assert_eq!(fac.deg(), 1);
        assert_eq!(fac.coef_at(0), fac.coef_at(1));
    }

    #[test]
    fn squarefree_test_2() {
        let p = 5;
        let poly = Polynomial::from_raw(vec![1, -2, 3, 1]);
        let ans = squarefree::<i64>(&poly, &p, 5);
        // squarefree factorization of x^3 + 3x^2 - 2x + 1 is (x+1)^3.
        // Since we don't distinguish two polynomials equal upto scalar-multiplication,
        // we need to check fac ~ x + 1 rather than fac == x + 1.
        assert_eq!(ans.len(), 1);
        let &(ref fac, e) = &ans[0];
        assert_eq!(e, 3);
        assert_eq!(fac.deg(), 1);
        assert_eq!(fac.coef_at(0), fac.coef_at(1));
    }

    #[test]
    fn squarefree_test_3() {
        let p = 5;
        let poly = Polynomial::from_raw(vec![1, 0, 0, 0, 0, 1]);
        let ans = squarefree::<i64>(&poly, &p, 5);
        // squarefree factorization of x^5 + 1 is (x+1)^5.
        // Since we don't distinguish two polynomials equal upto scalar-multiplication,
        // we need to check fac ~ x + 1 rather than fac == x + 1.
        assert_eq!(ans.len(), 1);
        let &(ref fac, e) = &ans[0];
        assert_eq!(e, 5);
        assert_eq!(fac.deg(), 1);
        assert_eq!(fac.coef_at(0), fac.coef_at(1));
    }

    #[test]
    fn squarefree_test_4() {
        let p = 5;
        let poly = Polynomial::from_raw(vec![2, 3, 1]);
        let ans = squarefree::<i64>(&poly, &p, 5);
        // squarefree factorization of x^2+3x+2 is x^2+3x+2.
        assert_eq!(ans.len(), 1);
        let &(ref fac, e) = &ans[0];
        assert_eq!(e, 1);
        assert_eq!(fac.deg(), 2);
    }
}
