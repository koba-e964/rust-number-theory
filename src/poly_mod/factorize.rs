#![allow(clippy::many_single_char_names, clippy::needless_range_loop, unused)]
use num::traits::{Num, NumAssign, NumOps, Zero};
use num::Integer;
use rand::distributions::uniform::SampleUniform;
use rand::{thread_rng, Rng};
use std::ops::Neg;

use crate::poly_mod::prim::{
    differential, divide_by_x_a, modinv, modpow, poly_divrem, poly_gcd, poly_mod, poly_mod_sub,
    poly_modpow, poly_of_mod,
};
use crate::polynomial::Polynomial;

/// If p is very large (so that p does not fit in usize), the parameter pusize is ignored.
/// In that case, the caller can pass any value.
pub fn squarefree<Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + From<i32>>(
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

/// Precondition: poly is a square-free polynomial mod p.
/// This function returns a vector of pairs (A_d, d),
/// where A_d is a product of distinct polynomials of degree d.
/// The returned array is sorted in d's ascending order.
fn degree<Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + From<i32>>(
    poly: &Polynomial<Int>,
    p: &Int,
) -> Vec<(Polynomial<Int>, usize)>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let x = Polynomial::from_raw(vec![Int::zero(), Int::one()]);
    let mut v = poly.clone();
    let mut w = x.clone();
    let mut d = 0;
    let mut result = vec![];
    while 2 * d + 2 <= v.deg() {
        d += 1;
        w = poly_modpow(&w, p, &v, p);
        let ad = poly_gcd(&poly_mod_sub(&w, &x, p), &v, p);
        if ad.deg() > 0 {
            result.push((ad.clone(), d));
            v = poly_divrem(&v, &ad, p).0;
            w = poly_divrem(&w, &v, p).1;
        }
    }
    if v.deg() > 0 {
        result.push((v.clone(), v.deg()));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // asserts a = c * b for some c in F(p)^\times.
    fn assert_associate(a: &Polynomial<i64>, b: &Polynomial<i64>, p: i64) {
        assert_eq!(a.deg(), b.deg());
        let al = a.coef_at(a.deg());
        let bl = b.coef_at(b.deg());
        let factor = al * modinv(&bl, &p) % p;
        assert_eq!(*a, b * &Polynomial::from_mono(factor));
    }

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
        assert_associate(&fac, &Polynomial::from_raw(vec![1, 1]), p);
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
        assert_associate(&fac, &Polynomial::from_raw(vec![1, 1]), p);
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
        assert_associate(&fac, &Polynomial::from_raw(vec![1, 1]), p);
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
        assert_associate(&fac, &Polynomial::from_raw(vec![2, 3, 1]), p);
    }

    #[test]
    fn degree_test_1() {
        let p = 3;
        let poly = Polynomial::from_raw(vec![2, 0, 0, 0, 1]);
        // distinct degree factorization of (x+1)(x+2)(x^2+1) is (x^2+2)(x^2+1)
        let result = degree::<i64>(&poly, &p);
        assert_eq!(result[0].1, 1);
        assert_associate(&result[0].0, &Polynomial::from_raw(vec![2, 0, 1]), p);
        assert_eq!(result[1].1, 2);
        assert_associate(&result[1].0, &Polynomial::from_raw(vec![1, 0, 1]), p);
    }

    #[test]
    fn degree_test_2() {
        let p = 3;
        let poly = Polynomial::from_raw(vec![1, 1, 0, 0, 1, 1]);
        // distinct degree factorization of (x+1)(x^4+1) is (x+1)(x^4+1)
        let result = degree::<i64>(&poly, &p);
        assert_eq!(result[0].1, 1);
        assert_associate(&result[0].0, &Polynomial::from_raw(vec![1, 1]), p);
        assert_eq!(result[1].1, 2);
        assert_associate(&result[1].0, &Polynomial::from_raw(vec![1, 0, 0, 0, 1]), p);
    }
}
