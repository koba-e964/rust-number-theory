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

fn final_split<
    Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + From<i32> + SampleUniform,
>(
    poly: &Polynomial<Int>,
    p: &Int,
    d: usize,
) -> Vec<Polynomial<Int>>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let mut result = vec![];
    if p.is_odd() {
        let mut rng = thread_rng();
        final_split_odd(poly, p, d, &mut result, &mut rng);
    } else {
        final_split_2(poly, d, &mut result);
    }
    result
}

fn final_split_odd<
    Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + From<i32> + SampleUniform,
>(
    poly: &Polynomial<Int>,
    p: &Int,
    d: usize,
    result: &mut Vec<Polynomial<Int>>,
    rng: &mut impl Rng,
) where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let k = poly.deg() / d;
    if k == 0 {
        unreachable!();
    }
    if k == 1 {
        result.push(poly.clone());
        return;
    }
    loop {
        let mut poly_raw = vec![Int::zero(); 2 * d];
        for i in 0..2 * d {
            poly_raw[i] = rng.gen_range(Int::zero()..p.clone());
        }
        let t = Polynomial::from_raw(poly_raw);
        // Iterating O(d) times to create p^d is okay:
        // we need O(d) computation to handle poly anyway.
        let mut e = Int::one();
        for _ in 0..d {
            e *= p.clone();
        }
        e -= Int::one();
        e /= Int::one() + Int::one();
        let mut tpow = poly_modpow(&t, &e, poly, p);
        tpow = poly_mod_sub(&tpow, &Polynomial::from_mono(Int::one()), p);
        let b = poly_gcd(&tpow, poly, p);
        if b.is_zero() || b.deg() == 0 || b.deg() == poly.deg() {
            continue;
        }
        final_split_odd(&b, p, d, result, rng);
        let div = poly_divrem(poly, &b, p).0;
        final_split_odd(&div, p, d, result, rng);
        return;
    }
}

fn final_split_2<Int: Clone + Integer + NumAssign + Num>(
    poly: &Polynomial<Int>,
    d: usize,
    result: &mut Vec<Polynomial<Int>>,
) where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let two = Int::one() + Int::one();
    let k = poly.deg() / d;
    if k == 0 {
        unreachable!();
    }
    if k == 1 {
        result.push(poly.clone());
        return;
    }
    let mut t = Polynomial::from_raw(vec![Int::zero(), Int::one()]);
    let x2 = Polynomial::from_raw(vec![Int::zero(), Int::zero(), Int::one()]);
    loop {
        let mut c = t.clone();
        for _ in 0..d - 1 {
            c = &(&c * &c) + &t;
            c = poly_mod(&c, &two);
            c = poly_divrem(&c, &poly, &two).1;
        }
        let b = poly_gcd(poly, &c, &two);
        if b.deg() == 0 || b.deg() == poly.deg() {
            t = &t * &x2;
            continue;
        }
        final_split_2(&b, d, result);
        let div = poly_divrem(poly, &b, &two).0;
        final_split_2(&div, d, result);
        return;
    }
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

    #[test]
    fn final_split_odd_test_1() {
        let p = 3;
        let poly = Polynomial::from_raw(vec![2, 0, 1]);
        // 2+x^2=(1+x)(2+x)
        let result = final_split::<i64>(&poly, &p, 1);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn final_split_odd_test_2() {
        let p = 3;
        let poly = Polynomial::from_raw(vec![0, 2, 0, 1]);
        // 2x+x^3=x(1+x)(2+x)
        let result = final_split::<i64>(&poly, &p, 1);
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn final_split_odd_test_3() {
        let p = 3;
        let poly = Polynomial::from_raw(vec![1, 0, 0, 0, 1]);
        // 1+x^4=(2+x+x^2)(2+2x+x^2)
        let result = final_split::<i64>(&poly, &p, 2);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn final_split_2_test_1() {
        let p = 2;
        let poly = Polynomial::from_raw(vec![1, 1, 1, 1, 1, 1, 1]);
        // 1 + x + ... + x^6 = (1+x+x^3)(1+x^2+x^3)
        let result = final_split::<i64>(&poly, &p, 3);
        assert_eq!(result.len(), 2);
    }
}
