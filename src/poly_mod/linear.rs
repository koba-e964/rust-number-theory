#![allow(clippy::needless_range_loop)]
use num::traits::{Num, NumAssign, NumOps};
use num::Integer;
use rand::distributions::uniform::SampleUniform;
use rand::{thread_rng, Rng};
use std::ops::Neg;

use crate::poly_mod::prim::{
    divide_by_x_a, modinv, modpow, poly_divrem, poly_gcd, poly_mod, poly_modpow, poly_of_mod,
};
use crate::polynomial::Polynomial;

/// Finds all linear factors of poly.
/// Precondition: p is a prime.
pub fn find_linear_factors<
    Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + SampleUniform,
>(
    poly: &Polynomial<Int>,
    p: Int,
) -> Vec<Int>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let mut result = vec![];
    let mut rng = thread_rng();
    let poly = poly_mod(poly, &p);
    let two = Int::one() + Int::one();
    if p == two {
        // we handle p = 2 case separately
        find_linear_factors_impl_mod2(&poly, &mut result);
        return result;
    }
    find_linear_factors_impl(&poly, p, &mut result, &mut rng);
    result
}

fn find_linear_factors_impl<
    Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + SampleUniform,
>(
    poly: &Polynomial<Int>,
    p: Int,
    result: &mut Vec<Int>,
    rng: &mut impl Rng,
) where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    if poly.deg() == 0 {
        return;
    }
    if poly.deg() == 1 {
        result.push((-poly.coef_at(0) * modinv(&poly.coef_at(1), &p)).mod_floor(&p));
        return;
    }
    let a = rng.gen_range(Int::zero()..p.clone());
    debug_assert!(modpow(&a, &p, &p) == a);
    let poly_orig = poly;
    let mut poly = poly.clone();
    if poly_of_mod(&poly, &a, &p) == Int::zero() {
        poly = divide_by_x_a(&poly, &a, &p);
        result.push(a.clone());
    }
    // gcd with (X + a)^{(p-1)/2} \pm 1
    let xa = Polynomial::from_raw(vec![(-a).mod_floor(&p), Int::one()]);
    let two = Int::one() + Int::one();
    let p1 = (p.clone() - Int::one()) / two;
    let xapow = poly_modpow(&xa, &p1, &poly, &p);
    let mut xapowp1 = xapow.clone() + Polynomial::from_mono(Int::one());
    if xapowp1.coef_at(0) >= p {
        xapowp1 = xapowp1 - Polynomial::from_mono(p.clone());
    }
    let gcd = poly_gcd(&xapowp1, &poly, &p);
    if gcd.deg() > 0 {
        let (quo, _) = poly_divrem(&poly, &gcd, &p);
        find_linear_factors_impl(&gcd, p.clone(), result, rng);
        poly = quo;
    }
    let mut xapowm1 = xapow + Polynomial::from_mono(p.clone() - Int::one());
    if xapowm1.coef_at(0) >= p {
        xapowm1 = xapowm1 - Polynomial::from_mono(p.clone());
    }
    let gcd = poly_gcd(&xapowm1, &poly, &p);
    if gcd.deg() > 0 {
        let (quo, _) = poly_divrem(&poly, &gcd, &p);
        find_linear_factors_impl(&gcd, p.clone(), result, rng);
        poly = quo;
    }
    // If poly == poly_orig, Since poly has no linear factors, we can simply discard poly.
    if poly_orig != &poly {
        find_linear_factors_impl(&poly, p, result, rng);
    }
}

fn find_linear_factors_impl_mod2<Int: Clone + Integer + NumAssign>(
    poly: &Polynomial<Int>,
    result: &mut Vec<Int>,
) {
    let vals = [Int::zero(), Int::one()];
    let two = Int::one() + Int::one();
    let mut poly = poly.clone();
    for i in 0..2 {
        while poly_of_mod(&poly, &vals[i], &two).is_zero() {
            poly = divide_by_x_a(&poly, &vals[i], &two);
            result.push(vals[i].clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_linear_factors_in_quadratic_1() {
        // X^2 + 1 mod 5
        let poly: Polynomial<i32> = Polynomial::from_raw(vec![1, 0, 1]);
        let p = 5;
        let mut factors = find_linear_factors::<i32>(&poly, p);
        factors.sort_unstable();
        assert_eq!(factors, vec![2, 3]);
    }

    #[test]
    fn find_linear_factors_in_quadratic_2() {
        // X^2 + 2 mod 5
        let poly: Polynomial<i32> = Polynomial::from_raw(vec![2, 0, 1]);
        let p = 5;
        let factors = find_linear_factors::<i32>(&poly, p);
        assert_eq!(factors, vec![]);
    }

    #[test]
    fn find_linear_factors_in_quadratic_3() {
        // An example found in Example 2.3 of https://www.cryptrec.go.jp/exreport/cryptrec-ex-0021-2001.pdf
        let poly = Polynomial::from_raw(vec![65696851, 38350500, -1304055, 1139835, 219113, 99535]);
        let p = 104743;
        let mut factors = find_linear_factors::<i64>(&poly, p);
        factors.sort_unstable();
        assert_eq!(factors, vec![15570, 20660, 69738]);
    }

    #[test]
    fn find_linear_factors_with_multiplicity_1() {
        // X^2 mod 23
        let poly: Polynomial<i32> = Polynomial::from_raw(vec![0, 0, 1]);
        let p = 23;
        let factors = find_linear_factors::<i32>(&poly, p);
        assert_eq!(factors, vec![0, 0]);
    }

    #[test]
    fn find_linear_factors_in_quadratic_mod_2_1() {
        // X^2 + 1 mod 2
        let poly = Polynomial::from_raw(vec![1, 0, 1]);
        let p = 2;
        let factors = find_linear_factors::<i32>(&poly, p);
        assert_eq!(factors, vec![1, 1]);
    }

    #[test]
    fn find_linear_factors_in_quadratic_mod_2_2() {
        // X^4 + X^3 + X^2 mod 2
        let poly = Polynomial::from_raw(vec![0, 0, 1, 1, 1]);
        let p = 2;
        let factors = find_linear_factors::<i32>(&poly, p);
        assert_eq!(factors, vec![0; 2]);
    }
}
