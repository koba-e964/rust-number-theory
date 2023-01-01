use core::ops::Neg;
use num::traits::{NumAssign, NumOps};
use num::Integer;

use crate::poly_mod::prim::poly_coprime_witness;
use crate::polynomial::Polynomial;

use super::prim::{poly_div, poly_divrem, poly_mod, poly_mul};

/// Performs Hensel lift.
///
/// Algorithm 3.5.5 in [Cohen].
///
/// The following preconditions must be met:
/// - p, q: integers (not necessarily prime), r := gcd(p, q)
/// - c, a, b, u, v: polynomials
/// - c = ab (mod q), au + bv = 1 (mod p)
/// - gcd(l(a), r) = 1, deg(u) < deg(b), deg(v) < deg(a), deg(c) = deg(a) + deg(b)
///
/// This function returns a triple (a_1, b_1, qr) which satisfies c = a_1 b_1 (mod qr).
pub fn hensel_lift<Int: Clone + Integer + NumAssign + Neg<Output = Int> + From<i32>>(
    p: &Int,
    q: &Int,
    c: &Polynomial<Int>,
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    u: &Polynomial<Int>,
    v: &Polynomial<Int>,
) -> (Polynomial<Int>, Polynomial<Int>, Int)
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    // 1. Euclidean division
    let r = p.gcd(q);
    let f = poly_mod(&poly_div(&(c - &(a * b)), q), &r);
    let (t, _) = poly_divrem(&(v * &f), a, &r);
    // 2. Terminate
    let a0 = v * &f - a * &t; // in Z[X]
    let b0 = u * &f + b * &t; // in Z[X]
    let qr = q * &r;
    let a1 = poly_mod(&(a + &poly_mul(&a0, q)), &qr);
    let b1 = poly_mod(&(b + &poly_mul(&b0, q)), &qr);
    (a1, b1, qr)
}

// TODO: quadratic_hensel_lift
// Algorithm 3.5.6 in [Cohen].

fn hensel_lift_multiple<Int: Clone + Integer + NumAssign + Neg<Output = Int> + From<i32>>(
    p: &Int,
    q: &Int,
    c: &Polynomial<Int>,
    factors: &[Polynomial<Int>],
) -> (Vec<Polynomial<Int>>, Int)
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let n = factors.len();
    if n == 0 {
        let r = p.gcd(q);
        return (vec![], q * &r);
    }
    let mut accumulated = vec![];
    let mut current: Polynomial<Int> = Polynomial::from_mono(Int::one());
    for i in 0..n {
        current = poly_mod(&(&current * &factors[i]), q);
        accumulated.push(current.clone());
    }
    let mut result = vec![];
    let mut product = c.clone();
    for i in (1..n).rev() {
        let (u, v) = poly_coprime_witness(&accumulated[i - 1], &factors[i], p);
        let (a1, b1, _) = hensel_lift(p, q, &product, &accumulated[i - 1], &factors[i], &u, &v);
        result.push(b1);
        product = a1;
    }
    result.push(product);
    result.reverse();
    let r = p.gcd(q);
    (result, q * &r)
}

/// Lifts c's factorization mod p into mod p^e.
///
/// factors cannot have duplicate polynomials.
///
/// It is not necessary that c = \prod factors holds; c = (constant) * \prod factors is enough.
pub fn lift_factorization<Int: Clone + Integer + NumAssign + Neg<Output = Int> + From<i32>>(
    p: &Int,
    e: u32,
    c: &Polynomial<Int>,
    factors: &[Polynomial<Int>],
) -> Vec<Polynomial<Int>>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    // TODO: improve from naive implementation
    let mut res = factors.to_vec();
    let mut cur = p.clone();
    let lc = c.coef_at(c.deg());
    for _ in 1..e {
        let next_cur = &cur * p;
        let invlc = lc.extended_gcd(&next_cur).x.mod_floor(&next_cur);
        let divided = poly_mod(&poly_mul(c, &invlc), &next_cur);
        let (sub, _) = hensel_lift_multiple(p, &cur, &divided, &res);
        res = sub;
        cur = next_cur;
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hensel_lift_works_0() {
        // An example found in [Cohen].
        // C(X) = X^2 + 2X + 3, A(X) = X - 3, B(X) = X - 4
        // C(X) = A(X)B(X) mod 9
        let c = Polynomial::<i32>::from_raw(vec![3, 2, 1]);
        let a = Polynomial::<i32>::from_raw(vec![-3, 1]);
        let b = Polynomial::<i32>::from_raw(vec![-4, 1]);
        let u = Polynomial::<i32>::from_mono(1i32);
        let v = Polynomial::<i32>::from_mono(-1i32);
        let (a1, b1, qr) = hensel_lift::<i32>(&9, &9, &c, &a, &b, &u, &v);
        // (X+60)(X+23) = X^2 + 2X + 3 (mod 81) is found
        assert_eq!(a1, Polynomial::from_raw(vec![60, 1]));
        assert_eq!(b1, Polynomial::from_raw(vec![23, 1]));
        assert_eq!(qr, 81);
    }

    #[test]
    fn lift_factorization_works_0() {
        // lift the factorization of X^3 - 2 mod 5 (= (X+2)(X^2 + 3X + 4)) to mod 5^3
        let p = 5;
        let a = Polynomial::from_raw(vec![2, 1]);
        let b = Polynomial::from_raw(vec![4, 3, 1]);
        let result =
            lift_factorization::<i32>(&p, 3, &Polynomial::from_raw(vec![-2, 0, 0, 1]), &[a, b]);
        // (X + 72)(X^2 + 53X + 59)
        assert_eq!(result[0], Polynomial::from_raw(vec![72, 1]));
        assert_eq!(result[1], Polynomial::from_raw(vec![59, 53, 1]));
    }
}
