use core::ops::Neg;
use num::traits::{Num, NumAssign, NumOps};
use num::Integer;

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
pub fn hensel_lift<Int: Clone + Integer + NumAssign + Num + Neg<Output = Int> + From<i32>>(
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
}
