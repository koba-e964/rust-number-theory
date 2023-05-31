#![allow(clippy::needless_range_loop)]
use num::integer::ExtendedGcd;
use num::traits::{NumAssign, NumOps, Zero};
use num::Integer;
use std::ops::Neg;

use crate::polynomial::Polynomial;

pub fn modpow<Int: Clone + Integer>(x: &Int, e: &Int, modulus: &Int) -> Int {
    let mut e = e.clone();
    let mut product = Int::one();
    let mut current = x.clone();
    let two = Int::one() + Int::one();
    while e > Int::zero() {
        if e.is_odd() {
            product = product * current.clone() % modulus.clone();
        }
        current = current.clone() * current % modulus.clone();
        e = e.div_floor(&two);
    }
    product
}

pub fn modinv<Int: Clone + Integer>(x: &Int, p: &Int) -> Int {
    let two = Int::one() + Int::one();
    modpow(x, &(p.clone() - two), p)
}

// Very few type implements From<usize>. From<i32> is much more convenient.
pub fn differential<Int: Clone + NumAssign + Integer + From<i32>>(
    f: &Polynomial<Int>,
    p: &Int,
) -> Polynomial<Int> {
    if f.is_zero() {
        return Polynomial::zero();
    }
    let f = f.differential();
    poly_mod(&f, p)
}

pub fn poly_of_mod<Int: Clone + NumAssign + Integer>(f: &Polynomial<Int>, a: &Int, p: &Int) -> Int {
    let mut sum = Int::zero();
    let deg = f.deg();
    for i in (0..deg + 1).rev() {
        sum *= a.clone();
        sum += f.coef_at(i);
        sum %= p.clone();
    }
    sum
}

pub fn poly_modpow<Int: Clone + NumAssign + Integer>(
    x: &Polynomial<Int>,
    e: &Int,
    g: &Polynomial<Int>,
    modulus: &Int,
) -> Polynomial<Int>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let mut e = e.clone();
    let mut product: Polynomial<Int> = Polynomial::from_mono(Int::one());
    let mut current = x.clone();
    let two = Int::one() + Int::one();
    while e > Int::zero() {
        if e.is_odd() {
            product = poly_divrem(&poly_mod(&(&product * &current), modulus), g, modulus).1;
        }
        current = poly_divrem(&poly_mod(&(&current * &current), modulus), g, modulus).1;
        e = e.div_floor(&two);
    }
    product
}

pub fn poly_mod<Int: Clone + NumAssign + Integer>(f: &Polynomial<Int>, p: &Int) -> Polynomial<Int> {
    let deg = f.deg();
    if deg == usize::MAX {
        // f = 0
        return f.clone();
    }
    let mut raw = vec![Int::zero(); deg + 1];
    for i in 0..deg + 1 {
        raw[i] = f.coef_at(i).mod_floor(p);
    }
    Polynomial::from_raw(raw)
}

pub fn poly_div<Int: Clone + NumAssign + Integer>(
    f: &Polynomial<Int>,
    divisor: &Int,
) -> Polynomial<Int> {
    let deg = f.deg();
    if deg == usize::MAX {
        // f = 0
        return f.clone();
    }
    let mut raw = vec![Int::zero(); deg + 1];
    for i in 0..deg + 1 {
        raw[i] = f.coef_at(i).div_floor(divisor);
    }
    Polynomial::from_raw(raw)
}

pub fn poly_mul<Int: Clone + NumAssign + Integer>(
    f: &Polynomial<Int>,
    multiplier: &Int,
) -> Polynomial<Int> {
    let deg = f.deg();
    if deg == usize::MAX {
        // f = 0
        return f.clone();
    }
    let mut raw = vec![Int::zero(); deg + 1];
    for i in 0..deg + 1 {
        raw[i] = f.coef_at(i) * multiplier.clone();
    }
    Polynomial::from_raw(raw)
}

// Copy-pasted from polynomial.rs
// TODO unify
pub fn poly_divrem<Int: Clone + NumAssign + Integer>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> (Polynomial<Int>, Polynomial<Int>)
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let a_deg = a.deg();
    let b_deg = b.deg();
    let zero: Int = Int::zero();
    if a.is_zero() || b.is_zero() || a_deg < b_deg {
        return (Polynomial::from_mono(zero), a.clone());
    }
    assert!(!b.dat[b_deg].is_zero());
    let lc = &b.dat[b_deg];
    let invlc = modinv(lc, p);
    let mut tmp = a.dat.clone();
    let mut quo = vec![zero; a_deg - b_deg + 1];
    // Naive division
    for i in (0..a_deg - b_deg + 1).rev() {
        let coef = (&tmp[i + b_deg] * &invlc).mod_floor(p);
        for j in 0..b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
            tmp[i + j] = tmp[i + j].mod_floor(p);
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}

pub fn poly_gcd<Int: Clone + NumAssign + Integer>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> Polynomial<Int>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let (_, rem) = poly_divrem(a, b, p);
    if rem.is_zero() {
        return b.clone();
    }
    poly_gcd(b, &rem, p)
}

/// Returns (g, u, v) where g = au + bv.
///
/// p must be >= 2.
pub fn poly_ext_gcd<Int: Clone + NumAssign + Integer + Neg<Output = Int>>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> (Polynomial<Int>, Polynomial<Int>, Polynomial<Int>)
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let (quo, rem) = poly_divrem(a, b, p);
    if rem.is_zero() {
        return (b.clone(), rem, Polynomial::from_mono(Int::one()));
    }
    let (g, u0, v0) = poly_ext_gcd(b, &rem, p);
    // g = b * u0 + (a - b * quo) * v0
    let v = poly_mod_sub(&u0, &poly_mod(&(&quo * &v0), p), p);
    (g, v0, v)
}

pub fn poly_mod_sub<Int: Clone + NumAssign + Integer + Neg<Output = Int>>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> Polynomial<Int> {
    poly_mod(&(a - b), p)
}

/// If gcd(a, b) = 1, finds u, v s.t. au + bv = 1.
pub fn poly_coprime_witness<Int: Clone + NumAssign + Integer + Neg<Output = Int>>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> (Polynomial<Int>, Polynomial<Int>)
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    let (g, u, v) = poly_ext_gcd(a, b, p);
    if g.deg() != 0 {
        panic!("The gcd should be non-zero constant, but it isn't");
    }
    let ExtendedGcd { x: inv, .. } = g.coef_at(0).extended_gcd(p);
    let inv = inv.mod_floor(p);
    (
        poly_mod(&poly_mul(&u, &inv), p),
        poly_mod(&poly_mul(&v, &inv), p),
    )
}

pub fn divide_by_x_a<Int: Clone + NumAssign + Integer>(
    poly: &Polynomial<Int>,
    a: &Int,
    p: &Int,
) -> Polynomial<Int> {
    let deg = poly.deg();
    let mut coefs = vec![Int::zero(); deg];
    let mut carry = Int::zero();
    for i in (0..deg).rev() {
        carry += poly.coef_at(i + 1);
        carry = carry.mod_floor(p);
        coefs[i] = carry.clone();
        carry *= a.clone();
    }
    carry += poly.coef_at(0);
    carry = carry.mod_floor(p);
    debug_assert!(carry == Int::zero());
    Polynomial::from_raw(coefs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_coprime_witness_works_0() {
        let p = 5;
        let a = Polynomial::from_raw(vec![2, 1]);
        let b = Polynomial::from_raw(vec![3, 1]);
        let (u, v) = poly_coprime_witness::<i32>(&a, &b, &p);
        assert_eq!(u, Polynomial::from_mono(4));
        assert_eq!(v, Polynomial::from_mono(1));
    }

    #[test]
    fn poly_coprime_witness_works_1() {
        let p = 5;
        let a = Polynomial::from_raw(vec![1, 1, 1]);
        let b = Polynomial::from_raw(vec![3, 1]);
        let (u, v) = poly_coprime_witness::<i32>(&a, &b, &p);
        // (X^2 + X + 1) * 3 + (X + 3) * (2X + 1) = 1 (mod 5)
        assert_eq!(u, Polynomial::from_mono(3));
        assert_eq!(v, Polynomial::from_raw(vec![1, 2]));
    }

    #[test]
    fn poly_coprime_witness_works_2() {
        let p = 125; // prime power
        let a = Polynomial::from_raw(vec![1, 1, 1]);
        let b = Polynomial::from_raw(vec![3, 1]);
        let (u, v) = poly_coprime_witness::<i32>(&a, &b, &p);
        // (X^2 + X + 1) * 18 + (X + 3) * (107X + 36) = 1 (mod 5^3)
        assert_eq!(u, Polynomial::from_mono(18));
        assert_eq!(v, Polynomial::from_raw(vec![36, 107]));
    }
}
