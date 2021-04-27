#![allow(clippy::needless_range_loop)]
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
            product = poly_divrem(&poly_mod(&(&product * &current), &modulus), g, &modulus).1;
        }
        current = poly_divrem(&poly_mod(&(&current * &current), &modulus), g, modulus).1;
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
        raw[i] = f.coef_at(i) % p.clone();
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

pub fn poly_mod_sub<Int: Clone + NumAssign + Integer + Neg<Output = Int>>(
    a: &Polynomial<Int>,
    b: &Polynomial<Int>,
    p: &Int,
) -> Polynomial<Int> {
    poly_mod(&(a - b), p)
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
