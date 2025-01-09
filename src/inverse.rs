use std::ops::{AddAssign, Div, Mul, Rem, Sub};

use num::{BigInt, Integer, One, Signed, Zero};

/// Perform extended gcd.
/// Returns (g, x, y) that satisfies g = gcd(a, b), g = xa + yb.
fn extgcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if false {
        extgcd_binary(a, b)
    } else {
        extgcd_division(a, b)
    }
}

#[allow(clippy::many_single_char_names)]
fn extgcd_division(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b == &BigInt::zero() {
        return (a.clone(), BigInt::one(), BigInt::zero());
    }
    let (q, r) = a.div_rem(b);
    let (g, x, y) = extgcd_division(b, &r);
    // bx + ry = g
    // bx + (a - bq)y = g
    // ya + (x - qy)b = g
    let newy = x - &q * &y;
    (g, y, newy)
}

fn extgcd_binary(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b == &BigInt::zero() {
        return (a.clone(), BigInt::one(), BigInt::zero());
    }
    if a == &BigInt::zero() {
        return (b.clone(), BigInt::zero(), BigInt::one());
    }
    if a < &BigInt::zero() {
        let (g, x, y) = extgcd_binary(&-a, b);
        return (g, -x, y);
    }
    if b < &BigInt::zero() {
        let (g, x, y) = extgcd_binary(a, &-b);
        return (g, x, -y);
    }
    let za = a.trailing_zeros().unwrap_or(0);
    let zb = b.trailing_zeros().unwrap_or(0);
    let z = za.min(zb);
    let a = a >> z;
    let b = b >> z;
    let (g, x, y) = if za <= zb {
        extgcd_1(&a, &b)
    } else {
        let (g, y, x) = extgcd_1(&b, &a);
        (g, x, y)
    };
    (g << z, x, y)
}

fn extgcd_1(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    debug_assert!(a.is_odd(), "{a}, {b}");
    if b.is_even() {
        let b1 = b >> 1;
        let (g, mut x, mut y) = extgcd_1(a, &b1);
        if y.is_even() {
            y >>= 1;
        } else {
            x -= b1;
            y += a;
            debug_assert!(y.is_even(), "{y}");
            y >>= 1;
        }
        return (g, x, y);
    }
    if a < b {
        let (g, x, y) = extgcd_1(b, a);
        return (g, y, x);
    }
    let ab = a - b;
    if ab == BigInt::zero() {
        return (b.clone(), BigInt::zero(), BigInt::one());
    }
    let (g, x, y) = extgcd_1(b, &ab);
    (g, y.clone(), x - &y)
}

/// Computes a^{-1} mod mo, or returns gcd(a, mo) if a and mo are not coprime.
/// mo should be positive.
pub fn inv(a: &BigInt, mo: &BigInt) -> Result<BigInt, BigInt> {
    let (g, x, _y) = extgcd(a, mo);
    if g.abs() != BigInt::one() {
        return Err(g.abs());
    }
    // g = \pm 1. Now xg = 1 (mod mo) holds.
    Ok(zmod::<BigInt>(&(&x * &g), mo))
}

/// Computes x % mo. The answer is always in [0, mo).
pub fn zmod<Int: Zero + Ord + for<'a> AddAssign<&'a Int>>(x: &Int, mo: &Int) -> Int
where
    for<'a> &'a Int: Mul<&'a Int, Output = Int>
        + Sub<&'a Int, Output = Int>
        + Div<&'a Int, Output = Int>
        + Rem<&'a Int, Output = Int>,
{
    let mut res = x % mo;
    if res < Int::zero() {
        res += mo;
    }
    res
}
