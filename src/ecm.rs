use num::bigint::RandBigInt;
use num::{BigInt, One, Signed, Zero};
use std::collections::HashMap;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, Rem, Sub};

use crate::prime;

/// Factorizes an integer.
/// This function calls other functions as subroutines.
pub fn factorize(x: &BigInt) -> Vec<(BigInt, u64)> {
    if x <= &BigInt::zero() {
        panic!("x <= 0: x = {}", x);
    }

    let b = 1000u64; // TODO

    let mut stack = vec![x.clone()];
    let mut map = HashMap::new();
    while let Some(now) = stack.pop() {
        if now <= BigInt::one() {
            continue;
        }
        if prime::is_prime(&now) {
            *map.entry(now).or_insert(0) += 1;
            continue;
        }
        let fac = ecm(
            &now,
            ECMConfig {
                b1: b,
                b2: 100 * b,
                verbose: false,
            },
        );
        let other = &now / &fac;
        stack.push(now);
        stack.push(other);
    }
    let mut result: Vec<(BigInt, u64)> = map.into_iter().collect();
    result.sort();
    result
}

/// Configuration for ECM.
#[allow(clippy::upper_case_acronyms)]
pub struct ECMConfig {
    pub b1: u64,
    pub b2: u64,
    pub verbose: bool,
}

/// Finds a factor.
#[allow(clippy::many_single_char_names)]
pub fn ecm(n: &BigInt, conf: ECMConfig) -> BigInt {
    debug_assert!(!prime::is_prime(n));

    let mut rng = rand::thread_rng();

    let mut count = 0u64;

    loop {
        count += 1;
        if conf.verbose {
            eprintln!("Trying curve {}, B1 = {}, B2 = {}", count, conf.b1, conf.b2);
        }
        // randomly pick a
        let a = rng.gen_bigint_range(&1.into(), n);
        let curve = Ell { a, n: n.clone() };

        // randomly select a point
        let x = rng.gen_bigint_range(&1.into(), n);
        let y = rng.gen_bigint_range(&1.into(), n);
        let pt = Point {
            x,
            y,
            z: BigInt::one(),
        };
        if let Err(fac) = ecm_oneshot(pt, curve, conf.b1, conf.b2) {
            if fac == BigInt::one() || &fac == n {
                continue;
            }
            debug_assert_eq!(n % &fac, BigInt::zero());
            return fac;
        }
    }
}

fn ecm_oneshot(mut pt: Point<BigInt>, curve: Ell<BigInt>, b1: u64, b2: u64) -> Result<(), BigInt> {
    for k in 1..b1 + 1 {
        pt = pt.mul(k.into(), &curve)?;
    }
    // Step 2: try all primes in range (b1, b2]
    let pt6 = pt.mul(6.into(), &curve)?;
    for &init in &[b1.saturating_sub(1) / 6 * 6 + 1, (b1 + 1) / 6 * 6 - 1] {
        let mut cur_e = init;
        let mut cur = pt.mul(init.into(), &curve)?;
        while cur_e <= b2 {
            cur = cur.add(&pt6, &curve)?;
            cur_e += 6;
        }
    }
    Ok(())
}

/// Projective coordinates
/// It is caller's responsibility to ensure Int does not overflow.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point<Int> {
    x: Int,
    y: Int,
    z: Int,
}

// y^2 = x^3 + ax + b (mod n)
// not (2 | n)
// In the actual computation, b is irrelevant.
#[derive(Debug)]
struct Ell<Int> {
    a: Int,
    n: Int,
}

impl<
        Int: Clone
            + Zero
            + One
            + Ord
            + for<'a> std::ops::AddAssign<&'a Int>
            + DivAssign<i64>
            + Sub<Output = Int>
            + Mul<i64, Output = Int>
            + Signed
            + Eq,
    > Point<Int>
where
    for<'a> &'a Int: Add<&'a Int, Output = Int>
        + Sub<&'a Int, Output = Int>
        + Mul<&'a Int, Output = Int>
        + Div<&'a Int, Output = Int>
        + Rem<&'a Int, Output = Int>
        + Rem<i64, Output = Int>
        + Mul<i64, Output = Int>,
{
    fn add(&self, other: &Self, curve: &Ell<Int>) -> Result<Self, Int> {
        if self.is_inf() {
            return Ok(other.clone());
        }
        if other.is_inf() {
            return Ok(self.clone());
        }
        let xdif = &self.x - &other.x;
        let n = &curve.n;
        if xdif == Int::zero() {
            if &(&self.y + &other.y) % n == Int::zero() {
                return Ok(Self::inf());
            }
            let lambda = &(&self.x * &self.x * 3) + &curve.a;
            let lambda = &lambda % n;
            let den = &(&self.y * 2) % n;
            let den2 = &(&den * &den) % n;
            let den3 = &(&den2 * &den) % n;
            let x3 = &lambda * &lambda - &(&self.x * 2) * &den2;
            let y3 = &lambda * &(&(&self.x * &den2) - &x3);
            let y3 = y3 - &self.y * &den3;
            return Self {
                x: zmod::<Int>(&(x3 * den), n),
                y: zmod::<Int>(&y3, n),
                z: den3,
            }
            .simplify(&curve);
        }
        let lambda = zmod::<Int>(&(&self.y - &other.y), n);
        let xdif2 = &(&xdif * &xdif) % n;
        let xdif3 = &(&xdif2 * &xdif) % n;
        let x3 = &lambda * &lambda - &(&self.x + &other.x) * &xdif2;
        let y3 = &lambda * &(&(&self.x * &xdif2) - &x3);
        let y3 = y3 - &self.y * &xdif3;
        Self {
            x: zmod::<Int>(&(x3 * xdif), n),
            y: zmod::<Int>(&y3, n),
            z: xdif3,
        }
        .simplify(curve)
    }
    fn mul(&self, mut e: Int, curve: &Ell<Int>) -> Result<Self, Int> {
        let mut sum = Self::inf();
        let mut cur = self.clone();
        while e > Int::zero() {
            if &e % 2 == Int::one() {
                sum = sum.add(&cur, curve)?;
            }
            cur = cur.add(&cur, curve)?;
            e /= 2;
        }
        Ok(sum)
    }
    fn inf() -> Self {
        Self {
            x: Int::zero(),
            y: Int::one(),
            z: Int::zero(),
        }
    }
    fn is_inf(&self) -> bool {
        self.z == Int::zero()
    }
    fn simplify(&self, curve: &Ell<Int>) -> Result<Self, Int> {
        if self.z == Int::zero() {
            return Ok(Self::inf());
        }
        let n = &curve.n;
        let invz = inv::<Int>(&self.z, n)?;
        let x = &(&self.x * &invz) % n;
        let y = &(&self.y * &invz) % n;
        Ok(Self {
            x,
            y,
            z: Int::one(),
        })
    }
}

/// Perform extended gcd.
/// Returns (g, x, y) that satisfies g = gcd(a, b), g = xa + yb.
#[allow(clippy::many_single_char_names)]
fn extgcd<Int: Zero + One + Eq + Sub<Output = Int> + Clone>(a: &Int, b: &Int) -> (Int, Int, Int)
where
    for<'a> &'a Int: Mul<&'a Int, Output = Int>
        + Sub<&'a Int, Output = Int>
        + Div<&'a Int, Output = Int>
        + Rem<&'a Int, Output = Int>,
{
    if b == &Int::zero() {
        return (a.clone(), Int::one(), Int::zero());
    }
    let q = a / b;
    let r = a - &(b * &q);
    let (g, x, y) = extgcd::<Int>(b, &r);
    // bx + ry = g
    // bx + (a - bq)y = g
    // ya + (x - qy)b = g
    let newy = x - &q * &y;
    (g, y, newy)
}

/// Computes a^{-1} mod mo, or returns gcd(a, mo) if a and mo are not coprime.
/// mo should be positive.
fn inv<
    Int: Zero + One + Eq + Ord + Sub<Output = Int> + Signed + for<'a> AddAssign<&'a Int> + Clone,
>(
    a: &Int,
    mo: &Int,
) -> Result<Int, Int>
where
    for<'a> &'a Int: Mul<&'a Int, Output = Int>
        + Sub<&'a Int, Output = Int>
        + Div<&'a Int, Output = Int>
        + Rem<&'a Int, Output = Int>,
{
    let (g, x, _y) = extgcd::<Int>(a, mo);
    if g.abs() != Int::one() {
        return Err(g.abs());
    }
    // g = \pm 1. Now xg = 1 (mod mo) holds.
    Ok(zmod::<Int>(&(x * g), mo))
}

/// Computes x % mo. The answer is always in [0, mo).
fn zmod<Int: Zero + Ord + for<'a> AddAssign<&'a Int>>(x: &Int, mo: &Int) -> Int
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_works_0() {
        let curve = Ell { a: 3, n: 5 };
        assert_eq!(
            Point { x: 2, y: 1, z: 1 }
                .add(&Point { x: 1, y: 4, z: 1 }, &curve)
                .unwrap(),
            Point { x: 1, y: 1, z: 1 }
        );
    }

    #[test]
    fn add_works_1() {
        // https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization#An_example
        let n = 455839i64;
        let curve = Ell { a: 5, n };
        assert_eq!(
            Point { x: 1, y: 1, z: 1 }
                .add(&Point { x: 1, y: 1, z: 1 }, &curve)
                .unwrap(),
            Point {
                x: 14,
                y: n - 53,
                z: 1
            }
        );
    }

    #[test]
    fn mul_works_0() {
        let curve = Ell { a: 3, n: 5 };
        let a = Point { x: 1, y: 1, z: 1 };
        assert!(a.mul(5, &curve).unwrap().is_inf());
    }

    #[test]
    fn mul_works_1() {
        // https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization#An_example
        let n = 455839i64;
        let curve = Ell { a: 5, n };
        let mut pt = Point { x: 1, y: 1, z: 1 };
        let mut i = 1;
        let factor = loop {
            pt = match pt.mul(i, &curve) {
                Ok(pt) => pt,
                Err(factor) => break factor,
            };
            i += 1;
        };
        assert_eq!(factor, 599);
    }

    #[test]
    fn ecm_works_0() {
        let a = BigInt::from(133);
        let factor = ecm(
            &a,
            ECMConfig {
                b1: 1000,
                b2: 100000,
                verbose: false,
            },
        );
        assert_eq!(a % factor, BigInt::zero());
    }
    #[test]
    fn ecm_works_1() {
        let large1 = BigInt::from(65_537u128);
        let large2 = BigInt::from(1_000_003u128);
        let a = large1 * large2;
        let factor = ecm(
            &a,
            ECMConfig {
                b1: 1000,
                b2: 100000,
                verbose: false,
            },
        );
        assert_eq!(a % factor, BigInt::zero());
    }

    // Too slow. This test takes about 16 sec in *release* build, about 1 min in debug build.
    #[allow(unused)]
    fn ecm_works_2() {
        let large1 = BigInt::from(1_000_000_007u128);
        let large2 = BigInt::from(1_000_000_009u128);
        let a = large1 * large2;
        let factor = ecm(
            &a,
            ECMConfig {
                b1: 1000,
                b2: 100000,
                verbose: false,
            },
        );
        eprintln!("factor={}", factor);
        assert_eq!(a % factor, BigInt::zero());
    }
}
