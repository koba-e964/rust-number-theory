use num::bigint::RandBigInt;
use num::{BigInt, One, Signed, Zero};
use std::collections::HashMap;
use std::ops::{AddAssign, Div, Mul, Rem, Sub};

use crate::perfect_power::perfect_power;
use crate::prime;

pub struct EcmStats {
    pub curve_count: u64,
}

/// Factorizes an integer.
/// This function calls other functions as subroutines.
pub fn factorize(x: &BigInt) -> Vec<(BigInt, u64)> {
    factorize_verbose(x, false).0
}
pub fn factorize_verbose(x: &BigInt, verbose: bool) -> (Vec<(BigInt, u64)>, EcmStats) {
    if x <= &BigInt::zero() {
        panic!("x <= 0: x = {}", x);
    }

    let b = select_b(x);

    let mut stack = vec![(x.clone(), 1)];
    let mut map = HashMap::new();
    let mut count = 0;
    while let Some((now, multiplicity)) = stack.pop() {
        if now <= BigInt::one() {
            continue;
        }
        if prime::is_prime(&now) {
            *map.entry(now).or_insert(0) += multiplicity;
            continue;
        }
        {
            let (b, k) = perfect_power(&now);
            if k >= 2 {
                stack.push((b, multiplicity * k as u64));
                continue;
            }
        }
        let (fac, nowcount) = ecm(
            &now,
            ECMConfig {
                b1: b,
                b2: 100 * b,
                verbose,
            },
        );
        count += nowcount;
        if fac == BigInt::one() {
            stack.push((now, multiplicity));
            continue;
        }
        let other = &now / &fac;
        stack.push((fac, multiplicity));
        stack.push((other, multiplicity));
    }
    let mut result: Vec<(BigInt, u64)> = map.into_iter().collect();
    result.sort();
    (result, EcmStats { curve_count: count })
}

/// Select appropriate B1.
fn select_b(n: &BigInt) -> u64 {
    if n <= &BigInt::from(1000u64) {
        return 4;
    }
    let lnx = n.bits() as f64 * 2.0f64.ln() / 2.0;
    let lnlnx = lnx.ln();
    let b = (lnx * lnlnx / 2.0).sqrt().exp(); // L(p)^{1/sqrt(2)}
    b as u64
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
pub fn ecm(n: &BigInt, conf: ECMConfig) -> (BigInt, u64) {
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
            if conf.verbose {
                eprintln!("Found factor after {count} trials");
            }
            return (fac, count);
        }
    }
}

fn ecm_oneshot(mut pt: Point, curve: Ell, b1: u64, b2: u64) -> Result<(), BigInt> {
    for k in 1..b1 + 1 {
        pt = pt.mul(k.into(), &curve)?;
        if pt.is_inf() {
            return Ok(());
        }
    }
    // Step 2: try all primes in range (b1, b2]
    for &init in &[b1.saturating_sub(1) / 6 * 6 + 1, (b1 + 1) / 6 * 6 - 1] {
        let mut cur_e = init;
        let p6 = {
            let p2 = pt.add(&pt, &curve)?;
            let p4 = p2.add(&p2, &curve)?;
            p2.add(&p4, &curve)?
        };
        pt = pt.mul(init.into(), &curve)?;
        if pt.is_inf() {
            return Ok(());
        }
        while cur_e <= b2 {
            cur_e += 6;
            pt = pt.add(&p6, &curve)?;
            if pt.is_inf() {
                return Ok(());
            }
        }
    }
    Ok(())
}

/// Projective coordinates
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point {
    x: BigInt,
    y: BigInt,
    z: BigInt,
}

// y^2 = x^3 + ax + b (mod n)
// not (2 | n)
// In the actual computation, b is irrelevant.
#[derive(Debug)]
struct Ell {
    a: BigInt,
    n: BigInt,
}

impl Point {
    fn add(&self, other: &Self, curve: &Ell) -> Result<Self, BigInt> {
        if self.is_inf() {
            return Ok(other.clone());
        }
        if other.is_inf() {
            return Ok(self.clone());
        }
        let xdif = &self.x - &other.x;
        let n = &curve.n;
        if xdif == BigInt::zero() {
            if &(&self.y + &other.y) % n == BigInt::zero() {
                return Ok(Self::inf());
            }
            let lambda = &(&self.x * &self.x * 3) + &curve.a;
            let lambda: BigInt = &lambda % n;
            let den: BigInt = &(&self.y * 2) % n;
            let den2 = &(&den * &den) % n;
            let den3 = &(&den2 * &den) % n;
            let x3 = &lambda * &lambda - (&self.x * 2) * &den2;
            let y3 = &lambda * &(&(&self.x * &den2) - &x3);
            let y3 = y3 - &self.y * &den3;
            return Self {
                x: zmod::<BigInt>(&(x3 * den), n),
                y: zmod::<BigInt>(&y3, n),
                z: den3,
            }
            .simplify(curve);
        }
        let lambda = zmod::<BigInt>(&(&self.y - &other.y), n);
        let xdif2 = &(&xdif * &xdif) % n;
        let xdif3 = &(&xdif2 * &xdif) % n;
        let x3 = &lambda * &lambda - &(&self.x + &other.x) * &xdif2;
        let y3 = &lambda * &(&(&self.x * &xdif2) - &x3);
        let y3 = y3 - &self.y * &xdif3;
        Self {
            x: zmod::<BigInt>(&(x3 * xdif), n),
            y: zmod::<BigInt>(&y3, n),
            z: xdif3,
        }
        .simplify(curve)
    }
    fn mul(&self, mut e: BigInt, curve: &Ell) -> Result<Self, BigInt> {
        let mut sum = Self::inf();
        let mut cur = self.clone();
        while e > BigInt::zero() {
            if &e % 2 == BigInt::one() {
                sum = sum.add(&cur, curve)?;
            }
            e /= 2;
            if e == BigInt::zero() {
                break;
            }
            cur = cur.add(&cur, curve)?;
        }
        Ok(sum)
    }
    fn inf() -> Self {
        Self {
            x: BigInt::zero(),
            y: BigInt::one(),
            z: BigInt::zero(),
        }
    }
    fn is_inf(&self) -> bool {
        self.z == BigInt::zero()
    }
    fn simplify(&self, curve: &Ell) -> Result<Self, BigInt> {
        if self.z == BigInt::zero() {
            return Ok(Self::inf());
        }
        let n = &curve.n;
        let invz = inv(&self.z, n)?;
        let x = &(&self.x * &invz) % n;
        let y = &(&self.y * &invz) % n;
        Ok(Self {
            x,
            y,
            z: BigInt::one(),
        })
    }
}

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
    let q = a / b;
    let r = a - &(b * &q);
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
        let (g, x, y) = extgcd(&-a, b);
        return (g, -x, y);
    }
    if b < &BigInt::zero() {
        let (g, x, y) = extgcd(a, &-b);
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
    debug_assert!(a % 2 == BigInt::one(), "{a}, {b}");
    if b % 2 == BigInt::zero() {
        let b1 = b >> 1;
        let (g, mut x, mut y) = extgcd_1(a, &b1);
        if &y % 2 == BigInt::zero() {
            y >>= 1;
        } else {
            x -= b1;
            y += a;
            debug_assert_eq!(&y % 2, BigInt::zero());
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
fn inv(a: &BigInt, mo: &BigInt) -> Result<BigInt, BigInt> {
    let (g, x, _y) = extgcd(a, mo);
    if g.abs() != BigInt::one() {
        return Err(g.abs());
    }
    // g = \pm 1. Now xg = 1 (mod mo) holds.
    Ok(zmod::<BigInt>(&(&x * &g), mo))
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
        let curve = Ell {
            a: 3.into(),
            n: 5.into(),
        };
        assert_eq!(
            Point {
                x: 2.into(),
                y: 1.into(),
                z: 1.into()
            }
            .add(
                &Point {
                    x: 1.into(),
                    y: 4.into(),
                    z: 1.into()
                },
                &curve
            )
            .unwrap(),
            Point {
                x: 1.into(),
                y: 1.into(),
                z: 1.into()
            }
        );
    }

    #[test]
    fn add_works_1() {
        // https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization#An_example
        let n = 455839i64;
        let curve = Ell {
            a: 5.into(),
            n: n.into(),
        };
        assert_eq!(
            Point {
                x: 1.into(),
                y: 1.into(),
                z: 1.into()
            }
            .add(
                &Point {
                    x: 1.into(),
                    y: 1.into(),
                    z: 1.into()
                },
                &curve
            )
            .unwrap(),
            Point {
                x: 14.into(),
                y: (n - 53).into(),
                z: 1.into()
            }
        );
    }

    #[test]
    fn mul_works_0() {
        let curve = Ell {
            a: 3.into(),
            n: 5.into(),
        };
        let a = Point {
            x: 1.into(),
            y: 1.into(),
            z: 1.into(),
        };
        assert!(a.mul(5.into(), &curve).unwrap().is_inf());
    }

    #[test]
    fn mul_works_1() {
        // https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization#An_example
        let n = 455839i64;
        let curve = Ell {
            a: 5.into(),
            n: n.into(),
        };
        let mut pt = Point {
            x: 1.into(),
            y: 1.into(),
            z: 1.into(),
        };
        let mut i = 1;
        let factor = loop {
            pt = match pt.mul(i.into(), &curve) {
                Ok(pt) => pt,
                Err(factor) => break factor,
            };
            i += 1;
        };
        assert_eq!(factor, 599.into());
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
        )
        .0;
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
        )
        .0;
        assert_eq!(a % factor, BigInt::zero());
    }

    // Too slow. This test takes about 16 sec in *release* build, about 1 min in debug build.
    #[test]
    #[ignore]
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
        )
        .0;
        eprintln!("factor={}", factor);
        assert_eq!(a % factor, BigInt::zero());
    }

    #[test]
    fn factorize_works_0() {
        let n = BigInt::from(1_000_000_007u128 * 1_000_000_007u128);
        let factors = factorize(&n);
        assert_eq!(factors.len(), 1);
    }
}
