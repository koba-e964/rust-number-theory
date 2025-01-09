use num::bigint::RandBigInt;
use num::{BigInt, One, Zero};
use std::collections::HashMap;

use crate::ecm::{select_b, ECMConfig, EcmStats};
use crate::inverse::{inv, zmod};
use crate::perfect_power::perfect_power;
use crate::prime;

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

/// Finds a factor.
#[allow(clippy::many_single_char_names)]
pub fn ecm(n: &BigInt, conf: ECMConfig) -> (BigInt, u64) {
    debug_assert!(!prime::is_prime(n));

    let mut rng = rand::thread_rng();

    let mut count = 0u64;
    let parallel_count = 10; // TODO: make this configurable

    loop {
        count += 1;
        if conf.verbose {
            eprintln!(
                "Trying curves {}, B1 = {}, B2 = {}",
                count, conf.b1, conf.b2
            );
        }
        // randomly pick curves
        let mut curves = vec![];
        for _ in 0..parallel_count {
            let a = rng.gen_bigint_range(&1.into(), n);
            let curve = Ell { a, n: n.clone() };
            curves.push(curve);
        }

        // randomly select points
        let mut points = vec![];
        for _ in 0..parallel_count {
            let x = rng.gen_bigint_range(&1.into(), n);
            let y = rng.gen_bigint_range(&1.into(), n);
            let pt = Point {
                x,
                y,
                z: BigInt::one(),
            };
            points.push(pt);
        }
        if let Err(fac) = ecm_oneshot_parallel(points, curves, conf.b1, conf.b2) {
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

fn ecm_oneshot_parallel(pts: Vec<Point>, curves: Vec<Ell>, b1: u64, b2: u64) -> Result<(), BigInt> {
    let k = pts.len();
    let mut joint = Vec::with_capacity(k);
    for i in 0..k {
        joint.push((pts[i].clone(), curves[i].clone()));
    }
    for mult in 1..b1 + 1 {
        let mut tmp = Point::many_muls(&joint, mult.into())?;
        for i in 0..k {
            joint[i].0 = std::mem::take(&mut tmp[i]);
        }
    }
    // Step 2: try all primes in range (b1, b2]
    for &init in &[b1.saturating_sub(1) / 6 * 6 + 1, (b1 + 1) / 6 * 6 - 1] {
        let mut cur_e = init;
        let p6 = {
            let mut tmp = Vec::with_capacity(k);
            for i in 0..k {
                tmp.push((joint[i].0.clone(), joint[i].0.clone(), joint[i].1.clone()));
            }
            let p2 = Point::many_adds(&tmp)?;
            for i in 0..k {
                tmp[i].0 = p2[i].clone();
                tmp[i].1 = p2[i].clone();
            }
            let p4 = Point::many_adds(&tmp)?;
            for i in 0..k {
                tmp[i].1 = p4[i].clone();
            }
            Point::many_adds(&tmp)?
        };
        let tmp = Point::many_muls(&joint, init.into())?;
        for i in 0..k {
            joint[i].0 = tmp[i].clone();
        }
        let mut tmp = Vec::with_capacity(k);
        for i in 0..k {
            tmp.push((joint[i].0.clone(), p6[i].clone(), joint[i].1.clone()));
        }
        while cur_e <= b2 {
            cur_e += 6;
            let result = Point::many_adds(&tmp)?;
            for i in 0..k {
                tmp[i].0 = result[i].clone();
            }
        }
    }
    Ok(())
}

/// Projective coordinates
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Point {
    x: BigInt,
    y: BigInt,
    z: BigInt,
}

// y^2 = x^3 + ax + b (mod n)
// not (2 | n)
// In the actual computation, b is irrelevant.
#[derive(Debug, Clone)]
struct Ell {
    a: BigInt,
    n: BigInt,
}

impl Point {
    fn many_adds(pts: &[(Self, Self, Ell)]) -> Result<Vec<Self>, BigInt> {
        let mut points = Vec::with_capacity(pts.len());
        for (p1, p2, curve) in pts {
            let xdif = &p1.x - &p2.x;
            let n = &curve.n;
            if xdif == BigInt::zero() {
                if &(&p1.y + &p2.y) % n == BigInt::zero() {
                    points.push(Self::inf());
                    continue;
                }
                let lambda = &(&p1.x * &p1.x * 3) + &curve.a;
                let lambda: BigInt = &lambda % n;
                let den: BigInt = &(&p1.y * 2) % n;
                let den2 = &(&den * &den) % n;
                let den3 = &(&den2 * &den) % n;
                let x3 = &lambda * &lambda - (&p1.x * 2) * &den2;
                let y3 = &lambda * &(&(&p1.x * &den2) - &x3);
                let y3 = y3 - &p1.y * &den3;
                points.push(Self {
                    x: zmod::<BigInt>(&(x3 * den), n),
                    y: zmod::<BigInt>(&y3, n),
                    z: den3,
                });
            } else {
                let lambda = zmod::<BigInt>(&(&p1.y - &p2.y), n);
                let xdif2 = &(&xdif * &xdif) % n;
                let xdif3 = &(&xdif2 * &xdif) % n;
                let x3 = &lambda * &lambda - &(&p1.x + &p2.x) * &xdif2;
                let y3 = &lambda * &(&(&p1.x * &xdif2) - &x3);
                let y3 = y3 - &p1.y * &xdif3;
                points.push(Self {
                    x: zmod::<BigInt>(&(x3 * xdif), n),
                    y: zmod::<BigInt>(&y3, n),
                    z: xdif3,
                });
            }
        }
        Self::many_simplify(&points, &pts[0].2.n)
    }
    fn many_muls(pts: &[(Self, Ell)], mut e: BigInt) -> Result<Vec<Self>, BigInt> {
        let k = pts.len();
        let mut sum = vec![Self::inf(); k];
        let mut cur = pts.to_vec();
        while e > BigInt::zero() {
            if &e % 2 == BigInt::one() {
                let mut dat = vec![];
                for i in 0..k {
                    dat.push((sum[i].clone(), cur[i].0.clone(), cur[i].1.clone()));
                }
                sum = Self::many_adds(&dat)?;
            }
            e /= 2;
            if e == BigInt::zero() {
                break;
            }
            let mut dat = vec![];
            for i in 0..k {
                dat.push((cur[i].0.clone(), cur[i].0.clone(), cur[i].1.clone()));
            }
            let tmp = Self::many_adds(&dat)?;
            for i in 0..k {
                cur[i].0 = tmp[i].clone();
            }
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
    fn many_simplify(pts: &[Self], n: &BigInt) -> Result<Vec<Self>, BigInt> {
        let k = pts.len();
        let mut zarr = vec![];
        let mut zprod = BigInt::one();
        for i in 0..k {
            zarr.push(pts[i].z.clone());
            if zarr[i] != BigInt::zero() {
                zprod = &zprod * &zarr[i];
            }
        }
        let mut zacc_l = vec![BigInt::one(); k + 1];
        let mut zacc_r = vec![BigInt::one(); k + 1];
        for i in 0..k {
            if zarr[i] != BigInt::zero() {
                zacc_l[i + 1] = &(&zacc_l[i] * &zarr[i]) % n;
            } else {
                zacc_l[i + 1] = zacc_l[i].clone();
            }
        }
        for i in (0..k).rev() {
            if zarr[i] != BigInt::zero() {
                zacc_r[i] = &(&zacc_r[i + 1] * &zarr[i]) % n;
            } else {
                zacc_r[i] = zacc_r[i + 1].clone();
            }
        }
        let invzprod = inv(&zprod, n)?;
        let mut result = Vec::with_capacity(k);
        for i in 0..k {
            if zarr[i] == BigInt::zero() {
                result.push(Self::inf());
                continue;
            }
            let invz = &zacc_l[i] * &zacc_r[i + 1] * &invzprod % n;
            let x = &(&pts[i].x * &invz) % n;
            let y = &(&pts[i].y * &invz) % n;
            result.push(Self {
                x,
                y,
                z: BigInt::one(),
            })
        }
        Ok(result)
    }
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
