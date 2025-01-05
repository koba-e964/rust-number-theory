use crate::polynomial::Polynomial;
use num::{traits::Pow, BigInt, BigRational, One, Zero};
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Algebraic {
    pub min_poly: Polynomial<BigInt>,
    pub expr: Polynomial<BigRational>,
}

impl Algebraic {
    /// minimal_poly should be irreducible in Z[x].
    pub fn new(minimal_poly: Polynomial<BigInt>) -> Self {
        Algebraic {
            min_poly: minimal_poly,
            expr: Polynomial::from_raw(vec![
                BigRational::from_integer(0.into()),
                BigRational::from_integer(1.into()),
            ]),
        }
    }
    pub fn with_expr(minimal_poly: Polynomial<BigInt>, expr: Polynomial<BigRational>) -> Self {
        debug_assert!(expr.deg() < minimal_poly.deg());
        Algebraic {
            min_poly: minimal_poly,
            expr,
        }
    }
    pub fn new_const(minimal_poly: Polynomial<BigInt>, x: BigRational) -> Self {
        Algebraic {
            min_poly: minimal_poly,
            expr: Polynomial::from_raw(vec![x]),
        }
    }
    pub fn from_int(minimal_poly: Polynomial<BigInt>, x: impl Into<BigInt>) -> Self {
        let u: BigInt = x.into();
        Algebraic {
            min_poly: minimal_poly,
            expr: Polynomial::from_raw(vec![u.into()]),
        }
    }
    /// Returns its minimal polynomial's degree.
    pub fn deg(&self) -> usize {
        self.min_poly.deg()
    }

    pub fn as_coefs(&self) -> Vec<BigRational> {
        let mut expr = self.expr.clone().dat;
        let deg = self.deg();
        expr.extend_from_slice(&vec![BigRational::from_integer(0.into()); deg - expr.len()]);
        expr
    }
}

// Operations on Algebraic assume that all numbers' min_poly are the same.

impl Add for Algebraic {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Algebraic {
            min_poly: self.min_poly,
            expr: self.expr + other.expr,
        }
    }
}
impl Add for &Algebraic {
    type Output = Algebraic;
    fn add(self, other: Self) -> Algebraic {
        Algebraic {
            min_poly: self.min_poly.clone(),
            expr: &self.expr + &other.expr,
        }
    }
}

impl Sub for Algebraic {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Algebraic {
            min_poly: self.min_poly,
            expr: self.expr - other.expr,
        }
    }
}
impl Sub for &Algebraic {
    type Output = Algebraic;
    fn sub(self, other: Self) -> Algebraic {
        Algebraic {
            min_poly: self.min_poly.clone(),
            expr: &self.expr - &other.expr,
        }
    }
}

/// Computes (a * b) mod c.
fn mul_with_mod<T: Into<BigRational> + Clone>(
    a: &Polynomial<BigRational>,
    b: &Polynomial<BigRational>,
    c: &Polynomial<T>,
) -> Polynomial<BigRational> {
    if a.is_zero() || b.is_zero() {
        return Polynomial::zero();
    }
    let n = c.deg();
    let a_deg = a.deg();
    let b_deg = b.deg();
    assert!(a_deg < n);
    assert!(b_deg < n);
    let mut result = vec![BigRational::from_integer(0.into()); n];
    // Stores b * x^i
    let mut cur = vec![BigRational::from_integer(0.into()); n + 1];
    cur[..b_deg + 1].clone_from_slice(&b.dat[..b_deg + 1]);
    let lc = c.dat[n].clone().into();
    for i in 0..a_deg + 1 {
        for j in 0..n {
            result[j] += &a.dat[i] * &cur[j];
        }
        if i < a_deg {
            // cur = cur * x
            for j in (0..n).rev() {
                cur.swap(j, j + 1);
            }
            // cur = cur % c
            let coef = &cur[n] / &lc;
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                cur[j] -= &coef * &c.dat[j].clone().into();
            }
            cur[n] = BigRational::from_integer(0.into());
        }
    }
    Polynomial::from_raw(result)
}

impl Mul for Algebraic {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Algebraic {
            expr: mul_with_mod(&self.expr, &other.expr, &self.min_poly),
            min_poly: self.min_poly,
        }
    }
}
impl Mul for &Algebraic {
    type Output = Algebraic;
    fn mul(self, other: Self) -> Algebraic {
        Algebraic {
            min_poly: self.min_poly.clone(),
            expr: mul_with_mod(&self.expr, &other.expr, &self.min_poly),
        }
    }
}

impl Pow<u64> for &Algebraic {
    type Output = Algebraic;
    fn pow(self, mut e: u64) -> Self::Output {
        let mut cur = self.clone();
        let mut prod = Algebraic::from_int(self.min_poly.clone(), 1);
        while e > 0 {
            if e % 2 == 1 {
                prod = &prod * &cur;
            }
            cur = &cur * &cur;
            e /= 2;
        }
        prod
    }
}

impl Pow<BigInt> for &Algebraic {
    type Output = Algebraic;
    fn pow(self, mut e: BigInt) -> Self::Output {
        let mut cur = self.clone();
        let mut prod = Algebraic::from_int(self.min_poly.clone(), 1);
        while e > BigInt::zero() {
            if &e % 2 == BigInt::one() {
                prod = &prod * &cur;
            }
            cur = &cur * &cur;
            e /= 2;
        }
        prod
    }
}

#[cfg(test)]
mod tests {
    use super::Algebraic;
    use crate::polynomial::Polynomial;
    #[test]
    fn test_alg_mul() {
        // Let theta be an algebraic number whose minimal polynomial is x^3 + x + 1.
        // Let eta = theta^2.
        // eta^3 + 2 * eta^2 + eta - 1 = 0
        let f = Polynomial::from_raw(vec![1.into(), 1.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(f.clone());
        let eta = &theta * &theta;
        let result = &eta * &(&eta * &eta) + &eta * &eta * Algebraic::from_int(f.clone(), 2) + eta
            - Algebraic::from_int(f.clone(), 1);
        assert_eq!(result, Algebraic::from_int(f, 0));
    }
}
