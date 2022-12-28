extern crate num;

use num::Signed;
use num::{pow, traits::NumAssign, BigInt, BigRational, Complex, Integer, One, Zero};
use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// The leading coefficient (a[a.len() - 1]) must not be 0.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Polynomial<R> {
    pub dat: Vec<R>,
}

impl<R> Polynomial<R> {
    // If 0, returns usize::max_value().
    pub fn deg(&self) -> usize {
        if self.dat.is_empty() {
            usize::max_value()
        } else {
            self.dat.len() - 1
        }
    }
    fn is_zero_primitive(&self) -> bool {
        self.dat.is_empty()
    }
}
impl<R: Zero> Polynomial<R> {
    pub fn from_raw(mut raw: Vec<R>) -> Self {
        let mut ma = 0;
        for i in (0..raw.len()).rev() {
            if !raw[i].is_zero() {
                ma = i + 1;
                break;
            }
        }
        raw.drain(ma..);
        Polynomial { dat: raw }
    }
    pub fn from_mono(v: impl Into<R>) -> Self {
        Polynomial::from_raw(vec![v.into()])
    }
}

impl<R: Zero + Clone> Polynomial<R> {
    /// Get the coefficient of x^index in the polynomial.
    pub fn coef_at(&self, index: usize) -> R {
        if index < self.dat.len() {
            return self.dat[index].clone();
        }
        R::zero()
    }
}

impl<Int: Clone + NumAssign + Integer + From<i32>> Polynomial<Int> {
    pub fn differential(&self) -> Polynomial<Int> {
        if self.is_zero_primitive() {
            return self.clone();
        }
        let deg = self.deg();
        let mut tmp = vec![0.into(); deg];
        #[allow(clippy::needless_range_loop)]
        for i in 0..deg {
            tmp[i] = self.dat[i + 1].clone();
            tmp[i] *= Int::from(i as i32 + 1);
        }
        Polynomial::from_raw(tmp)
    }
}

impl Polynomial<Complex<f64>> {
    pub fn differential_complex(&self) -> Polynomial<Complex<f64>> {
        if self.is_zero_primitive() {
            return self.clone();
        }
        let deg = self.deg();
        let mut tmp = vec![Complex::default(); deg];
        #[allow(clippy::needless_range_loop)]
        for i in 0..deg {
            tmp[i] = self.dat[i + 1] * (i + 1) as f64;
        }
        Polynomial::from_raw(tmp)
    }
}

impl<R: One + PartialEq> Polynomial<R> {
    fn is_monic(&self) -> bool {
        if self.is_zero_primitive() {
            return false;
        }
        let deg = self.deg();
        self.dat[deg].is_one()
    }
}

impl<R: AddAssign + Zero + MulAssign + Clone> Polynomial<R> {
    pub fn of(&self, x: &R) -> R {
        let mut sum = R::zero();
        let deg = self.deg();
        for i in (0..deg + 1).rev() {
            sum *= x.clone();
            sum += self.coef_at(i);
        }
        sum
    }
}

impl<'a, R: AddAssign + Clone + Zero> Add for &'a Polynomial<R> {
    type Output = Polynomial<R>;
    fn add(self, other: Self) -> Polynomial<R> {
        if self.dat.is_empty() {
            return other.clone();
        }
        if other.dat.is_empty() {
            return self.clone();
        }
        let self_deg = self.deg();
        let other_deg = other.deg();
        let ret_deg = std::cmp::max(self_deg, other_deg);
        let mut tmp = vec![R::zero(); ret_deg + 1];
        #[allow(clippy::needless_range_loop)]
        for i in 0..ret_deg + 1 {
            if i <= self_deg {
                tmp[i] += self.dat[i].clone();
            }
            if i <= other_deg {
                tmp[i] += other.dat[i].clone();
            }
        }
        Polynomial::from_raw(tmp)
    }
}
impl<R: AddAssign + Clone + Zero> Add for Polynomial<R> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if self.dat.is_empty() {
            return other;
        }
        if other.dat.is_empty() {
            return self;
        }
        &self + &other
    }
}

impl<R: Neg<Output = R>> Neg for Polynomial<R> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut dat = self.dat;
        for coef in dat.iter_mut() {
            unsafe {
                let mut tmp = std::ptr::read(coef as *mut R);
                tmp = -tmp;
                std::ptr::write(coef as *mut _, tmp);
            }
        }
        Polynomial { dat }
    }
}
impl<'a, R: Neg<Output = R> + Clone> Neg for &'a Polynomial<R> {
    type Output = Polynomial<R>;
    fn neg(self) -> Polynomial<R> {
        -self.clone()
    }
}

impl<'a, R: AddAssign + SubAssign + Neg<Output = R> + Clone + Zero> Sub for &'a Polynomial<R> {
    type Output = Polynomial<R>;
    fn sub(self, other: Self) -> Polynomial<R> {
        if self.dat.is_empty() {
            return -other;
        }
        if other.dat.is_empty() {
            return self.clone();
        }
        let self_deg = self.deg();
        let other_deg = other.deg();
        let ret_deg = std::cmp::max(self_deg, other_deg);
        let mut tmp = vec![R::zero(); ret_deg + 1];
        #[allow(clippy::needless_range_loop)]
        for i in 0..ret_deg + 1 {
            if i <= self_deg {
                tmp[i] += self.dat[i].clone();
            }
            if i <= other_deg {
                tmp[i] -= other.dat[i].clone();
            }
        }
        Polynomial::from_raw(tmp)
    }
}
impl<R: AddAssign + SubAssign + Neg<Output = R> + Clone + Zero> Sub for Polynomial<R> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        if self.dat.is_empty() {
            return -other;
        }
        if other.dat.is_empty() {
            return self;
        }
        &self - &other
    }
}

impl<'a, R: AddAssign + Clone + Zero> Mul for &'a Polynomial<R>
where
    for<'b> &'b R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Self) -> Polynomial<R> {
        if self.is_zero_primitive() || other.is_zero_primitive() {
            return Polynomial::from_raw(Vec::new());
        }
        let a_deg = self.deg();
        let b_deg = other.deg();
        let mut result = vec![R::zero(); a_deg + b_deg + 1];
        for i in 0..a_deg + 1 {
            for j in 0..b_deg + 1 {
                result[i + j] += &self.dat[i] * &other.dat[j];
            }
        }
        Polynomial::from_raw(result)
    }
}
impl<R: AddAssign + Clone + Zero> Mul for Polynomial<R>
where
    for<'a> &'a R: Mul<Output = R>,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        &self * &other
    }
}

impl<R: AddAssign + Clone + Zero> Zero for Polynomial<R> {
    fn is_zero(&self) -> bool {
        self.is_zero_primitive()
    }
    fn zero() -> Self {
        Polynomial { dat: Vec::new() }
    }
}

impl<R: Display + std::cmp::PartialEq + Zero + One + Signed> Debug for Polynomial<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        (self as &dyn std::fmt::Display).fmt(f)
    }
}
impl<R: Display + std::cmp::PartialEq + Zero + One + Signed> Display for Polynomial<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_zero_primitive() {
            return Debug::fmt(&0, f);
        }
        let d = self.deg();
        let mut term_appear = false;
        for i in (0..d + 1).rev() {
            if self.dat[i].is_zero() {
                continue;
            }
            if term_appear {
                write!(f, " + ")?;
            }
            if self.dat[i].is_one() && i > 0 {
            } else {
                if self.dat[i].is_positive() {
                    write!(f, "{}", self.dat[i])?;
                } else {
                    write!(f, "({})", self.dat[i])?;
                }
            }
            if i >= 2 {
                write!(f, "X^{}", i)?;
            } else if i == 1 {
                write!(f, "X")?;
            }
            term_appear = true;
        }
        Ok(())
    }
}

// b should be monic, or it panics.
pub fn div_rem_bigint(
    a: &Polynomial<BigInt>,
    b: &Polynomial<BigInt>,
) -> (Polynomial<BigInt>, Polynomial<BigInt>) {
    assert!(b.is_monic());
    // Falls back to pseudo_div_rem_bigint.
    pseudo_div_rem_bigint(a, b)
}

/// Performs pseudo-division.
pub fn pseudo_div_rem_bigint(
    a: &Polynomial<BigInt>,
    b: &Polynomial<BigInt>,
) -> (Polynomial<BigInt>, Polynomial<BigInt>) {
    let a_deg = a.deg();
    let b_deg = b.deg();
    if a.is_zero() || b.is_zero() || a_deg < b_deg {
        return (Polynomial::from_mono(0), a.clone());
    }
    let lcb = b.dat[b_deg].clone();
    assert!(!b.is_zero());
    let mut tmp = a.dat.clone();
    let mut quo = vec![0.into(); a_deg - b_deg + 1];
    let diff = a_deg - b_deg;
    // Multiply a (tmp) by lcb^{diff + 1}
    let factor = pow(lcb.clone(), diff + 1);
    for coef in tmp.iter_mut() {
        *coef *= &factor;
    }

    // Naive division
    for i in (0..diff + 1).rev() {
        let coef = &tmp[i + b_deg] / &lcb;
        for j in 0..b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}

pub fn div_rem_bigrational(
    a: &Polynomial<BigRational>,
    b: &Polynomial<BigRational>,
) -> (Polynomial<BigRational>, Polynomial<BigRational>) {
    let a_deg = a.deg();
    let b_deg = b.deg();
    let zero: BigInt = 0.into();
    if a.is_zero() || b.is_zero() || a_deg < b_deg {
        return (Polynomial::from_mono(zero), a.clone());
    }
    assert!(!b.dat[b_deg].is_zero());
    let lc = &b.dat[b_deg];
    let mut tmp = a.dat.clone();
    let mut quo = vec![zero.into(); a_deg - b_deg + 1];
    // Naive division
    for i in (0..a_deg - b_deg + 1).rev() {
        let coef = &tmp[i + b_deg] / lc;
        for j in 0..b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}

#[cfg(test)]
mod tests {
    use super::{div_rem_bigint, div_rem_bigrational, pseudo_div_rem_bigint, Polynomial};
    use num::{BigInt, Zero};
    #[test]
    fn test_sub_zero() {
        let p1: Polynomial<BigInt> = Polynomial::zero();
        let p2: Polynomial<BigInt> = Polynomial::from_raw(vec![1.into(), 2.into()]);
        assert_eq!(
            p1 - p2,
            Polynomial::from_raw(vec![(-1).into(), (-2).into()])
        );
    }
    #[test]
    fn test_mul_bigint() {
        // x^2 + x + 1
        let p1: Polynomial<BigInt> = Polynomial::from_raw(vec![1.into(), 1.into(), 1.into()]);
        // x - 1
        let p2: Polynomial<BigInt> = Polynomial::from_raw(vec![(-1).into(), 1.into()]);
        assert_eq!(
            p1 * p2,
            Polynomial::from_raw(vec![(-1).into(), 0.into(), 0.into(), 1.into()])
        );
    }
    #[test]
    fn test_div_rem_bigint() {
        // x^4 + x^2 + 1
        let p1 = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into(), 0.into(), 1.into()]);
        // x^2 + 2x + 3
        let p2 = Polynomial::from_raw(vec![3.into(), 2.into(), 1.into()]);
        let (q, r) = div_rem_bigint(&p1, &p2);
        assert_eq!(
            q,
            Polynomial::from_raw(vec![2.into(), (-2).into(), 1.into()])
        );
        assert_eq!(r, Polynomial::from_raw(vec![(-5).into(), 2.into()]));
    }
    #[test]
    fn test_div_rem_bigrational() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Vec<BigInt> = vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()];
        let p1 = Polynomial::from_raw(p1.into_iter().map(|x| x.into()).collect());
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Vec<BigInt> = vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()];
        let p2 = Polynomial::from_raw(p2.into_iter().map(|x| x.into()).collect());
        let (q, r) = div_rem_bigrational(&p1, &p2);
        let q_ex: Vec<(BigInt, BigInt)> = vec![(33.into(), 49.into()), (9.into(), 7.into())];
        let q_ex = Polynomial::from_raw(q_ex.into_iter().map(|x| x.into()).collect());
        let r_ex: Vec<(BigInt, BigInt)> = vec![
            (47.into(), 49.into()),
            ((-576).into(), 49.into()),
            ((-478).into(), 49.into()),
            ((-411).into(), 49.into()),
        ];
        let r_ex = Polynomial::from_raw(r_ex.into_iter().map(|x| x.into()).collect());
        assert_eq!(q, q_ex);
        assert_eq!(r, r_ex);
    }
    #[test]
    fn test_pseudo_div_rem_bigint() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Vec<BigInt> = vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()];
        let p1 = Polynomial::from_raw(p1);
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Vec<BigInt> = vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()];
        let p2 = Polynomial::from_raw(p2);
        let (q, r) = pseudo_div_rem_bigint(&p1, &p2);
        let q_ex: Vec<BigInt> = vec![33.into(), 63.into()];
        let q_ex = Polynomial::from_raw(q_ex);
        let r_ex: Vec<BigInt> = vec![47.into(), (-576).into(), (-478).into(), (-411).into()];
        let r_ex = Polynomial::from_raw(r_ex);
        assert_eq!(q, q_ex);
        assert_eq!(r, r_ex);
    }
}
