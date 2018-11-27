extern crate num;

use std::ops::{Add, AddAssign};
use std::fmt::{Debug, Display};
use num::{BigInt, BigRational, One, Zero, pow};

// The leading coefficient (a[a.len() - 1]) must not be 0.
#[derive(Clone, PartialEq, Eq)]
pub struct Polynomial<R> {
    pub dat: Vec<R>
}

impl<R> Polynomial<R> {
    // If 0, returns usize::max_value().
    pub fn deg(&self) -> usize {
        if self.dat.len() == 0 {
            usize::max_value()
        } else {
            self.dat.len() - 1
        }
    }
    fn is_zero_primitive(&self) -> bool {
        self.dat.len() == 0
    }
}
impl<R: Zero> Polynomial<R> {
    pub fn from_raw(mut raw: Vec<R>) -> Self {
        let mut ma = 0;
        for i in (0 .. raw.len()).rev() {
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

impl<R: One + PartialEq> Polynomial<R> {
    fn is_monic(&self) -> bool {
        if self.is_zero_primitive() { return false; }
        let deg = self.deg();
        self.dat[deg].is_one()
    }
}

impl<'a, R: AddAssign + Clone + Zero> Add for &'a Polynomial<R> {
    type Output = Polynomial<R>;
    fn add(self, other: Self) -> Polynomial<R> {
        if self.dat.len() == 0 {
            return other.clone();
        }
        if other.dat.len() == 0 {
            return self.clone();
        }
        let self_deg = self.deg();
        let other_deg = other.deg();
        let ret_deg = std::cmp::max(self_deg, other_deg);
        let mut tmp = vec![R::zero(); ret_deg];
        for i in 0 .. ret_deg + 1 {
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
    type Output = Polynomial<R>;
    fn add(self, other: Self) -> Self {
        &self + &other
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

impl<R: Display + std::cmp::PartialEq + Zero> Debug for Polynomial<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        (self as &std::fmt::Display).fmt(f)
    }
}
impl<R: Display + std::cmp::PartialEq + Zero> Display for Polynomial<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_zero_primitive() {
            return Debug::fmt(&0, f);
        }
        let d = self.deg();
        let mut term_appear = false;
        for i in (0 .. d + 1).rev() {
            if self.dat[i].is_zero() { continue; }
            if term_appear {
                write!(f, " + ")?;
            }
            write!(f, "({})", self.dat[i])?;
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
pub fn div_rem_bigint(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>)
                      -> (Polynomial<BigInt>, Polynomial<BigInt>) {
    assert!(b.is_monic());
    // Falls back to pseudo_div_rem_bigint.
    pseudo_div_rem_bigint(a, b)
}

/// Performs pseudo-division.
pub fn pseudo_div_rem_bigint(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>)
                             -> (Polynomial<BigInt>, Polynomial<BigInt>) {
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
    for i in 0 .. a_deg + 1 {
        tmp[i] *= &factor;
    }

    // Naive division
    for i in (0 .. diff + 1).rev() {
        let mut coef = &tmp[i + b_deg] / &lcb;
        for j in 0 .. b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}

pub fn div_rem_bigrational(a: &Polynomial<BigRational>, b: &Polynomial<BigRational>)
               -> (Polynomial<BigRational>, Polynomial<BigRational>) {
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
    for i in (0 .. a_deg - b_deg + 1).rev() {
        let mut coef = &tmp[i + b_deg] / lc;
        for j in 0 .. b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}

#[cfg(test)]
mod tests {
    use num::BigInt;
    use super::{Polynomial, div_rem_bigint, div_rem_bigrational, pseudo_div_rem_bigint};
    #[test]
    fn test_div_rem_bigint() {
        // x^4 + x^2 + 1
        let p1 = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into(), 0.into(), 1.into()]);
        // x^2 + 2x + 3
        let p2 = Polynomial::from_raw(vec![3.into(), 2.into(), 1.into()]);
        let (q, r) = div_rem_bigint(&p1, &p2);
        assert_eq!(q, Polynomial::from_raw(vec![2.into(), (-2).into(), 1.into()]));
        assert_eq!(r, Polynomial::from_raw(vec![(-5).into(), 2.into()]));
    }
    #[test]
    fn test_div_rem_bigrational() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Vec<BigInt> =
            vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()];
        let p1 = Polynomial::from_raw(
            p1.into_iter().map(|x| x.into()).collect());
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Vec<BigInt> =
            vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()];
        let p2 = Polynomial::from_raw(
            p2.into_iter().map(|x| x.into()).collect());
        let (q, r) = div_rem_bigrational(&p1, &p2);
        let q_ex: Vec<(BigInt, BigInt)> =
            vec![(33.into(), 49.into()), (9.into(), 7.into())];
        let q_ex = Polynomial::from_raw(
            q_ex.into_iter().map(|x| x.into()).collect());
        let r_ex: Vec<(BigInt, BigInt)> =
            vec![(47.into(), 49.into()),
                 ((-576).into(), 49.into()),
                 ((-478).into(), 49.into()),
                 ((-411).into(), 49.into())];
        let r_ex = Polynomial::from_raw(
            r_ex.into_iter().map(|x| x.into()).collect());
        assert_eq!(q, q_ex);
        assert_eq!(r, r_ex);
    }
    #[test]
    fn test_pseudo_div_rem_bigint() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Vec<BigInt> =
            vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()];
        let p1 = Polynomial::from_raw(p1);
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Vec<BigInt> =
            vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()];
        let p2 = Polynomial::from_raw(p2);
        let (q, r) = pseudo_div_rem_bigint(&p1, &p2);
        let q_ex: Vec<BigInt> =
            vec![33.into(), 63.into()];
        let q_ex = Polynomial::from_raw(q_ex);
        let r_ex: Vec<BigInt> =
            vec![47.into(),
                 (-576).into(),
                 (-478).into(),
                 (-411).into()];
        let r_ex = Polynomial::from_raw(r_ex);
        assert_eq!(q, q_ex);
        assert_eq!(r, r_ex);
    }
}
