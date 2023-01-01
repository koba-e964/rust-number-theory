extern crate num;

use crate::{
    poly_mod::{poly_div, poly_mul},
    polynomial::{div_rem_bigrational, pseudo_div_rem_bigint, Polynomial},
};
use num::{pow, BigInt, BigRational, Integer, One, Zero};

/// Using naive arithmetic
pub fn resultant_rational(a: &Polynomial<BigRational>, b: &Polynomial<BigRational>) -> BigRational {
    if a.is_zero() || b.is_zero() {
        return BigRational::zero();
    }
    let a_deg = a.deg();
    let b_deg = b.deg();
    // b is a constant
    if b_deg == 0 {
        return pow(b.dat[0].clone(), a_deg);
    }
    let (_, r) = div_rem_bigrational(a, b);
    let r_deg = r.deg();
    if r_deg == usize::max_value() {
        return BigRational::zero();
    }
    let mut sub = resultant_rational(b, &r);
    sub *= pow(b.dat[b_deg].clone(), a_deg - r_deg);
    if a_deg % 2 == 1 && b_deg % 2 == 1 {
        sub = -sub;
    }
    sub
}

#[allow(clippy::many_single_char_names)]
fn resultant_smart_gcd(f: &Polynomial<BigInt>, g: &Polynomial<BigInt>) -> Polynomial<BigInt> {
    if f.is_zero() {
        return g.clone();
    }
    let contf = f.content();
    let contg = g.content();
    let d = contf.gcd(&contg);
    let mut f = poly_div(f, &contf);
    let mut g = poly_div(g, &contg);
    let mut a = BigInt::one();
    let mut b = BigInt::one();
    loop {
        if g.is_zero() {
            break;
        }
        let f_deg = f.deg();
        let g_deg = g.deg();
        if g_deg == 0 {
            f = Polynomial::<BigInt>::from_mono(BigInt::one());
            break;
        }
        if f_deg < g_deg {
            std::mem::swap(&mut f, &mut g);
            continue;
        }
        let delta = f_deg - g_deg;
        let (_, h) = pseudo_div_rem_bigint(&f, &g);
        f = g;
        g = h;
        // Divide g by a * b^delta
        let factor = &a * pow(b.clone(), delta);
        for i in 0..g.dat.len() {
            g.dat[i] /= &factor;
        }
        a = f.dat[f.deg()].clone();
        b = pow(a.clone(), delta) * &b / pow(b, delta);
    }
    let contf = f.content();
    poly_mul(&poly_div(&f, &contf), &d)
}

#[allow(clippy::many_single_char_names)]
fn resultant_smart(f: &Polynomial<BigInt>, g: &Polynomial<BigInt>) -> BigInt {
    if f.is_zero() {
        return BigInt::zero();
    }
    let mut a = BigInt::one();
    let mut b = BigInt::one();
    let mut f = f.clone();
    let mut g = g.clone();
    let mut s = 1;
    loop {
        if g.is_zero() {
            return BigInt::zero();
        }
        let f_deg = f.deg();
        let g_deg = g.deg();
        if f_deg % 2 == 1 && g_deg % 2 == 1 {
            s = -s;
        }
        if g_deg == 0 {
            break;
        }
        if f_deg < g_deg {
            std::mem::swap(&mut f, &mut g);
            continue;
        }
        let delta = f_deg - g_deg;
        let (_, h) = pseudo_div_rem_bigint(&f, &g);
        f = g;
        g = h;
        // Divide g by a * b^delta
        let factor = &a * pow(b.clone(), delta);
        for i in 0..g.dat.len() {
            g.dat[i] /= &factor;
        }
        a = f.dat[f.deg()].clone();
        b = pow(a.clone(), delta) * &b / pow(b, delta);
    }
    debug_assert_eq!(g.deg(), 0);
    debug_assert!(f.deg() >= 1);
    let mut result = pow(g.dat.swap_remove(0), f.deg());
    result /= pow(b, f.deg() - 1);
    if s == -1 {
        result = -result;
    }
    result
}

pub fn resultant(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>) -> BigInt {
    resultant_smart(a, b)
}

pub fn resultant_gcd(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>) -> Polynomial<BigInt> {
    resultant_smart_gcd(a, b)
}

#[cfg(test)]
mod tests {
    use super::resultant;
    use crate::polynomial::Polynomial;
    use num::BigInt;
    #[test]
    fn test_resultant() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Polynomial<BigInt> = Polynomial::from_raw(vec![
            5.into(),
            0.into(),
            2.into(),
            0.into(),
            6.into(),
            9.into(),
        ]);
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Polynomial<BigInt> =
            Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
        assert_eq!(resultant(&p1, &p2), 335159672.into());
    }
    #[test]
    fn test_resultant_2() {
        // 2x^2 + 5x + 2
        let p: Polynomial<BigInt> = Polynomial::from_raw(vec![2.into(), 5.into(), 2.into()]);
        // x^2 + 2
        let q: Polynomial<BigInt> = Polynomial::from_raw(vec![2.into(), 0.into(), 1.into()]);
        assert_eq!(resultant(&p, &q), 54.into());
    }
    // delta = 2
    #[test]
    fn test_resultant_3() {
        // x^4 + x^2 + 2
        let p: Polynomial<BigInt> =
            Polynomial::from_raw(vec![2.into(), 0.into(), 1.into(), 0.into(), 1.into()]);
        // x^2 + 1
        let q: Polynomial<BigInt> = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into()]);
        assert_eq!(resultant(&p, &q), 4.into());
    }
    // delta = 3
    #[test]
    fn test_resultant_4() {
        // x^6 + x^3 + 2
        let p: Polynomial<BigInt> = Polynomial::from_raw(vec![
            2.into(),
            0.into(),
            0.into(),
            1.into(),
            0.into(),
            0.into(),
            1.into(),
        ]);
        // x^3 + 1
        let q: Polynomial<BigInt> =
            Polynomial::from_raw(vec![1.into(), 0.into(), 0.into(), 1.into()]);
        assert_eq!(resultant(&p, &q), 8.into());
    }
}
