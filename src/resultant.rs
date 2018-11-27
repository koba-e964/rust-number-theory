extern crate num;

use num::{BigInt, BigRational, Zero, pow};
use polynomial::{Polynomial, div_rem_bigrational};

// TODO using naive arithmetic
pub fn resultant_rational(a: &Polynomial<BigRational>, b: &Polynomial<BigRational>)
             -> BigRational {
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
    if r_deg == usize::max_value() { return BigRational::zero(); }
    let mut sub = resultant_rational(b, &r);
    sub *= pow(b.dat[b_deg].clone(), a_deg - r_deg);
    if a_deg % 2 == 1 && b_deg % 2 == 1 {
        sub = -sub;
    }
    sub
}
// Very naive way to calculate resultant(a, b).
pub fn resultant(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>)
                 -> BigInt {
    let ar = Polynomial::from_raw(
        a.dat.iter().map(|x| x.clone().into()).collect());
    let br = Polynomial::from_raw(
        b.dat.iter().map(|x| x.clone().into()).collect());
    let result = resultant_rational(&ar, &br);
    assert!(result.is_integer());
    result.to_integer()
}


#[cfg(test)]
mod tests {
    use num::BigInt;
    use polynomial::Polynomial;
    use super::resultant;
    #[test]
    fn test_resultant() {
        // 9x^5 + 6x^4 + 2x^2 + 5
        let p1: Polynomial<BigInt> =
            Polynomial::from_raw(vec![5.into(), 0.into(), 2.into(), 0.into(), 6.into(), 9.into()]);
        // 7x^4 + x^3 + 6x^2 + 6x + 6
        let p2: Polynomial<BigInt> =
            Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
        assert_eq!(resultant(&p1, &p2), 335159672.into());
    }
}
