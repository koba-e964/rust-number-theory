use crate::numerical_roots::find_roots;
use crate::polynomial::Polynomial;
use num::{BigInt, Complex, Zero};

pub fn discriminant(f: &Polynomial<BigInt>) -> BigInt {
    assert!(!f.is_zero());
    let mut res = crate::resultant::resultant(f, &f.differential());
    let m = f.deg();
    if m % 4 == 2 || m % 4 == 3 {
        res = -res;
    }
    res / &f.dat[m]
}

pub fn discriminant_inaccurate(f: &Polynomial<Complex<f64>>) -> Complex<f64> {
    let roots = find_roots(f.clone());
    let n = roots.len();
    let mut prod = Complex::new(1.0, 0.0);
    for i in 0..n {
        for j in 0..i {
            prod *= (roots[i] - roots[j]).powi(2);
        }
    }
    for _ in 0..2 * n - 2 {
        prod *= f.coef_at(f.deg());
    }
    prod
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::Polynomial;
    use num::BigInt;
    #[test]
    fn test_discriminant_linear() {
        // 1771x + 24
        let p: Polynomial<BigInt> = Polynomial::from_raw(vec![24.into(), 1771.into()]);
        assert_eq!(discriminant(&p), 1.into());
    }
    #[test]
    fn test_discriminant_quadratic() {
        // 31x^2 + 1771x + 24
        let p: Polynomial<BigInt> = Polynomial::from_raw(vec![24.into(), 1771.into(), 31.into()]);
        assert_eq!(discriminant(&p), (1771 * 1771 - 4 * 24 * 31).into());
    }
    #[test]
    fn test_discriminant_cubic() {
        // x^3 + 9x + 1
        let p: Polynomial<BigInt> =
            Polynomial::from_raw(vec![1.into(), 9.into(), 0.into(), 1.into()]);
        assert_eq!(discriminant(&p), (-2943).into());
    }
    #[test]
    fn test_discriminant_cubic_2() {
        // 2x^3 + x^2 - 2x + 3
        let p: Polynomial<BigInt> =
            Polynomial::from_raw(vec![3.into(), (-2).into(), 1.into(), 2.into()]);
        assert_eq!(discriminant(&p), (-1132).into());
    }

    #[test]
    fn test_discriminant_cubic_inaccurate_2() {
        // 2x^3 + x^2 - 2x + 3
        let p: Polynomial<Complex<f64>> =
            Polynomial::from_raw(vec![3.0.into(), (-2.0).into(), 1.0.into(), 2.0.into()]);
        let d = discriminant_inaccurate(&p);
        let diff = d + 1132.0;
        assert!(diff.norm() <= 1.0e-6);
    }
}
