use num::Complex;
use rand::Rng;

use crate::polynomial::Polynomial;

const EPS: f64 = 1.0e-18;

pub fn find_roots_reim(mut poly: Polynomial<Complex<f64>>) -> (Vec<f64>, Vec<Complex<f64>>) {
    let mut re = vec![];
    let mut im = vec![];
    // Randomly pick an initial starting point
    let mut rng = rand::thread_rng();
    let mut trial = 3;
    while poly.deg() > 0 && trial > 0 {
        let der = poly.differential();
        let r = rng.gen_range(0.0..2.0);
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let x = Complex::from_polar(r, theta);
        if let Some(x) = find_once(&poly, &der, x) {
            if x.im.abs() <= EPS {
                re.push(x.re);
            } else {
                im.push(x);
                poly = divide_by_x_a(&poly, x.conj());
            }
            poly = divide_by_x_a(&poly, x);
        } else {
            trial -= 1;
        }
    }
    assert_eq!(poly.deg(), 0);
    (re, im)
}

pub fn find_roots(mut poly: Polynomial<Complex<f64>>) -> Vec<Complex<f64>> {
    let mut roots = vec![];
    // Randomly pick an initial starting point
    let mut rng = rand::thread_rng();
    let mut trial = 3;
    while poly.deg() > 0 && trial > 0 {
        let der = poly.differential();
        let r = rng.gen_range(0.0..2.0);
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let x = Complex::from_polar(r, theta);
        if let Some(x) = find_once(&poly, &der, x) {
            roots.push(x);
            poly = divide_by_x_a(&poly, x);
        } else {
            trial -= 1;
        }
    }
    assert_eq!(poly.deg(), 0);
    roots
}

fn find_once(
    poly: &Polynomial<Complex<f64>>,
    der: &Polynomial<Complex<f64>>,
    mut x: Complex<f64>,
) -> Option<Complex<f64>> {
    let mut iter = 0;
    loop {
        let val = poly.of(&x);
        let diff = der.of(&x);
        x -= val / diff;
        if val.norm_sqr() <= EPS {
            break;
        }
        if iter >= 100 {
            return None;
        }
        iter += 1;
    }
    Some(x)
}

fn divide_by_x_a(p: &Polynomial<Complex<f64>>, a: Complex<f64>) -> Polynomial<Complex<f64>> {
    let deg = p.deg();
    let mut coefs = vec![Complex::default(); deg];
    let mut carry = Complex::default();
    for i in (0..deg).rev() {
        carry *= a;
        carry += p.coef_at(i + 1);
        coefs[i] = carry;
    }
    Polynomial::from_raw(coefs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_roots_works() {
        let p = vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        ];
        let p = Polynomial::from_raw(p);
        let roots = find_roots(p);
        assert_eq!(roots.len(), 2);
    }
}
