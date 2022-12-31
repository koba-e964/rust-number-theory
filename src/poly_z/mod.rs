use num::{BigInt, One, Signed, Zero};
use number_theory_elementary::Primes;

use crate::{
    poly_mod::{self, lift_factorization},
    polynomial::{pseudo_div_rem_bigint, Polynomial},
};

/// Factorizes a polynomial.
///
/// `a` must be squarefree and monic. (TODO: remove these constraints)
pub fn factorize(a: &Polynomial<BigInt>) -> Vec<(Polynomial<BigInt>, usize)> {
    // TODO: implement 1 Reduce to squarefree and primitive
    // a must be monic

    // Theorem 3.5.1 in [Cohen]
    let n = a.deg();
    assert_eq!(a.coef_at(n), BigInt::one(), "lc(a) == 1 must hold");
    let mut sum = a.coef_at(n).abs();
    for i in 0..n + 1 {
        sum += a.coef_at(i).abs();
    }
    let two = BigInt::one() + BigInt::one();
    let mut bound = sum;
    for _ in 0..n - 1 {
        bound = &bound * &two;
    }
    bound = &(&bound * &two) * &a.coef_at(n);

    let mut p = BigInt::zero();
    let mut pusize = 0;
    for now in Primes::new() {
        // check if (a, a') = 1 in F_p[X]
        let nowint: BigInt = (now as i32).into();
        let a = poly_mod::poly_mod(a, &nowint);
        let a_p = poly_mod::differential(&a, &nowint);
        let gcd = poly_mod::poly_gcd::<BigInt>(&a, &a_p, &nowint);
        if gcd.deg() == 0 {
            p = nowint;
            pusize = now;
            break;
        }
    }
    let mut e = 0;
    let mut pe = BigInt::one();
    while pe <= bound {
        pe = &pe * &p;
        e += 1;
    }
    let pe2 = &pe / &BigInt::from(2); // floor(p^e/2)
    let factors = poly_mod::factorize_mod_p::<BigInt>(a, &p, pusize);
    assert!(factors.iter().all(|&(_, e)| e == 1));
    let factors: Vec<Polynomial<BigInt>> = factors.into_iter().map(|(poly, _)| poly).collect();
    let mut lifted = lift_factorization::<BigInt>(&p, e, a, &factors);
    // 5. Try combination
    let mut d = 1;
    let mut a = a.clone();
    let mut result = vec![];
    'outer: while 2 * d <= lifted.len() {
        assert!(lifted.len() <= 25);
        for bits in 0usize..1 << lifted.len() {
            if bits.count_ones() as usize != d {
                continue;
            }
            let mut prod: Polynomial<BigInt> = Polynomial::from_mono(BigInt::one());
            for i in 0..lifted.len() {
                if (bits & 1 << i) != 0 {
                    prod = poly_mod::poly_mod(&(&prod * &lifted[i]), &pe);
                }
            }
            // modify prod so that all coefficients are in [-p^e/2, p^e/2)
            let bias = Polynomial::from_raw(vec![pe2.clone(); prod.deg() + 1]);
            prod = poly_mod::poly_mod(&(&prod + &bias), &pe) - bias;
            let (quo, rem) = pseudo_div_rem_bigint(&a, &prod);
            if !rem.is_zero() {
                continue;
            }
            result.push((prod, 1));
            a = quo;
            for i in (0..lifted.len()).rev() {
                if (bits & 1 << i) != 0 {
                    lifted.remove(i);
                }
            }
            continue 'outer;
        }
        d += 1;
    }
    result.push((a, 1));
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorize_works_0() {
        let a = Polynomial::<BigInt>::from_raw(vec![(-3).into(), 2.into(), 1.into()]);
        let result = factorize(&a);
        // X^2+2X-3 = (X+3)(X-1)
        assert_eq!(result.len(), 2);
        let factor1: Polynomial<BigInt> = Polynomial::from_raw(vec![3.into(), 1.into()]);
        let factor2: Polynomial<BigInt> = Polynomial::from_raw(vec![(-1).into(), 1.into()]);
        assert!(
            result == vec![(factor1.clone(), 1), (factor2.clone(), 1)]
                || result == vec![(factor2, 1), (factor1, 1)]
        );
    }
}
