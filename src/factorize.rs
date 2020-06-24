extern crate num;

use num::{BigInt, Zero};

pub fn factorize(n: &BigInt) -> Vec<(BigInt, u64)> {
    assert!(*n >= 1.into());
    // TODO naive
    let mut p: BigInt = 2.into();
    let mut n = n.clone();
    let mut fac = Vec::new();
    while &p * &p <= n {
        let mut e = 0;
        while (&n % &p).is_zero() {
            e += 1;
            n /= &p;
        }
        if e > 0 {
            fac.push((p.clone(), e));
        }
        p += 1;
    }
    if n > 1.into() {
        fac.push((n, 1));
    }
    fac
}

#[cfg(test)]
mod tests {
    use super::factorize;
    #[test]
    fn test_factorize() {
        let mut res = factorize(&10.into());
        res.sort_unstable();
        assert_eq!(res, [(2.into(), 1), (5.into(), 1)]);
    }
    #[test]
    fn test_factorize_2() {
        let mut res = factorize(&36355439941184i64.into());
        res.sort_unstable();
        assert_eq!(
            res,
            [
                (2.into(), 6),
                (7.into(), 1),
                (13.into(), 1),
                (149.into(), 1),
                (41894959.into(), 1)
            ]
        );
    }
}
