use num::{BigInt, FromPrimitive, One, ToPrimitive, Zero};
use std::cmp::max;

// Types used in this algorithm
type K = f64;
type SqNorm = f64; // squared norms
type R = BigInt;

// helper functions
fn inner(a: &[K], b: &[K]) -> K {
    debug_assert_eq!(a.len(), b.len());
    let mut sum: K = 0.0;
    for i in 0..a.len() {
        sum += a[i] * b[i];
    }
    sum
}
fn norm_sqr(a: &[K]) -> SqNorm {
    let mut sum: K = 0.0;
    for i in 0..a.len() {
        sum += a[i] * a[i];
    }
    sum
}
fn to_int(a: K) -> R {
    BigInt::from_f64((a + 0.5).floor()).unwrap()
}
fn to_real(r: R) -> K {
    r.to_f64().unwrap()
}

/// Algorithm 2.6.3 in \[Cohen\]
/// The returned value (lll, h) satisfies lll = h * basis (as matrices).
///
/// \[Cohen\]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
pub fn lll(basis: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<BigInt>>) {
    let n = basis.len();
    let m = basis[0].len();
    let mut k = 1;
    let mut kmax = 0;
    let mut basis = basis.to_vec();
    let mut bstar = basis.clone();
    let mut b = vec![0.0; n];
    let mut mu = vec![vec![0.0; n]; n];
    let mut h = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        h[i][i] = BigInt::one();
    }

    // subroutines
    macro_rules! swap {
        ($k:expr) => {
            let k: usize = $k;
            basis.swap(k, k + 1);
            h.swap(k, k + 1);
            mu.swap(k, k + 1);
            let tmpmu = mu[k][k]; // old mu[k+1][k]
            mu[k][k] = 0.0;
            let tmpb = b[k + 1] + tmpmu * tmpmu * b[k];
            mu[k + 1][k] = tmpmu * b[k] / tmpb;
            let tmpbasis = bstar[k].clone();
            bstar[k] = bstar[k + 1].clone();
            for l in 0..n {
                bstar[k][l] += tmpmu * tmpbasis[l];
            }
            for l in 0..n {
                bstar[k + 1][l] = -mu[k + 1][k] * bstar[k + 1][l] + b[k + 1] / tmpb * tmpbasis[l];
            }
            b[k + 1] = b[k + 1] * b[k] / tmpb;
            b[k] = tmpb;
            for i in k + 2..kmax + 1 {
                let t = mu[i][k + 1];
                mu[i][k + 1] = mu[i][k] - tmpmu * t;
                mu[i][k] = t + mu[k + 1][k] * mu[i][k];
            }
        };
    }

    macro_rules! red {
        ($k:expr, $l:expr) => {
            let k = $k;
            let l = $l;
            if mu[k][l].abs() >= 0.5 {
                let q = to_int(mu[k][l]);
                let qr = to_real(q.clone());
                for u in 0..n {
                    basis[k][u] -= qr * basis[l][u];
                    let tmp = &q * &h[l][u];
                    h[k][u] -= tmp;
                }
                mu[k][l] -= qr;
            }
        };
    }

    b[0] = norm_sqr(&basis[0]);
    loop {
        // Step 2
        if kmax < k {
            kmax = k;
            bstar[k] = basis[k].clone();
            for j in 0..k {
                mu[k][j] = inner(&basis[k], &bstar[j]) / b[j];
                for l in 0..m {
                    bstar[k][l] -= bstar[j][l] * mu[k][j];
                }
            }
            b[k] = norm_sqr(&bstar[k]);
        }
        loop {
            red!(k, k - 1);
            if b[k] < (0.75 - mu[k][k - 1] * mu[k][k - 1]) * b[k - 1] {
                swap!(k - 1);
                k = max(1, k - 1);
            } else {
                break;
            }
        }
        for l in (0..k - 1).rev() {
            red!(k, l);
        }
        k += 1;
        if k >= n {
            break;
        }
    }
    (basis, h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lll_works_1() {
        let basis = vec![
            vec![1.0, 1.0, 1.0],
            vec![-1.0, 0.0, 2.0],
            vec![3.0, 5.0, 6.0],
        ];
        let (lll, h) = lll(&basis);
        // check if lll = h * basis.
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += h[i][k].to_f64().unwrap() * basis[k][j];
                }
                assert!((lll[i][j] - sum).abs() <= 1.0e-9);
            }
        }
    }
}
