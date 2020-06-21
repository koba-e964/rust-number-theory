//! Computes the Hermite normal form (HNF) of a given matrix.
use num::{BigInt, Signed, Zero};
use std::cmp::{max, min};
use std::fmt::Display;

/// A matrix guaranteed to be in HNF.
/// Matrices are represented as a sequence of rows: a[i][j] represents (i, j)-coefficient of the matrix.
#[derive(Clone, Debug)]
pub struct HNF(pub Vec<Vec<BigInt>>);

impl Display for HNF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        let &HNF(inner) = &self;
        let n = inner[0].len();
        for row in inner {
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                write!(f, "{}{}", row[j], if j + 1 == n { "\n" } else { " " })?;
            }
        }
        Ok(())
    }
}

/// Algorithm 2.4.4 in [Cohen]
/// Computes the HNF of a given matrix.
///
/// [Cohen]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
#[allow(clippy::many_single_char_names)]
pub fn hnf(a: &[Vec<BigInt>]) -> HNF {
    let mut a = a.to_vec();
    // Step 1: Initialize
    let m = a.len(); // #rows
    let n = a[0].len(); // #columns
    let zero = BigInt::zero();
    let mut k = n - 1;
    let l = max(m, n) - n;
    for i in (l..m).rev() {
        loop {
            // Step 2: Row finished?
            let allzero = (0..k).all(|j| a[i][j] == BigInt::zero());
            if allzero {
                if a[i][k] < BigInt::zero() {
                    for j in 0..m {
                        a[i][j] *= -1;
                    }
                }
                break;
            } else {
                // Step 3: Choose non-zero entry
                let ind = (0..k).position(|j| a[i][j] != zero).unwrap();
                let mut mi = (a[i][ind].clone(), ind);
                for j in ind + 1..k + 1 {
                    if a[i][j] != zero {
                        mi = min(mi, (a[i][j].abs(), j));
                    }
                }
                let j0 = mi.1;
                for row in a.iter_mut() {
                    row.swap(j0, k);
                }
                // Step 4: Reduce
                let b = a[i][k].clone();
                assert_ne!(b, zero);
                for j in 0..k {
                    let q = floor_div(&a[i][j], &b);
                    // A_j -= q * A_k
                    for row in a.iter_mut() {
                        let val = &row[k] * &q;
                        row[j] -= val;
                    }
                }
            }
        }
        // Step 5: Final reductions
        if a[i][k] == zero {
            k += 1;
        } else {
            let b = a[i][k].clone();
            for j in k + 1..n {
                let q = floor_div(&a[i][j], &b);
                // A_j -= q * A_k
                for row in a.iter_mut() {
                    let val = &row[k] * &q;
                    row[j] -= val;
                }
            }
        }
        // Step 6: Finished?
        if i == l {
            break;
        }
        k -= 1;
    }
    // Step 6: Finished?
    let mut w = vec![Vec::new(); m];
    for i in 0..m {
        w[i] = a[i][k..n].to_vec();
    }
    HNF(w)
}

/// Computes floor(a / b).
fn floor_div(a: &BigInt, b: &BigInt) -> BigInt {
    let mut q = a / b;
    if a < &(&q * b) {
        q -= 1;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hnf_returns_nonnegative_matrix() {
        // 3 1
        // 1 1
        let a: Vec<Vec<BigInt>> = vec![vec![3.into(), 1.into()], vec![1.into(), 1.into()]];
        let hnf = hnf(&a);
        assert_eq!(
            hnf.0,
            vec![vec![2.into(), 1.into()], vec![0.into(), 1.into()]],
        )
    }

    #[test]
    fn hnf_works_with_zero_columns() {
        // 0 3 1
        // 0 1 1
        let a: Vec<Vec<BigInt>> = vec![
            vec![0.into(), 3.into(), 1.into()],
            vec![0.into(), 1.into(), 1.into()],
        ];
        let hnf = hnf(&a);
        assert_eq!(
            hnf.0,
            vec![vec![2.into(), 1.into()], vec![0.into(), 1.into()]],
        );
    }
}
