//! Computes the Hermite normal form (HNF) of a given matrix.
use num::{BigInt, Signed, Zero};
use std::cmp::{max, min};
use std::fmt::Display;

/// A matrix guaranteed to be in HNF.
///
/// Unlike the convention in [Cohen], an HNF is represented as a *lower*-triangular matrix and constructed by using only elementary *row* operations.
/// This convention is justified as follows: it is natural to treat an HNF as a sequence of basis vectors.
/// Sequences of vectors are indexed as a[i][j]; in this indexing,
/// a[i] is treated as a row vector in the ordinary matrix indexing.
/// Therefore, it is better understood as a sequence of row vectors than a sequence of column vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HNF(pub Vec<Vec<BigInt>>);

impl HNF {
    /// Construct the sum module of the given two modules. The resulting module is returned as an HNF.
    pub fn union(a: &HNF, b: &HNF) -> HNF {
        assert_eq!(a.0[0].len(), b.0[0].len());
        let na = a.0.len();
        let nb = b.0.len();
        let m = a.0[0].len();
        // Concatenate and find HNF.
        let mut mat = vec![vec![BigInt::zero(); m]; na + nb];
        #[allow(clippy::needless_range_loop)]
        for i in 0..na {
            mat[i].clone_from_slice(&a.0[i]);
        }
        #[allow(clippy::needless_range_loop)]
        for i in 0..nb {
            mat[i + na].clone_from_slice(&b.0[i]);
        }
        hnf(&mat)
    }
}

impl Display for HNF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
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
    let n = a.len(); // #rows
    let m = a[0].len(); // #columns
    let zero = BigInt::zero();
    let mut k = n - 1;
    let l = max(m, n) - n;
    for i in (l..m).rev() {
        loop {
            // Step 2: Row finished?
            let allzero = (0..k).all(|j| a[j][i] == BigInt::zero());
            if allzero {
                if a[k][i] < BigInt::zero() {
                    // A_k = -A_k
                    for entry in a[k].iter_mut() {
                        *entry *= -1;
                    }
                }
                break;
            } else {
                // Step 3: Choose non-zero entry
                let ind = (0..k).position(|j| a[j][i] != zero).unwrap();
                let mut mi = (a[ind][i].clone(), ind);
                #[allow(clippy::needless_range_loop)]
                for j in ind + 1..k + 1 {
                    if a[j][i] != zero {
                        mi = min(mi, (a[j][i].abs(), j));
                    }
                }
                let j0 = mi.1;
                a.swap(j0, k);
                // Step 4: Reduce
                let b = a[k][i].clone();
                assert_ne!(b, zero);
                for j in 0..k {
                    let q = floor_div(&a[j][i], &b);
                    // A_j -= q * A_k
                    #[allow(clippy::needless_range_loop)]
                    for u in 0..m {
                        let val = &a[k][u] * &q;
                        a[j][u] -= val;
                    }
                }
            }
        }
        // Step 5: Final reductions
        if a[k][i] == zero {
            k += 1;
        } else {
            let b = a[k][i].clone();
            for j in k + 1..n {
                let q = floor_div(&a[j][i], &b);
                // A_j -= q * A_k
                #[allow(clippy::needless_range_loop)]
                for u in 0..m {
                    let val = &a[k][u] * &q;
                    a[j][u] -= val;
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
    let w = a[k..].to_vec();
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
            vec![vec![2.into(), 0.into()], vec![1.into(), 1.into()]],
        )
    }

    #[test]
    fn hnf_works_with_zero_columns() {
        // 0 0
        // 3 1
        // 1 1
        let a: Vec<Vec<BigInt>> = vec![
            vec![0.into(), 0.into()],
            vec![3.into(), 1.into()],
            vec![1.into(), 1.into()],
        ];
        let hnf = hnf(&a);
        assert_eq!(
            hnf.0,
            vec![vec![2.into(), 0.into()], vec![1.into(), 1.into()]],
        );
    }
}
