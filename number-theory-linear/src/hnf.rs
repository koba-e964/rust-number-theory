//! Computes the Hermite normal form (HNF) of a given matrix.
use num::{BigInt, One, Signed, Zero};
use std::cmp::min;
use std::fmt::Display;

/// A matrix guaranteed to be in HNF.
///
/// Unlike the convention in [Cohen], an HNF is represented as a *lower*-triangular matrix and constructed by using only elementary *row* operations.
/// This convention is justified as follows: it is natural to treat an HNF as a sequence of basis vectors.
/// Sequences of vectors are indexed as a[i][j]; in this indexing,
/// a[i] is treated as a row vector in the ordinary matrix indexing.
/// Therefore, it is better understood as a sequence of row vectors than a sequence of column vectors.
#[allow(clippy::upper_case_acronyms)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HNF(pub(crate) Vec<Vec<BigInt>>);

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
        HNF::new(&mat)
    }
    /// Algorithm 2.4.4 in [Cohen]
    /// Computes the HNF of a given matrix.
    ///
    /// [Cohen]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
    pub fn new(a: &[Vec<BigInt>]) -> HNF {
        hnf_with_ker(a).0
    }

    /// Returns a collection of row vectors u s.t. uA = 0.
    /// Returned elements are guaranteed to be linearly independent.
    pub fn kernel(a: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
        hnf_with_ker(a).1
    }

    pub fn as_vecs(&self) -> Vec<Vec<BigInt>> {
        self.0.clone()
    }

    pub fn into_vecs(self) -> Vec<Vec<BigInt>> {
        self.0
    }

    pub fn determinant(&self) -> BigInt {
        let mut prod = BigInt::from(1);
        if self.dim() != self.deg() {
            return BigInt::zero();
        }
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.0.len() {
            prod *= &self.0[i][i];
        }
        prod
    }

    pub fn dim(&self) -> usize {
        self.0.len()
    }

    pub fn deg(&self) -> usize {
        if self.dim() == 0 {
            0
        } else {
            self.0[0].len()
        }
    }
}

impl Display for HNF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        let &HNF(inner) = &self;
        if inner.is_empty() {
            return write!(f, "()");
        }
        let n = inner[0].len();
        for row in inner {
            assert_eq!(row.len(), n);
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                write!(f, "{}{}", row[j], if j + 1 == n { "\n" } else { " " })?;
            }
        }
        Ok(())
    }
}

impl AsRef<[Vec<BigInt>]> for HNF {
    fn as_ref(&self) -> &[Vec<BigInt>] {
        &self.0
    }
}

#[allow(clippy::many_single_char_names)]
pub fn hnf_with_ker(a: &[Vec<BigInt>]) -> (HNF, Vec<Vec<BigInt>>) {
    let (w, u, k) = hnf_with_u(a);
    let u = u[..k].to_vec();
    (w, u)
}

/// Algorithm 2.4.4 in [Cohen]
/// Given a n * m matrix A, Computes the HNF B of A, an n * n matrix U s.t. B = UA and dim ker A.
///
/// [Cohen]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
#[allow(clippy::many_single_char_names)]
pub fn hnf_with_u(a: &[Vec<BigInt>]) -> (HNF, Vec<Vec<BigInt>>, usize) {
    if a.is_empty() {
        return (HNF(vec![]), vec![], 0);
    }
    let mut a = a.to_vec();
    // Step 1: Initialize
    let n = a.len(); // #rows
    let m = a[0].len(); // #columns
    let zero = BigInt::zero();
    let mut k = n - 1;
    let mut u = vec![vec![zero.clone(); n]; n];
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        u[i][i] = BigInt::one();
    }
    for i in (0..m).rev() {
        loop {
            // Step 2: Row finished?
            let allzero = (0..k).all(|j| a[j][i] == BigInt::zero());
            if allzero {
                if a[k][i] < BigInt::zero() {
                    // A_k = -A_k
                    // U_k = -U_k
                    for entry in a[k].iter_mut() {
                        *entry *= -1;
                    }
                    for entry in u[k].iter_mut() {
                        *entry *= -1;
                    }
                }
                break;
            } else {
                // Step 3: Choose non-zero entry
                let ind = (0..k).position(|j| a[j][i] != zero).unwrap();
                let mut mi = (a[ind][i].abs(), ind);
                #[allow(clippy::needless_range_loop)]
                for j in ind + 1..k + 1 {
                    if a[j][i] != zero {
                        mi = min(mi, (a[j][i].abs(), j));
                    }
                }
                let j0 = mi.1;
                a.swap(j0, k);
                u.swap(j0, k);
                // Step 4: Reduce
                let b = a[k][i].clone();
                assert_ne!(b, zero);
                for j in 0..k {
                    let q = floor_div(&a[j][i], &b);
                    // A_j -= q * A_k
                    // U_j -= q * U_k
                    #[allow(clippy::needless_range_loop)]
                    for v in 0..m {
                        let val = &a[k][v] * &q;
                        a[j][v] -= val;
                    }
                    #[allow(clippy::needless_range_loop)]
                    for v in 0..n {
                        let val = &u[k][v] * &q;
                        u[j][v] -= val;
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
                // U_j -= q * U_k
                #[allow(clippy::needless_range_loop)]
                for v in 0..m {
                    let val = &a[k][v] * &q;
                    a[j][v] -= val;
                }
                #[allow(clippy::needless_range_loop)]
                for v in 0..n {
                    let val = &u[k][v] * &q;
                    u[j][v] -= val;
                }
            }
        }
        // Step 6: Finished?
        // Modified, because the condition i == l (:= max(m, n) - n) appears to be incorrect.
        if k == 0 || i == 0 {
            break;
        }
        k -= 1;
    }
    // Step 6: Finished?
    let w = a[k..].to_vec();
    (HNF(w), u, k)
}

/// Computes floor(a / b).
fn floor_div(a: &BigInt, b: &BigInt) -> BigInt {
    if b < &BigInt::zero() {
        return floor_div(&(-a), &(-b));
    }
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
        let hnf = HNF::new(&a);
        assert_eq!(
            hnf.0,
            vec![vec![2.into(), 0.into()], vec![1.into(), 1.into()]],
        )
    }

    #[test]
    fn hnf_works_with_zero_rows() {
        // 0 0
        // 3 1
        // 1 1
        let a: Vec<Vec<BigInt>> = vec![
            vec![0.into(), 0.into()],
            vec![3.into(), 1.into()],
            vec![1.into(), 1.into()],
        ];
        let hnf = HNF::new(&a);
        assert_eq!(
            hnf.0,
            vec![vec![2.into(), 0.into()], vec![1.into(), 1.into()]],
        );
    }

    #[test]
    fn hnf_works_with_zero_columns() {
        // 2 0 1 0
        // 1 1 0 0
        let a: Vec<Vec<BigInt>> = vec![
            vec![2.into(), 0.into(), 1.into(), 0.into()],
            vec![1.into(), 1.into(), 0.into(), 0.into()],
        ];
        let hnf = HNF::new(&a);
        assert_eq!(
            hnf.0,
            vec![
                vec![1.into(), 1.into(), 0.into(), 0.into()],
                vec![2.into(), 0.into(), 1.into(), 0.into()],
            ]
        );
    }

    #[test]
    fn hnf_terminates() {
        // 1516
        // -154
        // -336
        // -1423
        let a: Vec<Vec<BigInt>> = vec![
            vec![1516.into()],
            vec![(-154).into()],
            vec![(-336).into()],
            vec![(-1423).into()],
        ];
        let hnf = HNF::new(&a);
        assert_eq!(hnf.0.len(), 1);
    }

    #[test]
    fn kernel_works() {
        // 5 0
        // 7 0
        // 2 0
        let a: Vec<Vec<BigInt>> = vec![
            vec![5.into(), 0.into()],
            vec![7.into(), 0.into()],
            vec![2.into(), 0.into()],
        ];
        let kern = HNF::kernel(&a);
        let kern = HNF::new(&kern);
        // ker A is a 2-dimensional Z-module, spanned by e.g. {(-1 1 -1), (2 0 -5)}.
        // Whatever its "reduction" is, it has two rows, each of which A annihilates.
        assert_eq!(kern.0.len(), 2);
        for i in 0..2 {
            assert_eq!(
                &kern.0[i][0] * 5 + &kern.0[i][1] * 7 + &kern.0[i][2] * 2,
                0.into(),
            );
        }
    }
}
