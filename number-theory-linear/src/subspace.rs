use num::rational::Ratio;
use num::traits::{Inv, NumAssign};
use num::{Integer, One, Zero};

#[derive(Eq, PartialEq, Clone, Debug)]
pub enum IIMError {
    LinearlyDependent,
    NotInImage,
}

#[derive(Eq, PartialEq, Clone, Debug)]
pub enum SupplementError {
    InsufficientRank,
}

/// Algorithm 2.3.5 (Inverse Image Matrix) in [Cohen].
///
/// mmat: n * m
///
/// vmat: r * m
///
/// Finds xmat: r * n s.t. vmat = xmat mmat or reports that such xmat does not exist.
///
/// Note that matrices are transposed in this implementation with respect to the algorithm in [Cohen].
///
/// [Cohen]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
#[allow(clippy::needless_range_loop)]
pub fn iim<Int: Clone + Integer + NumAssign>(
    mmat: &[Vec<Ratio<Int>>],
    vmat: &[Vec<Ratio<Int>>],
) -> Result<Vec<Vec<Ratio<Int>>>, IIMError> {
    let n = mmat.len();
    let r = vmat.len();
    let m = mmat[0].len();
    let mut mmat = mmat.to_vec();
    let mut bmat = vmat.to_vec();
    assert_eq!(vmat[0].len(), m);
    for j in 0..n {
        let i = if let Some(i) = (j..m).find(|&i| !mmat[j][i].is_zero()) {
            i
        } else {
            return Err(IIMError::LinearlyDependent);
        };
        if i > j {
            for k in 0..n {
                mmat[k].swap(i, j);
            }
            for k in 0..r {
                bmat[k].swap(i, j);
            }
        }
        let d = mmat[j][j].clone().inv();
        let mut c = vec![Ratio::<Int>::zero(); m];
        for k in j + 1..m {
            c[k] = &d * &mmat[j][k];
        }
        for l in j + 1..n {
            for k in j + 1..m {
                let tmp = &c[k] * &mmat[l][j];
                mmat[l][k] -= tmp;
            }
        }
        for l in 0..r {
            for k in j + 1..m {
                let tmp = &c[k] * &bmat[l][j];
                bmat[l][k] -= tmp;
            }
        }
    }
    // Step 6: Solve triangular system
    // mmat is lower-triangular
    let mut xmat = vec![vec![Ratio::<Int>::zero(); n]; r];
    for i in (0..n).rev() {
        for k in 0..r {
            let mut tmp = bmat[k][i].clone();
            for j in i + 1..n {
                tmp -= &mmat[j][i] * &xmat[k][j];
            }
            xmat[k][i] = tmp / &mmat[i][i];
        }
    }
    // Step 7: Check rest of matrix
    for k in n..m {
        for i in 0..r {
            let mut xmsum = Ratio::<Int>::zero();
            for j in 0..n {
                xmsum += &xmat[i][j] * &mmat[j][k]
            }
            if xmsum != bmat[i][k] {
                return Err(IIMError::NotInImage);
            }
        }
    }
    Ok(xmat)
}

/// Algorithm 2.3.6 (Supplement a Basis) in [Cohen].
///
/// mmat: k * n
#[allow(clippy::needless_range_loop)]
pub fn supplement_basis<Int: Clone + Integer + NumAssign>(
    mmat: &[Vec<Ratio<Int>>],
) -> Result<Vec<Vec<Ratio<Int>>>, SupplementError> {
    let k = mmat.len();
    let n = mmat[0].len();
    let mmat_orig = mmat;
    let mut mmat = mmat_orig.to_vec();
    let mut bmat = vec![vec![Ratio::<Int>::zero(); n]; n];
    for i in 0..n {
        bmat[i][i] = Ratio::<Int>::one();
    }
    for s in 0..k {
        let mut t = n;
        for i in s..n {
            if !mmat[s][i].is_zero() {
                t = i;
                break;
            }
        }
        if t >= n {
            return Err(SupplementError::InsufficientRank);
        }
        let d = mmat[s][t].clone().inv();
        bmat[t] = bmat[s].clone();
        bmat[s] = mmat_orig[s].clone();
        for j in s + 1..k {
            mmat[j].swap(s, t);
            let coef_js = &mmat[j][s] * &d;
            for i in 0..n {
                if i != s && i != t {
                    let tmp = &mmat[s][i] * &coef_js;
                    mmat[j][i] -= tmp;
                }
            }
        }
    }
    Ok(bmat)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iim_works_0() {
        let mmat = vec![vec![1.into(), 0.into()], vec![2.into(), 1.into()]];
        let vmat = vec![vec![3.into(), 1.into()]];
        // vmat = (1 1) mmat
        assert_eq!(iim(&mmat, &vmat), Ok(vec![vec![1.into(), 1.into()]]));
    }

    #[test]
    fn iim_works_1() {
        let mmat = vec![vec![1.into(), 0.into()], vec![2.into(), 0.into()]];
        let vmat = vec![vec![3.into(), 1.into()]];
        // mmat's rows are linearly dependent: 2(1 0) - (2 0) = 0
        assert_eq!(iim(&mmat, &vmat), Err(IIMError::LinearlyDependent));
    }

    #[test]
    fn iim_works_2() {
        let mmat = vec![
            vec![1.into(), 0.into(), 1.into()],
            vec![2.into(), 0.into(), 3.into()],
        ];
        let vmat = vec![vec![3.into(), 1.into(), 4.into()]];
        // There is no matrix xmat s.t. vmat = xmat mmat
        assert_eq!(iim(&mmat, &vmat), Err(IIMError::NotInImage));
    }

    #[test]
    fn iim_works_3() {
        let mmat = vec![
            vec![1.into(), 0.into(), 1.into()],
            vec![2.into(), 0.into(), 3.into()],
        ];
        let vmat = vec![vec![1.into(), 0.into(), 1.into()]];
        // vmat = (1 0) mmat
        assert_eq!(iim(&mmat, &vmat), Ok(vec![vec![1.into(), 0.into()]]));
        let vmat = vec![vec![1.into(), 0.into(), (-1).into()]];
        // vmat = (5 -2) mmat
        assert_eq!(iim(&mmat, &vmat), Ok(vec![vec![5.into(), (-2).into()]]));
    }

    #[test]
    fn supplement_basis_0() {
        let mmat = vec![
            vec![1.into(), 0.into(), 1.into()],
            vec![2.into(), 0.into(), 3.into()],
        ];
        let result = supplement_basis(&mmat).unwrap();
        assert_eq!(result[..2], mmat);
    }

    #[test]
    fn supplement_basis_1() {
        let mmat = vec![vec![1.into(), 0.into()], vec![2.into(), 1.into()]];
        assert_eq!(supplement_basis(&mmat), Ok(mmat));
    }

    #[test]
    fn supplement_basis_2() {
        let mmat = vec![
            vec![1.into(), 0.into(), 1.into()],
            vec![2.into(), 0.into(), 2.into()],
        ];
        let result = supplement_basis(&mmat);
        assert_eq!(result, Err(SupplementError::InsufficientRank));
    }
}
