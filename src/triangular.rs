use num::{BigInt, BigRational, Zero};

/// TODO: exploit the triangularity of b.
#[allow(clippy::needless_range_loop)]
pub fn mul_inv_from_right_exact(a: &[Vec<BigInt>], b: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    let n = a.len();
    let mut brat = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            brat[i][j] = b[i][j].clone().into();
        }
    }
    let invb = crate::matrix::inv(&brat);
    let mut ans = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = BigRational::zero();
            for k in 0..n {
                sum += &invb[k][j] * &a[i][k];
            }
            assert!(sum.is_integer());
            ans[i][j] = sum.to_integer();
        }
    }
    ans
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mul_inv_from_right_exact_works() {
        let a = vec![vec![10.into(), 0.into()], vec![0.into(), (-2).into()]];
        let b = vec![vec![10.into(), 0.into()], vec![5.into(), 1.into()]];
        let c = mul_inv_from_right_exact(&a, &b);
        let c_expected = vec![vec![1.into(), 0.into()], vec![1.into(), (-2).into()]];
        assert_eq!(c, c_expected);
    }
}
