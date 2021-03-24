use num::rational::Ratio;
use num::traits::NumAssign;
use num::{Integer, Zero};

use crate::{matrix, MatrixNotInvertible};

/// TODO: exploit the triangularity of b.
#[allow(clippy::needless_range_loop)]
pub fn mul_inv_from_right_exact<Int: Clone + Integer + NumAssign>(
    a: &[Vec<Int>],
    b: &[Vec<Int>],
) -> Result<Vec<Vec<Int>>, MatrixNotInvertible> {
    let n = a.len();
    let mut brat = vec![vec![Ratio::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            brat[i][j] = b[i][j].clone().into();
        }
    }
    let invb = matrix::inv(&brat)?;
    let mut ans = vec![vec![Int::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = Ratio::zero();
            for k in 0..n {
                sum += &invb[k][j] * &a[i][k];
            }
            assert!(sum.is_integer());
            ans[i][j] = sum.to_integer();
        }
    }
    Ok(ans)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mul_inv_from_right_exact_works() {
        let a = vec![vec![10, 0], vec![0, -2]];
        let b = vec![vec![10, 0], vec![5, 1]];
        let c = mul_inv_from_right_exact(&a, &b);
        let c_expected = vec![vec![1, 0], vec![1, -2]];
        assert_eq!(c, Ok(c_expected));
    }
}
