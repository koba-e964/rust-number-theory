extern crate num;

use num::{BigRational, One, Zero};

/// Given a vector a of length n, consisting of vectors of length n,
/// returns a's determinant, regarding a as a square matrix.
/// Complexity: O(n^3)
pub fn determinant(a: &[Vec<BigRational>]) -> BigRational {
    let n = a.len();
    let mut a = a.to_vec();
    let mut result = BigRational::one();
    for i in 0..n {
        let mut idx = None;
        for j in i..n {
            if a[j][i] != BigRational::zero() {
                idx = Some(j);
                break;
            }
        }
        let idx = match idx {
            None => return BigRational::zero(),
            Some(idx) => idx,
        };
        a.swap(i, idx);
        for j in i + 1..n {
            let factor = &a[j][i] / &a[i][i];
            for k in i..n {
                let tmp = &factor * &a[i][k];
                a[j][k] -= tmp;
            }
        }
        result *= &a[i][i];
    }
    result
}

#[cfg(test)]
mod tests {
    use super::determinant;
    use num::{BigInt, BigRational};
    fn to_rat(n: i64) -> BigRational {
        let n: BigInt = n.into();
        n.into()
    }
    #[test]
    fn test_determinant() {
        // det((2, -1; 5, -4)) = -3
        let mat = vec![vec![to_rat(2), to_rat(-1)], vec![to_rat(5), to_rat(-4)]];
        assert_eq!(determinant(&mat), to_rat(-3));
    }
    #[test]
    fn test_determinant_2() {
        // det((3, 2, 1; -1, 2, 2; -2, -3, 2)) = 33
        let mat = vec![
            vec![to_rat(3), to_rat(2), to_rat(1)],
            vec![to_rat(-1), to_rat(2), to_rat(2)],
            vec![to_rat(-2), to_rat(-3), to_rat(2)],
        ];
        assert_eq!(determinant(&mat), to_rat(33));
    }
}
