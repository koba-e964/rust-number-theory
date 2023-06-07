use num::rational::Ratio;
use num::traits::NumAssign;
use num::{Integer, Zero};

use crate::MatrixNotInvertible;

/// Solves ax = b.
///
/// If a is not invertible, this function returns Err(MatrixNotInvertible).
///
/// Complexity: O(n^3)
pub fn solve_linear_system<Int: Clone + Integer + NumAssign>(
    a: &[Vec<Ratio<Int>>],
    b: &[Ratio<Int>],
) -> Result<Vec<Ratio<Int>>, MatrixNotInvertible> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    let n = a.len();
    assert_eq!(b.len(), n);

    // Performs a simple Gaussian elimination process.
    let mut col = 0;
    for row in 0..n {
        let mut nxt = n;
        for i in col..n {
            if a[row][i] != Ratio::zero() {
                nxt = i;
                break;
            }
        }
        if nxt == n {
            return Err(MatrixNotInvertible);
        }
        for row in a.iter_mut() {
            row.swap(col, nxt);
        }
        b.swap(col, nxt);
        let arc = a[row][col].clone();
        for row in a.iter_mut() {
            row[col] /= &arc;
        }
        b[col] /= &arc;
        for i in 0..n {
            if i == col {
                continue;
            }
            let coef = a[row][i].clone();
            for row in a.iter_mut() {
                let val = &coef * &row[col];
                row[i] -= val;
            }
            let val = &coef * &b[col];
            b[i] -= val;
        }
        col += 1;
    }
    Ok(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    use num::One;

    #[test]
    fn solve_linear_system_works() {
        let one = Ratio::<i64>::one();
        let a: Vec<Vec<Ratio<i64>>> = vec![vec![one, one * 2], vec![one * 3, one * 4]];
        let b = vec![one * 5, one * 8];
        let ans = solve_linear_system(&a, &b).unwrap();
        // (2 1) (1 2; 3 4) = (5 8)
        assert_eq!(ans, vec![one * 2, one]);
    }
}
