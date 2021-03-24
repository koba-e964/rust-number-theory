use num::rational::Ratio;
use num::traits::{Inv, NumAssign};
use num::{Integer, One, Zero};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct MatrixNotInvertible;

#[allow(clippy::needless_range_loop)]
pub fn inv<Int: Clone + Integer + NumAssign>(
    a: &[Vec<Ratio<Int>>],
) -> Result<Vec<Vec<Ratio<Int>>>, MatrixNotInvertible> {
    let n = a.len();
    let mut a = a.to_vec();
    let mut b = vec![vec![Ratio::zero(); n]; n];
    for i in 0..n {
        b[i][i] = Ratio::one();
    }
    for i in 0..n {
        let mut idx = None;
        for j in i..n {
            if a[j][i] != Ratio::zero() {
                idx = Some(j);
                break;
            }
        }
        let idx = match idx {
            None => return Err(MatrixNotInvertible),
            Some(idx) => idx,
        };
        a.swap(i, idx);
        b.swap(i, idx);
        let factor = a[i][i].clone().inv();
        for k in 0..n {
            a[i][k] *= &factor;
            b[i][k] *= &factor;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let factor = &a[j][i].clone();
            for k in 0..n {
                let tmp = factor * &a[i][k];
                a[j][k] -= tmp;
                let tmp = factor * &b[i][k];
                b[j][k] -= tmp;
            }
        }
    }
    Ok(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn inv_works() {
        let a = vec![vec![5.into(), 2.into()], vec![2.into(), 1.into()]];
        let inv = inv(&a);
        let inv_expected = vec![vec![1.into(), (-2).into()], vec![(-2).into(), 5.into()]];
        assert_eq!(inv, Ok(inv_expected));
    }
}
