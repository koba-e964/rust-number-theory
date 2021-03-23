use num::{traits::Inv, BigRational, One, Zero};

#[allow(clippy::needless_range_loop)]
pub fn inv(a: &[Vec<BigRational>]) -> Vec<Vec<BigRational>> {
    let n = a.len();
    let mut a = a.to_vec();
    let mut b = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        b[i][i] = BigRational::one();
    }
    for i in 0..n {
        let mut idx = None;
        for j in i..n {
            if a[j][i] != BigRational::zero() {
                idx = Some(j);
                break;
            }
        }
        let idx = match idx {
            None => panic!(),
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
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::BigInt;
    #[test]
    fn inv_works() {
        let a = vec![
            vec![BigInt::from(5).into(), BigInt::from(2).into()],
            vec![BigInt::from(2).into(), BigInt::from(1).into()],
        ];
        let inv = inv(&a);
        let inv_expected = vec![
            vec![BigInt::from(1).into(), BigInt::from(-2).into()],
            vec![BigInt::from(-2).into(), BigInt::from(5).into()],
        ];
        assert_eq!(inv, inv_expected);
    }
}
