use num::{BigInt, BigRational, One, Signed, Zero};

use crate::ideal::{FracIdeal, Ideal};
use number_theory_linear::{determinant, hnf::HNF, matrix};

/// Multiplication table of a ring of integers (or orders).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultTable {
    table: Vec<Vec<Vec<BigInt>>>,
}

impl MultTable {
    pub fn new(table: Vec<Vec<Vec<BigInt>>>) -> Self {
        MultTable { table }
    }
    pub fn deg(&self) -> usize {
        self.table.len()
    }
    #[allow(clippy::needless_range_loop)]
    pub fn mul(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        debug_assert_eq!(a.len(), b.len());
        let n = a.len();
        debug_assert_eq!(n, self.deg());
        let mut result = vec![0.into(); n];
        for i in 0..n {
            for j in 0..n {
                let prod = &a[i] * &b[j];
                for k in 0..n {
                    result[k] += &prod * &self.table[i][j][k];
                }
            }
        }
        result
    }

    #[allow(clippy::needless_range_loop)]
    pub fn trace(&self, a: &[BigInt]) -> BigInt {
        let n = self.deg();
        let mut sum = BigInt::zero();
        for i in 0..n {
            for j in 0..n {
                sum += &a[i] * &self.table[j][i][j];
            }
        }
        sum
    }

    #[allow(clippy::needless_range_loop)]
    pub fn norm(&self, a: &[BigInt]) -> BigInt {
        let n = self.deg();
        let mut sum = vec![vec![BigRational::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    sum[j][k] += BigRational::from(&a[i] * &self.table[i][j][k]);
                }
            }
        }
        determinant(&sum).to_integer()
    }

    /// Returns (i, d) s.t. a^{-1} = i / d
    #[allow(clippy::needless_range_loop)]
    pub fn inv(&self, a: &[BigInt]) -> (Vec<BigInt>, BigInt) {
        let norm = self.norm(a).abs();
        let n = self.deg();
        let mut sum = vec![vec![BigRational::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    sum[j][k] += BigRational::from(&a[i] * &self.table[i][j][k]);
                }
            }
        }
        let inv = matrix::inv(&sum).unwrap();
        let mut ans = vec![BigInt::zero(); n];
        // Assuming w[0] is always 1.
        for i in 0..n {
            ans[i] = (&inv[0][i] * &norm).to_integer()
        }
        (ans, norm)
    }

    /// Computes the inverse of the different.
    #[allow(clippy::needless_range_loop)]
    pub fn get_inv_diff(&self) -> FracIdeal<'_> {
        let n = self.deg();
        let mut tr_mat = vec![vec![BigRational::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                let v = &self.table[i][j];
                tr_mat[i][j] = self.trace(v).into();
            }
        }
        let d = matrix::inv(&tr_mat).unwrap();
        let mut denom_lcm = BigInt::one();

        for i in 0..n {
            for j in 0..n {
                denom_lcm = num::integer::lcm(denom_lcm, d[i][j].denom().clone());
            }
        }
        let mut int = vec![vec![BigInt::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                int[i][j] = (&d[i][j] * &denom_lcm).to_integer();
            }
        }
        FracIdeal::new(denom_lcm, Ideal::new(HNF::new(&int), self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebraic::Algebraic;
    use crate::order::Order;
    use crate::polynomial::Polynomial;
    use num::Signed;

    #[test]
    fn mult_table_works() {
        // Multiplication table for Z[i]
        let table = vec![
            vec![vec![1.into(), 0.into()], vec![0.into(), 1.into()]],
            vec![vec![0.into(), 1.into()], vec![(-1).into(), 0.into()]],
        ];
        let table = MultTable::new(table);
        let a = vec![2.into(), 3.into()];
        let b = vec![4.into(), 1.into()];
        let prod = table.mul(&a, &b);
        assert_eq!(prod, vec![5.into(), 14.into()]); // (2+3i) * (4+i) = 5 + 14i
    }

    #[test]
    fn inv_works_1() {
        // Z[sqrt(-5)]
        let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        // 1 + sqrt(-5)
        let a = vec![1.into(), 1.into()];
        let (inv, norm) = mult_table.inv(&a);
        assert_eq!(inv, vec![1.into(), (-1).into()]);
        assert_eq!(norm, 6.into());
    }

    #[test]
    fn inv_works_2() {
        // Z[sqrt(5)]
        let p = Polynomial::from_raw(vec![(-5).into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        // 1 + sqrt(5)
        let a = vec![1.into(), 1.into()];
        let (inv, norm) = mult_table.inv(&a);
        // the inverse is (1 - sqrt(5)) / (-4) = (-1 + sqrt(5)) / 4. Even if the norm is negative, the returned denominator is positive.
        assert_eq!(inv, vec![(-1).into(), 1.into()]);
        assert_eq!(norm, 4.into());
    }

    #[test]
    fn get_inv_diff_works() {
        // Z[sqrt(-5)]
        let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        // The inverse of the different: (sqrt(-5)) / 10
        let inv_diff = mult_table.get_inv_diff();
        let den = inv_diff.denom();
        let num = inv_diff.numer().norm();
        assert_eq!(den * den / num, o.discriminant(&theta).abs());
    }

    #[test]
    fn get_inv_diff_works_2() {
        // Q(cbrt(-19))
        let p = Polynomial::from_raw(vec![19.into(), 0.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let o = crate::integral_basis::find_integral_basis(&theta);
        let mult_table = o.get_mult_table(&theta);
        // Because D(Q(cbrt(-19))) = -1083 = -3 * 19^2, inv_diff's norm should be 1/1083.
        let inv_diff = mult_table.get_inv_diff();
        let den = inv_diff.denom();
        let num = inv_diff.numer().norm();
        assert_eq!(den * den * den / num, o.discriminant(&theta).abs());
    }
}
