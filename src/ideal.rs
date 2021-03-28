use crate::algebraic::Algebraic;
use crate::mult_table::MultTable;
use num::{BigInt, BigRational, One, Zero};
use number_theory_linear::hnf::HNF;
use std::ops::{Add, Mul};

/// An ideal represented by an HNF. Basis is of Z_K (the integral basis), not of K.
/// Specifically, let w_1, ..., w_n be a basis of Z_K. Then this object represents a Z-lattice {a_1 w_1 + ... + a_n w_n},
/// where (a_1, ..., a_n) ranges over all row vectors in hnf. It is required that the resulting Z-lattice is an ideal.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Ideal<'mul> {
    hnf: HNF,
    mult_table: &'mul MultTable,
}

impl<'mul> Ideal<'mul> {
    /// Creates an ideal from an HNF and a multiplication table.
    /// This is a low-level API; for convenient ideal creation you might want to use `principal` and `+`.
    pub fn new(hnf: HNF, mult_table: &'mul MultTable) -> Self {
        Ideal { hnf, mult_table }
    }
    pub fn norm(&self) -> BigInt {
        self.hnf.determinant()
    }
    /// Creates a principal ideal generated by elem.
    pub fn principal(elem: &[BigInt], mult_table: &'mul MultTable) -> Self {
        let deg = mult_table.deg();
        debug_assert_eq!(elem.len(), deg);
        let mut rows = vec![];
        for i in 0..deg {
            let mut wi = vec![BigInt::zero(); deg];
            wi[i] = BigInt::one();
            rows.push(mult_table.mul(elem, &wi));
        }
        Ideal {
            hnf: HNF::hnf(&rows),
            mult_table,
        }
    }

    pub fn deg(&self) -> usize {
        self.mult_table.deg()
    }

    /// Finds a such that (a) = self /\ Z.
    pub fn cap_z(&self) -> BigInt {
        self.hnf.0[0][0].clone()
    }

    /// Given an ideal and the inverse of the different, finds the former's inverse.
    #[allow(clippy::needless_range_loop)]
    pub fn inv(&self, inv_diff: &FracIdeal<'mul>) -> FracIdeal<'mul> {
        let n = self.deg();

        let a = self.cap_z();
        let c = self * inv_diff.numer();

        let mut ab = vec![vec![BigInt::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                ab[i][j] = &a * &inv_diff.numer().hnf.0[i][j];
            }
        }
        let c = c.hnf.as_vecs();
        let d = number_theory_linear::triangular::mul_inv_from_right_exact(&ab, &c).unwrap();

        let mut trd = vec![vec![BigInt::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                trd[i][j] = d[j][i].clone();
            }
        }
        FracIdeal::new(a, Ideal::new(HNF::hnf(&trd), self.mult_table))
    }

    pub fn contains(&self, num: &[BigInt]) -> bool {
        // Solved in an elephant way: checking if this + (num) == this
        // TODO improve
        let num_ideal = Self::principal(num, self.mult_table);
        let new_ideal = self + &num_ideal;
        new_ideal == *self
    }
}

impl<'a, 'mul> Add for &'a Ideal<'mul> {
    type Output = Ideal<'mul>;

    fn add(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(self.mult_table, rhs.mult_table);
        let mut basis_b = rhs.hnf.as_vecs();
        let mut res = self.hnf.as_vecs();
        res.append(&mut basis_b);
        let hnf = HNF::hnf(&res);
        Ideal {
            hnf,
            mult_table: self.mult_table,
        }
    }
}

impl<'a, 'mul> Mul for &'a Ideal<'mul> {
    type Output = Ideal<'mul>;

    /// O(deg^5)
    fn mul(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(self.mult_table, rhs.mult_table);
        let basis_a = self.hnf.as_vecs();
        let basis_b = rhs.hnf.as_vecs();
        let mut res = vec![];
        for v in &basis_a {
            for w in &basis_b {
                let prod = self.mult_table.mul(v, w);
                res.push(prod)
            }
        }
        let hnf = HNF::hnf(&res);
        Ideal {
            hnf,
            mult_table: self.mult_table,
        }
    }
}

/// A fractional ideal.
#[derive(Debug, Clone)]
pub struct FracIdeal<'mul> {
    denom: BigInt,
    numer: Ideal<'mul>,
}

impl<'mul> FracIdeal<'mul> {
    pub fn new(denom: BigInt, numer: Ideal<'mul>) -> Self {
        FracIdeal { denom, numer }
    }

    pub fn denom(&self) -> &BigInt {
        &self.denom
    }

    pub fn numer(&self) -> &Ideal<'mul> {
        &self.numer
    }
}

/// A fractional ideal represented in two-element form.
/// Cf. http://www.kurims.kyoto-u.ac.jp/EMIS/journals/JTNB/2004-1/Belabas.pdf, 6.13
pub struct TwoElementFracIdeal(BigRational, Algebraic);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::order::Order;
    use crate::polynomial::Polynomial;

    #[test]
    fn ideal_norm_test() {
        // Z[sqrt(-5)], (2, 1 + sqrt(-5))
        let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let hnf = HNF::hnf(&[
            vec![1.into(), 1.into()],
            vec![5.into(), 1.into()],
            vec![2.into(), 0.into()],
            vec![0.into(), 2.into()],
        ]);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        let x = Ideal::new(hnf, &mult_table);
        assert_eq!(x.norm(), 2.into());
    }

    #[test]
    fn ideal_mul_test() {
        // Z[sqrt(-5)], (2, 1 + sqrt(-5))
        let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let hnf = HNF::hnf(&[
            vec![1.into(), 1.into()],
            vec![5.into(), 1.into()],
            vec![2.into(), 0.into()],
            vec![0.into(), 2.into()],
        ]);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        let x = Ideal::new(hnf, &mult_table);
        assert_eq!(x.norm(), 2.into());
        let two = &x * &x; // (2, 1 + sqrt(-5))^2 = (2)
        assert_eq!(two, Ideal::principal(&[2.into(), 0.into()], &mult_table));
    }

    #[test]
    fn ideal_inv_test() {
        // Z[sqrt(-5)], (2, 1 + sqrt(-5))
        let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let hnf = HNF::hnf(&[
            vec![1.into(), 1.into()],
            vec![5.into(), 1.into()],
            vec![2.into(), 0.into()],
            vec![0.into(), 2.into()],
        ]);
        let o = Order::singly_gen(&theta);
        let mult_table = o.get_mult_table(&theta);
        let x = Ideal::new(hnf, &mult_table);
        // The inverse of the different: (sqrt(-5)) / 10
        let inv_diff = mult_table.get_inv_diff();
        let x_inv = x.inv(&inv_diff); // (2, 1 + sqrt(-5)) / 2
        assert_eq!(x_inv.numer(), &x);
        assert_eq!(x_inv.denom(), &2.into());
    }
}
