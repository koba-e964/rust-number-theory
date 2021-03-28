#![allow(clippy::needless_range_loop)]

use num::traits::Pow;
use num::{BigInt, BigRational, One, Zero};
use std::fmt::{Debug, Display};

use crate::algebraic::Algebraic;
use crate::discriminant::discriminant;
use crate::mult_table::MultTable;
use crate::polynomial::Polynomial;
use number_theory_linear::hnf::HNF;
use number_theory_linear::{determinant, gauss_elim};

/// Order. Constructed from n vectors independent over Q.
#[derive(Clone)]
pub struct Order {
    pub basis: Vec<Vec<BigRational>>,
}

impl Order {
    /// Returns an order's degree.
    pub fn deg(&self) -> usize {
        self.basis.len()
    }

    pub fn discriminant(&self, theta: &Algebraic) -> BigInt {
        let deg = theta.min_poly.deg();
        let det = determinant(&self.basis);
        let disc = discriminant(&theta.min_poly);
        let lc = theta.min_poly.coef_at(deg);
        let value = BigRational::from_integer(disc) * &det * &det;
        let value = value / BigRational::from_integer(lc.pow(2 * (deg - 1)));
        assert!(value.is_integer());
        value.to_integer()
    }

    /// Returns Z[theta].
    pub fn singly_gen(theta: &Algebraic) -> Self {
        let deg = theta.deg();
        let min_poly = &theta.min_poly;
        let mut basis = vec![vec![BigRational::zero(); deg]; deg];
        let mut cur = Algebraic::new_const(min_poly.clone(), BigRational::one());
        for row in basis.iter_mut() {
            row[..cur.expr.dat.len()].clone_from_slice(&cur.expr.dat);
            cur = &cur * theta;
        }
        Order { basis }
    }

    /// Reduce this order to HNF.
    pub fn hnf_reduce(&self) -> Order {
        let mut lcm = BigInt::one();
        let deg = self.basis.len();
        for row in &self.basis {
            for elem in row {
                lcm = num::integer::lcm(lcm, elem.denom().clone());
            }
        }
        let mut basis = vec![vec![BigInt::zero(); deg]; deg];
        for i in 0..deg {
            for j in 0..deg {
                basis[i][j] = (&self.basis[i][j] * &lcm).to_integer();
            }
        }
        let hnf = HNF::hnf(&basis);
        let mut result = vec![vec![BigRational::zero(); deg]; deg];
        for i in 0..deg {
            for j in 0..deg {
                result[i][j] = BigRational::new(hnf.0[i][j].clone(), lcm.clone());
            }
        }
        Order { basis: result }
    }

    pub fn get_mult_table(&self, theta: &Algebraic) -> MultTable {
        let deg = self.deg();
        let mut table = vec![vec![vec![BigInt::zero(); deg]; deg]; deg];
        // This code snipped is copy-pasted from round2.
        // TODO: unify
        for i in 0..deg {
            let oi = Self::create_num(&self.basis[i], &theta);
            for j in 0..deg {
                let oj = Self::create_num(&self.basis[j], &theta);
                let prod = &oi * &oj;
                let mut b = vec![BigRational::zero(); deg];
                for k in 0..deg {
                    b[k] = prod.expr.coef_at(k);
                }
                let inv = gauss_elim(&self.basis, &b).expect("O is not linearly independent");
                for k in 0..deg {
                    assert!(inv[k].is_integer());
                    table[i][j][k] = inv[k].to_integer();
                }
            }
        }
        MultTable::new(table)
    }
    // This code snippet is copy-pasted from round2.
    // TODO: unify
    fn create_num(a: &[BigRational], theta: &Algebraic) -> Algebraic {
        Algebraic {
            min_poly: theta.min_poly.clone(),
            expr: Polynomial::from_raw(a.to_vec()),
        }
    }

    /// Converts an Algebraic (represented as a polynomial of theta) to a coefficient vector (with this Z-basis).
    pub fn to_z_basis(&self, a: &Algebraic) -> Vec<BigRational> {
        let deg = self.deg();
        let mut b = vec![BigRational::zero(); deg];
        for k in 0..deg {
            b[k] = a.expr.coef_at(k);
        }
        gauss_elim(&self.basis, &b).expect("O is not linearly independent")
    }

    pub fn to_z_basis_int(&self, a: &Algebraic) -> Vec<BigInt> {
        let deg = self.deg();
        let inv = self.to_z_basis(a);
        let mut returned = vec![BigInt::zero(); deg];
        for i in 0..deg {
            assert!(inv[i].is_integer());
            returned[i] = inv[i].clone().to_integer();
        }
        returned
    }
}

impl Display for Order {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        let &Order { basis } = &self;
        if basis.is_empty() {
            return write!(f, "()");
        }
        let n = basis[0].len();
        for row in basis {
            for j in 0..n {
                write!(f, "{}{}", row[j], if j + 1 == n { "\n" } else { " " })?;
            }
        }
        Ok(())
    }
}

impl Debug for Order {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

/// For an order a and its suborder b, computes (a : b).
pub fn index(a: &Order, b: &Order) -> BigInt {
    let quot = determinant(&b.basis) / determinant(&a.basis);
    if !quot.is_integer() {
        panic!(
            "Error: quot is not an integer: {} = {} / {}",
            quot,
            determinant(&b.basis),
            determinant(&a.basis),
        );
    }
    quot.to_integer()
}

/// theta's minimal polynomial should be monic.
pub fn trivial_order_monic(theta: &Algebraic) -> Order {
    let deg = theta.deg();
    let mut basis = vec![vec![BigRational::zero(); deg]; deg];
    for i in 0..deg {
        basis[i][i] = BigRational::one();
    }
    Order { basis }
}

/// Returns Z[theta] \cap Z[1 / theta].
pub fn non_monic_initial_order(theta: &Algebraic) -> Order {
    let deg = theta.deg();
    let mut basis = vec![vec![BigRational::zero(); deg]; deg];
    basis[0][0] = BigRational::one();
    for i in 1..deg {
        for j in 1..i + 1 {
            basis[i][j] = BigRational::from_integer(theta.min_poly.coef_at(deg - (i - j)));
        }
    }
    Order { basis }
}

/// Creates the minimum order that contains both a and b.
pub fn union(a: &Order, b: &Order) -> Order {
    assert_eq!(a.basis[0].len(), b.basis[0].len());
    let m = a.basis[0].len();
    let na = a.basis.len();
    let nb = b.basis.len();
    // Concatenate and find HNF.
    // lcm of denominators in a and b.
    let mut lcm = BigInt::one();
    for rowa in &a.basis {
        for elem in rowa {
            let den = elem.denom();
            lcm = num::integer::lcm(lcm, den.clone());
        }
    }
    for rowb in &b.basis {
        for elem in rowb {
            let den = elem.denom();
            lcm = num::integer::lcm(lcm, den.clone());
        }
    }
    let mut basis_a = vec![vec![BigInt::zero(); m]; na];
    let mut basis_b = vec![vec![BigInt::zero(); m]; nb];
    for i in 0..na {
        for j in 0..m {
            let val = (&a.basis[i][j] * &lcm).to_integer();
            basis_a[i][j] = val;
        }
    }
    for i in 0..nb {
        for j in 0..m {
            let val = (&b.basis[i][j] * &lcm).to_integer();
            basis_b[i][j] = val;
        }
    }
    let hnf = HNF::union(&HNF(basis_a), &HNF(basis_b));
    let n = hnf.0[0].len();
    let mut neword = vec![vec![BigRational::zero(); m]; n];
    for i in 0..n {
        for j in 0..m {
            neword[i][j] = BigRational::new(hnf.0[i][j].clone(), lcm.clone());
        }
    }
    Order { basis: neword }
}

#[cfg(test)]
mod tests1 {
    use super::*;
    use crate::algebraic::Algebraic;
    use crate::polynomial::Polynomial;

    /// Returns 1 + 6i.
    fn theta() -> Algebraic {
        Algebraic::new(Polynomial::from_raw(vec![37.into(), (-2).into(), 1.into()]))
    }

    /// For x = 1 + 6i, returns Z[6i] = Z[x].
    fn o() -> Order {
        trivial_order_monic(&theta())
    }
    /// For x = 1 + 6i, returns Z[3i] = Z[(1 + x) / 2].
    fn o1() -> Order {
        let inv2 = BigRational::new(1.into(), 2.into());
        Order {
            basis: vec![
                vec![BigRational::one(), BigRational::zero()],
                vec![inv2.clone(), inv2],
            ],
        }
    }
    /// For x = 1 + 6i, returns Z[2i] = Z[(2 + x) / 3].
    fn o2() -> Order {
        let inv3 = BigRational::new(1.into(), 3.into());
        Order {
            basis: vec![
                vec![BigRational::one(), BigRational::zero()],
                vec![inv3.clone() * BigInt::from(2), inv3],
            ],
        }
    }

    #[test]
    fn index_works() {
        assert_eq!(index(&o1(), &o()), 2.into());
        assert_eq!(index(&o2(), &o()), 3.into());
    }

    #[test]
    fn union_works() {
        let union = union(&o1(), &o2());
        assert_eq!(index(&union, &o()), 6.into());
    }

    #[test]
    fn discriminant_works() {
        assert_eq!(o().discriminant(&theta()), (-144i64).into());
        assert_eq!(o1().discriminant(&theta()), (-36i64).into());
        assert_eq!(o2().discriminant(&theta()), (-16i64).into());
        assert_eq!(union(&o1(), &o2()).discriminant(&theta()), (-4i64).into());
    }

    #[test]
    fn get_mult_table_works() {
        let mult_table = o().get_mult_table(&theta());
        assert_eq!(
            mult_table.mul(&[1.into(), 1.into()], &[1.into(), 1.into()]),
            vec![(-36).into(), 4.into()]
        );
    }
}

#[cfg(test)]
mod tests2 {
    use super::*;
    use crate::algebraic::Algebraic;
    use crate::polynomial::Polynomial;

    fn min_poly() -> Polynomial<BigInt> {
        Polynomial::from_raw(vec![
            5.into(),
            6.into(),
            (-7).into(),
            6.into(),
            (-7).into(),
            6.into(),
        ])
    }

    /// Returns a root of f(x) := 6x^5 - 7x^4 + 6x^3 - 7x^2 + 6x + 5.
    fn theta() -> Algebraic {
        Algebraic::new(min_poly())
    }

    /// Let theta be a root of f(x) := 6x^5 - 7x^4 + 6x^3 - 7x^2 + 6x + 5.
    /// Returns Z[theta] \cap Z[1 / theta].
    fn o() -> Order {
        non_monic_initial_order(&theta())
    }
    fn o1() -> Order {
        let o2_old: Vec<_> = vec![
            vec![1, 0, 0, 0, 0],
            vec![0, 6, 0, 0, 0],
            vec![0, 5, 6, 0, 0],
        ]
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|x| BigRational::new(x.into(), 1.into()))
                .collect()
        })
        .collect();
        let o2_new: Vec<_> = vec![vec![1, 9, 11, 6, 0], vec![1, 2, 10, 5, 6]]
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|x| BigRational::new(x.into(), 2.into()))
                    .collect()
            })
            .collect();
        Order {
            basis: {
                let mut v = o2_old;
                let mut o2_new = o2_new;
                v.append(&mut o2_new);
                v
            },
        }
    }
    #[test]
    fn index_works() {
        assert_eq!(index(&o1(), &o()), 4.into());
    }
    #[test]
    fn discriminant_works() {
        assert_eq!(o().discriminant(&theta()), 9851980752i64.into());
    }
    #[test]
    fn singly_gen_works() {
        let obig = Order::singly_gen(&(&theta() * &Algebraic::from_int(min_poly(), 6)));
        let o = non_monic_initial_order(&theta());
        eprintln!("D(Obig) = {}", obig.discriminant(&theta()));
        assert_eq!(index(&o, &obig), BigInt::from(6).pow(6u32));
    }
}
