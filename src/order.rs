use num::{BigInt, BigRational, One, Zero};

use algebraic::Algebraic;
use determinant::determinant;
use hnf;

/// Order. Constructed from n vectors independent over Q.
#[derive(Clone, Debug)]
pub struct Order {
    pub basis: Vec<Vec<BigRational>>,
}

impl Order {
    /// Returns an order's degree.
    pub fn deg(&self) -> usize {
        self.basis.len()
    }
}

/// For an order a and its suborder b, computes (a : b).
pub fn index(a: &Order, b: &Order) -> BigInt {
    (determinant(&b.basis) / determinant(&a.basis)).to_integer()
}

pub fn trivial_order(theta: &Algebraic) -> Order {
    let deg = theta.deg();
    let mut basis = vec![vec![BigRational::zero(); deg]; deg];
    #[allow(clippy::needless_range_loop)]
    for i in 0..deg {
        basis[i][i] = BigRational::one();
    }
    Order { basis }
}

/// Creates the minimum order that contains both a and b.
/// TODO
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
    let mut mat = vec![vec![BigInt::zero(); na + nb]; m];
    #[allow(clippy::needless_range_loop)]
    for j in 0..na {
        for i in 0..m {
            let val = (&a.basis[i][j] * &lcm).to_integer();
            mat[j][i] = val;
        }
    }
    #[allow(clippy::needless_range_loop)]
    for j in 0..nb {
        for i in 0..m {
            let val = (&b.basis[i][j] * &lcm).to_integer();
            mat[j][i + na] = val;
        }
    }
    let hnf = hnf::hnf(&mat);
    let n = hnf.0[0].len();
    let mut neword = vec![vec![BigRational::zero(); m]; n];
    for i in 0..m {
        #[allow(clippy::needless_range_loop)]
        for j in 0..n {
            neword[j][i] = BigRational::new(hnf.0[i][j].clone(), lcm.clone());
        }
    }
    Order { basis: neword }
}

#[cfg(test)]
mod tests1 {
    use super::*;
    use algebraic::Algebraic;
    use polynomial::Polynomial;

    /// For x = 1 + 6i, returns Z[6i] = Z[x].
    fn o() -> Order {
        let theta = Algebraic::new(Polynomial::from_raw(vec![37.into(), (-2).into(), 1.into()]));
        trivial_order(&theta)
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
}

#[cfg(test)]
mod tests2 {
    use super::*;

    /// Let theta be a root of f(x) := 6x^5 - 7x^4 + 6x^3 - 7x^2 + 6x + 5.
    /// Returns Z[6 theta, 5 / theta].
    fn o() -> Order {
        Order {
            basis: vec![
                vec![1, 0, 0, 0, 0],
                vec![0, 6, 0, 0, 0],
                vec![0, 5, 6, 0, 0],
                vec![0, 4, 5, 6, 0],
                vec![0, 3, 4, 5, 6],
            ]
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|x| BigRational::new(x.into(), 1.into()))
                    .collect()
            })
            .collect(),
        }
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
}
