use num::{BigInt, BigRational, One, Zero};
use std::ops::{AddAssign, Mul, RemAssign};

use crate::algebraic::Algebraic;
use crate::order::{index, Order};
use crate::polynomial::Polynomial;
use number_theory_linear::gauss_elim;
use number_theory_linear::hnf::HNF;

/// Performs round 2 algorithm to an order once.
#[allow(clippy::needless_range_loop)]
pub fn one_step(theta: &Algebraic, o: &Order, p: &BigInt) -> (Order, u64) {
    let mut o = o.clone();
    let deg = theta.deg();
    // I_p
    // pow should be a power of p and >= deg.
    let mut pow = BigInt::one();
    while pow < BigInt::from(deg) {
        pow *= p;
    }

    // Find the multiplication table first.
    // o[i] * o[j] = \sum_k table[i][j][k] * o[k]
    // Usually the content of table is held mod p or mod p^2.
    // Time complexity: O(d^5) operations
    let p2 = p * p;
    let mut table = vec![vec![vec![BigInt::zero(); deg]; deg]; deg];
    let mut table2 = vec![vec![vec![BigInt::zero(); deg]; deg]; deg];
    for i in 0..deg {
        let oi = create_num(&o.basis[i], &theta);
        for j in 0..deg {
            let oj = create_num(&o.basis[j], &theta);
            let prod = &oi * &oj;
            let mut b = vec![BigRational::zero(); deg];
            for k in 0..deg {
                b[k] = prod.expr.coef_at(k);
            }
            let inv = gauss_elim(&o.basis, &b).expect("O is not linearly independent");
            for k in 0..deg {
                assert!(inv[k].is_integer());
                table2[i][j][k] = inv[k].to_integer() % &p2;
                table[i][j][k] = &table2[i][j][k] % p;
            }
        }
    }

    // phi(w_i)
    let mut phiw: Vec<Vec<BigInt>> = vec![];
    for i in 0..deg {
        let mut val: Vec<BigInt> = vec![BigInt::zero(); deg];
        val[i] = BigInt::one();
        phiw.push(pow_mod_p::<BigInt>(&val, &pow, &table, &p));
    }
    // I_p + pO in terms of O's basis
    let mut basis = vec![vec![BigInt::from(0); deg]; 2 * deg];
    for i in 0..deg {
        for j in 0..deg {
            let val = phiw[i][j].clone();
            basis[i][j] = val;
        }
    }
    for i in 0..deg {
        basis[i + deg][i] = p.clone();
    }
    // I_p in terms of O's basis
    let mut i_p = HNF::hnf(&HNF::kernel(&basis));
    let i_p_len = i_p.0.len();
    for row in i_p.0.iter_mut() {
        row.truncate(deg);
    }

    // U_p
    let mut u_p = i_p.clone();

    for i in 0..i_p_len {
        // U_p eta[i] + pI_p in terms of O's basis
        let mut tmp_basis = vec![vec![]; u_p.0.len() + i_p_len];
        for j in 0..u_p.0.len() {
            // Find eta[i] * eta[j] mod p^2
            let prod = mul_mod_p::<BigInt>(&i_p.0[i], &u_p.0[j], &table2, &p2);
            tmp_basis[j] = prod;
        }
        for i in 0..i_p_len {
            tmp_basis[i + u_p.0.len()] = vec![BigInt::zero(); deg];
            for j in 0..deg {
                tmp_basis[i + u_p.0.len()][j] = &i_p.0[i][j] * p;
            }
        }
        // new_u_p is in terms of U_p + pI_p
        let mut new_u_p = HNF::hnf(&HNF::kernel(&tmp_basis));
        for row in new_u_p.0.iter_mut() {
            row.truncate(u_p.0.len());
        }
        // In terms of O's basis
        let mut tmp_basis = vec![vec![BigInt::zero(); deg]; new_u_p.0.len()];
        for i in 0..new_u_p.0.len() {
            for j in 0..u_p.0.len() {
                for k in 0..deg {
                    tmp_basis[i][k] += &u_p.0[j][k] * &new_u_p.0[i][j];
                }
            }
        }
        let new_u_p = HNF::hnf(&tmp_basis);
        u_p = new_u_p;
    }

    assert!(u_p.0.len() <= deg);
    // U_p in O's basis
    let mut new_o_basis = vec![vec![BigInt::zero(); deg]; u_p.0.len() + deg];
    for i in 0..u_p.0.len() {
        new_o_basis[i].clone_from_slice(&u_p.0[i]);
    }
    for i in 0..deg {
        new_o_basis[u_p.0.len() + i][i] = p.clone();
    }
    let u_p = HNF::hnf(&new_o_basis);
    assert_eq!(u_p.0.len(), deg);

    // Basis conversion: new O in terms of Q(theta)'s basis (theta^i)
    // TODO: reduce
    let mut new_o_basis = vec![vec![BigRational::zero(); deg]; deg];
    for i in 0..deg {
        for j in 0..deg {
            for k in 0..deg {
                new_o_basis[i][k] +=
                    BigRational::new(u_p.0[i][j].clone(), p.clone()) * &o.basis[j][k];
            }
        }
    }
    let new_o = Order { basis: new_o_basis };
    let mut index = index(&new_o, &o);
    o = new_o.hnf_reduce();
    let mut howmany = 0;
    while index > BigInt::one() {
        assert_eq!(&index % p, BigInt::zero());
        index /= p;
        howmany += 1;
    }
    (o, howmany)
}

fn create_num(a: &[BigRational], theta: &Algebraic) -> Algebraic {
    Algebraic {
        min_poly: theta.min_poly.clone(),
        expr: Polynomial::from_raw(a.to_vec()),
    }
}

fn pow_mod_p<Int>(a: &[Int], e: &BigInt, table: &[Vec<Vec<Int>>], p: &Int) -> Vec<Int>
where
    Int: AddAssign + Zero + for<'a> RemAssign<&'a Int> + Clone,
    for<'a> &'a Int: Mul<&'a Int, Output = Int>,
{
    // To avoid mentioning the identity element, multiply by a beforhand.
    let mut e = e - 1;
    let mut prod = a.to_vec();
    let mut cur = a.to_vec();
    while e > BigInt::zero() {
        if &e % 2 == BigInt::one() {
            // Type annotation is necessary; if type annotation is missing, compilation will result in
            // overflow evaluating the requirement `&'b num::rational::Ratio<_>: std::ops::Mul`.
            prod = mul_mod_p::<Int>(&prod, &cur, table, p);
        }
        cur = mul_mod_p::<Int>(&cur, &cur, table, p);
        e /= 2;
    }
    prod
}

/// Complexity: O(n^3) operations
#[allow(clippy::needless_range_loop)]
fn mul_mod_p<Int>(a: &[Int], b: &[Int], table: &[Vec<Vec<Int>>], p: &Int) -> Vec<Int>
where
    Int: AddAssign + Zero + for<'a> RemAssign<&'a Int> + Clone,
    for<'a> &'a Int: Mul<&'a Int, Output = Int>,
{
    let n = a.len();
    let mut result = vec![Int::zero(); n];
    for i in 0..n {
        for j in 0..n {
            let coef = &a[i] * &b[j];
            for k in 0..n {
                result[k] += &coef * &table[i][j][k];
            }
        }
    }
    for i in 0..n {
        result[i] %= p;
    }
    result
}
