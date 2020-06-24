use num::traits::Pow;
use num::{BigInt, BigRational, One, Zero};

use crate::algebraic::Algebraic;
use crate::hnf::HNF;
use crate::order::{index, Order};
use crate::polynomial::Polynomial;

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

    // phi(w_i)
    let mut phiw = vec![];
    for i in 0..deg {
        let val = create_num(&o.basis[i], &theta);
        phiw.push(val.pow(pow.clone()));
    }
    let mut lcm = BigInt::one();
    for i in 0..deg {
        for j in 0..deg {
            let val = o.basis[i][j].denom().clone();
            lcm = num::integer::lcm(lcm, val);
        }
    }
    // I_p + pO
    let mut basis = vec![vec![BigInt::from(0); deg]; 2 * deg];
    for i in 0..deg {
        for j in 0..deg {
            let val = phiw[i].expr.coef_at(j) * &lcm.clone();
            assert!(val.is_integer());
            basis[i][j] = val.to_integer();
        }
    }
    for i in 0..deg {
        for j in 0..deg {
            let val = &o.basis[i][j] * &lcm;
            if !val.is_integer() {
                eprintln!("val = {}, lcm = {}", val, lcm);
            }
            assert!(val.is_integer());
            let val = val.to_integer();
            basis[i + deg][j] = &val * p;
        }
    }
    // I_p in terms of O's basis
    let mut hnf = HNF::hnf(&HNF::kernel(&basis));
    for row in hnf.0.iter_mut() {
        row.truncate(deg);
    }
    // transformation: I_p in terms of Q(theta) * lcm
    let i_p_len = hnf.0.len();
    let mut i_p = vec![vec![BigInt::zero(); deg]; i_p_len];
    for i in 0..i_p_len {
        for j in 0..deg {
            for k in 0..deg {
                let value = &o.basis[j][k] * &lcm;
                assert!(value.is_integer());
                i_p[i][k] += &hnf.0[i][j] * &value.to_integer();
            }
        }
    }

    // U_p
    let mut u_p = HNF::hnf(&i_p);

    for i in 0..i_p_len {
        let ei = create_num_int(&i_p[i], &theta);
        let mut tmp_basis = vec![vec![BigInt::zero(); deg]; u_p.0.len() + i_p_len];
        for i in 0..u_p.0.len() {
            let prod = &create_num_int(&u_p.0[i], &theta) * &ei;
            for j in 0..deg {
                let value = &prod.expr.coef_at(j) / &lcm;
                assert!(value.is_integer());
                tmp_basis[i][j] = value.to_integer();
            }
        }
        for i in 0..i_p_len {
            for j in 0..deg {
                tmp_basis[i + u_p.0.len()][j] = &i_p[i][j] * p;
            }
        }
        // new_u_p is in terms of U_p + pI_p
        let mut new_u_p = HNF::hnf(&HNF::kernel(&tmp_basis));
        for row in new_u_p.0.iter_mut() {
            row.truncate(u_p.0.len());
        }
        // In terms of theta
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
    // U_p
    let mut new_o_basis = vec![vec![BigInt::zero(); deg]; u_p.0.len() + deg];
    for i in 0..u_p.0.len() {
        new_o_basis[i].clone_from_slice(&u_p.0[i]);
    }
    for i in 0..deg {
        for j in 0..deg {
            new_o_basis[u_p.0.len() + i][j] = (&o.basis[i][j] * &lcm * p).to_integer();
        }
    }
    let u_p = HNF::hnf(&new_o_basis);
    assert_eq!(u_p.0.len(), deg);

    let mut new_o_basis = vec![vec![BigRational::zero(); deg]; deg];
    for i in 0..deg {
        for j in 0..deg {
            new_o_basis[i][j] = BigRational::new(u_p.0[i][j].clone(), &lcm * p);
        }
    }
    let new_o = Order { basis: new_o_basis };
    let mut index = index(&new_o, &o);
    o = new_o;
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
fn create_num_int(a: &[BigInt], theta: &Algebraic) -> Algebraic {
    Algebraic {
        min_poly: theta.min_poly.clone(),
        expr: Polynomial::from_raw(a.iter().map(|x| BigRational::from(x.clone())).collect()),
    }
}
