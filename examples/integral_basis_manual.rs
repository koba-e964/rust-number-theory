#![allow(clippy::needless_range_loop)]

extern crate num;
extern crate rust_number_theory;

use num::traits::Pow;
use num::{BigInt, BigRational, One, Signed, Zero};

use rust_number_theory::algebraic::Algebraic;
use rust_number_theory::hnf::HNF;
use rust_number_theory::order::{index, non_monic_initial_order, Order};
use rust_number_theory::polynomial::Polynomial;

fn factorize(x: &BigInt) -> Vec<(BigInt, u64)> {
    assert!(x > &BigInt::zero());
    let mut x = x.clone();
    let mut v: BigInt = 2.into();
    let mut ans = vec![];
    while &v * &v <= x {
        let mut e = 0;
        while (&x % &v).is_zero() {
            x /= &v;
            e += 1;
        }
        if e > 0 {
            ans.push((v.clone(), e));
        }
        v += 1;
    }
    if x > BigInt::one() {
        ans.push((x, 1));
    }
    ans
}

fn main() {
    let p;
    let case_nr = 1;
    if case_nr == 1 {
        p = Polynomial::from_raw(vec![
            5.into(),
            6.into(),
            (-7).into(),
            6.into(),
            (-7).into(),
            6.into(),
        ]);
    } else if case_nr == 2 {
        p = Polynomial::from_raw(vec![1.into(), 9.into(), 0.into(), 1.into()]);
    } else if case_nr == 3 {
        p = Polynomial::from_raw(vec![37.into(), 2.into(), 1.into()]);
    } else {
        panic!();
    }
    let theta = Algebraic::new(p);
    find_integral_basis(&theta);
}

fn find_integral_basis(theta: &Algebraic) {
    let deg = theta.deg();
    let mut o = non_monic_initial_order(&theta);
    let disc = o.discriminant(&theta);
    let disc_fac = factorize(&disc.abs());
    eprintln!("disc={}", disc);
    for &(ref p, e) in &disc_fac {
        eprintln!("disc_fac: {}^{}", p, e);
    }
    for &(ref p, mut e) in &disc_fac {
        while e >= 2 {
            eprintln!("Considering p = {}: {}^{} | disc", p, p, e);
            // I_p
            let pow = p.pow(3u32); // note: in the real implementation, pow should be a power of p.

            // phi(w_i)
            let mut phiw = vec![];
            for i in 0..deg {
                let val = create_num(&o.basis[i], &theta);
                phiw.push(power(&val, &pow));
            }
            let mut lcm = BigInt::one();
            for i in 0..deg {
                for j in 0..deg {
                    let val = phiw[i].expr.coef_at(j).denom().clone();
                    lcm = num::integer::lcm(lcm, val);
                }
            }
            // I_p + pO
            let mut basis = vec![vec![BigInt::from(0); deg]; 2 * deg];
            for i in 0..deg {
                for j in 0..deg {
                    let val = (phiw[i].expr.coef_at(j) * &lcm.clone()).to_integer();
                    basis[i][j] = val;
                }
            }
            for i in 0..deg {
                for j in 0..deg {
                    let val = (&o.basis[i][j] * &lcm).to_integer();
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
                        i_p[i][k] += &hnf.0[i][j] * &(&o.basis[j][k] * &lcm).to_integer();
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
            eprintln!("(new_o : o) = {}", index);
            eprintln!("o = {}", o);
            eprintln!("new_o = {}", new_o);
            if index == BigInt::one() {
                eprintln!("done: p = {}", p);
                break;
            }
            o = new_o;
            while index > BigInt::one() {
                assert_eq!(&index % p, BigInt::zero());
                index /= p;
                e -= 2;
            }
        }
    }
    eprintln!("D(Z_K) = {}", o.discriminant(&theta));
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

// Super tenuki
fn power(a: &Algebraic, e: &BigInt) -> Algebraic {
    let mut prod: Algebraic = Algebraic::from_int(a.min_poly.clone(), 1);
    let mut count = BigInt::zero();
    while &count < e {
        prod = &prod * a;
        count += 1;
    }
    prod
}
