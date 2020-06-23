#![allow(clippy::needless_range_loop)]

extern crate num;
extern crate rust_number_theory;

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
    if v > BigInt::one() {
        ans.push((x, 1));
    }
    ans
}

#[allow(unused)]
fn main_sub() {
    let p = Polynomial::from_raw(vec![
        5.into(),
        6.into(),
        (-7).into(),
        6.into(),
        (-7).into(),
        6.into(),
    ]);
    let deg = p.deg();
    let theta = Algebraic::new(p);
    let mut o = non_monic_initial_order(&theta);
    let disc = o.discriminant(&theta);
    let disc_fac = factorize(&disc);
    eprintln!("disc={}", disc);
    for &(ref p, e) in &disc_fac {
        eprintln!("disc_fac: {}^{}", p, e);
    }
    #[allow(unused_assignments)]
    #[allow(clippy::never_loop)]
    for &(ref p, mut e) in &disc_fac {
        while e >= 2 {
            eprintln!("Considering p = {}: {}^{} | disc", p, p, e);
            e -= 2;
            // I_p
            let pow = 5.into(); // note: in the real implementation, pow should be a power of p.
                                // phi(w_i)
            let mut phiw = vec![];
            for i in 0..deg {
                let val = create_num(&o.basis[i], &theta);
                phiw.push(power(&val, &pow));
            }
            eprintln!("phiw = {:?}", phiw);
            let mut lcm = BigInt::one();
            for i in 0..deg {
                for j in 0..deg {
                    let val = phiw[i].expr.coef_at(j).denom().clone();
                    lcm = num::integer::lcm(lcm, val);
                }
            }
            eprintln!("lcm = {}", lcm);
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
            eprintln!("basis = {}", HNF(basis.clone()));
            let mut hnf = HNF::hnf(&HNF::kernel(&basis));
            for i in 0..deg {
                for j in 0..deg {}
                hnf.0[i].truncate(deg);
            }
            // Hack: In this case, I_p = U_p.
            eprintln!("hnf (I_p) = {}", hnf);
            let mut o_basis = vec![vec![BigInt::zero(); deg]; deg];
            for i in 0..deg {
                for j in 0..deg {
                    let val = (&o.basis[i][j] * &lcm).to_integer();
                    o_basis[i][j] = &val * p;
                }
            }
            // inner product

            let o_basis = HNF::union(&HNF(o_basis), &hnf);
            eprintln!("o_basis = {}", o_basis);
            let mut new_o_basis = vec![vec![BigRational::zero(); deg]; deg];
            for i in 0..deg {
                for j in 0..deg {
                    new_o_basis[i][j] = BigRational::new(o_basis.0[i][j].clone(), &lcm * p);
                }
            }
            let new_o = Order { basis: new_o_basis };
            eprintln!("(new_o : o) = {}", index(&new_o, &o));
            eprintln!("o = {}", o);
            eprintln!("new_o = {}", new_o);
            o = new_o;
            return;
        }
    }
}

fn main() {
    let p = Polynomial::from_raw(vec![1.into(), 9.into(), 0.into(), 1.into()]);
    let deg = p.deg();
    let theta = Algebraic::new(p);
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
            e -= 2;
            // I_p
            let pow = 3.into(); // note: in the real implementation, pow should be a power of p.

            // phi(w_i)
            let mut phiw = vec![];
            for i in 0..deg {
                let val = create_num(&o.basis[i], &theta);
                phiw.push(power(&val, &pow));
            }
            eprintln!("phiw = {:?}", phiw);
            let mut lcm = BigInt::one();
            for i in 0..deg {
                for j in 0..deg {
                    let val = phiw[i].expr.coef_at(j).denom().clone();
                    lcm = num::integer::lcm(lcm, val);
                }
            }
            eprintln!("lcm = {}", lcm);
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
            eprintln!("basis = {}", HNF(basis.clone()));
            let mut hnf = HNF::hnf(&HNF::kernel(&basis));
            for i in 0..deg {
                hnf.0[i].truncate(deg);
            }
            eprintln!("hnf (I_p) = {}", hnf);
            let i_p_len = hnf.0.len();

            // U_p
            let mut u_p = hnf.clone();

            for i in 0..i_p_len {
                let ei = create_num_int(&hnf.0[i], &theta);
                let mut tmp_basis = vec![vec![BigInt::zero(); deg]; u_p.0.len() + hnf.0.len()];
                for i in 0..u_p.0.len() {
                    let prod = &create_num_int(&u_p.0[i], &theta) * &ei;
                    for j in 0..deg {
                        tmp_basis[i][j] = (&prod.expr.coef_at(j) / &lcm).to_integer();
                    }
                }
                for i in 0..i_p_len {
                    for j in 0..deg {
                        tmp_basis[i + u_p.0.len()][j] = &hnf.0[i][j] * p;
                    }
                }
                // new_u_p is in terms of U_p + pI_p
                let mut new_u_p = HNF::hnf(&HNF::kernel(&tmp_basis));
                for i in 0..deg {
                    new_u_p.0[i].truncate(deg);
                }
                // In terms of theta
                let mut tmp_basis = vec![vec![BigInt::zero(); i_p_len]; new_u_p.0.len()];
                for i in 0..new_u_p.0.len() {
                    for j in 0..u_p.0.len() {
                        for k in 0..deg {
                            tmp_basis[i][k] += &u_p.0[j][k] * &new_u_p.0[i][j];
                        }
                    }
                }
                let new_u_p = HNF::hnf(&tmp_basis);
                eprintln!("new U_p = {}", new_u_p);
                u_p = new_u_p;
            }

            eprintln!("hnf (U_p) = {}", u_p);

            eprintln!("U_p = {}", u_p);
            assert_eq!(u_p.0.len(), deg);
            let mut new_o_basis = vec![vec![BigRational::zero(); deg]; deg];
            for i in 0..deg {
                for j in 0..deg {
                    new_o_basis[i][j] = BigRational::new(u_p.0[i][j].clone(), &lcm * p);
                }
            }
            let new_o = Order { basis: new_o_basis };
            eprintln!("(new_o : o) = {}", index(&new_o, &o));
            eprintln!("o = {}", o);
            eprintln!("new_o = {}", new_o);
            o = new_o;
        }
    }
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
