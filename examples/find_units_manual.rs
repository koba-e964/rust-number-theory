#![allow(clippy::needless_range_loop)]

use num::{bigint::Sign, BigInt, BigRational, Complex, One, ToPrimitive, Zero};
use number_theory_linear::determinant_real;
use number_theory_linear::hnf::{self, HNF};
use rand::Rng;
use rust_number_theory::{
    algebraic::Algebraic,
    ideal::Ideal,
    integral_basis::find_integral_basis,
    mult_table::MultTable,
    numerical_roots::find_roots_reim,
    order::{self, Order},
    poly_mod::find_linear_factors,
    polynomial::Polynomial,
};
use std::collections::HashMap;

/// (Prime ideal, its residue class degree)
type PrimeIdeal<'mul> = (Ideal<'mul>, usize);

// Returns the factorization of (p). If this function fails to compute it, this function returns None.
fn factor_prime<'mul>(
    p: &BigInt,
    poly: &Polynomial<BigInt>,
    mult_table: &'mul MultTable,
    o: &Order,
    theta: &Algebraic,
) -> Option<Vec<PrimeIdeal<'mul>>> {
    let pow_basis = Order::singly_gen(theta);
    let index = order::index(o, &pow_basis);
    if (index % p).is_zero() {
        return None;
    }
    let mut factors = find_linear_factors::<BigInt>(poly, p.clone());
    if factors.len() != poly.deg() {
        return None;
    }
    factors.sort();
    factors.dedup();
    let deg = poly.deg();
    let mut ideals = vec![];
    for a in factors {
        // (p, theta - a)
        let mut pnum = vec![BigInt::zero(); deg];
        pnum[0] = p.clone();
        let ideal = Ideal::principal(&pnum, mult_table);
        let theta_a = Algebraic::with_expr(
            theta.min_poly.clone(),
            Polynomial::from_raw(vec![-BigRational::from(a.clone()), BigInt::one().into()]),
        );
        let theta_a = o.to_z_basis_int(&theta_a);
        let ideal = &ideal + &Ideal::principal(&theta_a, mult_table);
        eprintln!(
            "p = {}, a = {}, theta_a = {:?}, ideal = {:?} (norm= {})",
            p,
            a,
            theta_a,
            ideal,
            ideal.norm(),
        );
        ideals.push((ideal, 1));
    }
    Some(ideals)
}

fn factorize_with_known_primes<'mul>(
    num: &[BigInt],
    map: &HashMap<BigInt, Vec<PrimeIdeal<'mul>>>,
    mult_table: &'mul MultTable,
) -> Option<Vec<(BigInt, usize)>> {
    let norm = mult_table.norm(&num);
    if norm.is_zero() {
        return None;
    }
    let mut remaining = norm.clone();
    let mut factors = vec![];
    for (p, ps) in map {
        let mut e = 0;
        while (&remaining % p).is_zero() {
            e += 1;
            remaining /= p;
        }
        let mut dividing = vec![];
        let mut fsum = 0;
        for i in 0..ps.len() {
            let &(ref pideal, f) = &ps[i];
            if pideal.contains(num) {
                dividing.push((pideal.clone(), i));
                fsum += f;
            }
        }
        if fsum == e {
            // Each prime ideal divides num exactly once.
            for (_, idx) in dividing {
                factors.push((p.clone(), idx));
            }
        } else if dividing.len() == 1 {
            // Only one prime ideal on (p) divides num. (num) = (that prime ideal)^e.
            for _ in 0..e / fsum {
                factors.push((p.clone(), dividing[0].1));
            }
        } else {
            return None;
        }
    }
    if remaining.pow(2).is_one() {
        eprintln!("num = {:?}, norm = {}, factors = {:?}", num, norm, factors);
        Some(factors)
    } else {
        None
    }
}

fn main() {
    let mut rng = rand::thread_rng();

    let poly_vec: Vec<BigInt> = vec![(-41).into(), 0.into(), 1.into()];
    let poly = Polynomial::from_raw(poly_vec.clone());
    let poly_complex =
        Polynomial::from_raw(poly_vec.into_iter().map(|b| b.to_f64().unwrap()).collect());
    let deg = poly.deg();
    let theta = Algebraic::new(poly.clone());
    let o = find_integral_basis(&theta);
    eprintln!("o = {:?}", o);
    let mult_table = o.get_mult_table(&theta);
    let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43];
    let mut map = HashMap::new();
    let mut offsets = HashMap::new();
    let mut offset = 0;
    for &p in &primes {
        let p = BigInt::from(p);
        if let Some(ps) = factor_prime(&p, &poly, &mult_table, &o, &theta) {
            map.insert(p.clone(), ps.clone());
            offsets.insert(p.clone(), offset);
            offset += ps.len();
        }
    }
    // Find integers with factorization with small primes
    let w = offset;
    let mut rows = vec![];
    let mut nums = vec![];
    for a in 0..20 {
        for b in -10..10 {
            let num: Vec<BigInt> = vec![a.into(), b.into()];
            if let Some(factors) = factorize_with_known_primes(&num, &map, &mult_table) {
                let mut row = vec![BigInt::zero(); w];
                for (p, idx) in factors {
                    let offset = offsets[&p];
                    row[offset + idx] += 1;
                }
                eprintln!("row = {:?}, num = {:?}", row, num);
                rows.push(row);
                nums.push(num);
            }
        }
    }
    let h = rows.len();
    let ker = HNF::kernel(&rows);
    let (principal, _u, _k) = hnf::hnf_with_u(&rows);
    let cl = principal.determinant();
    eprintln!("tentative Cl(K) = {}", cl);
    for p in &principal.0 {
        eprintln!("{:?}", p);
    }
    let mut unit_cand = vec![];
    for entry in ker {
        // Because inverting an integer is a costly operation, we will invert only once in the last step.
        let mut num = vec![BigInt::zero(); deg];
        num[0] += 1;
        let mut den = num.clone();
        for i in 0..h {
            if entry[i].sign() == Sign::Plus {
                for _ in num::range(BigInt::zero(), entry[i].clone()) {
                    num = mult_table.mul(&num, &nums[i]);
                }
            }
            if entry[i].sign() == Sign::Minus {
                for _ in num::range(BigInt::zero(), -entry[i].clone()) {
                    den = mult_table.mul(&den, &nums[i]);
                }
            }
        }
        let (deninv, denden) = mult_table.inv(&den);
        let mut res = mult_table.mul(&num, &deninv);
        for i in 0..deg {
            assert_eq!(&res[i] % &denden, BigInt::zero());
            res[i] /= &denden;
        }
        eprintln!("res = {:?}", res);
        unit_cand.push(res);
    }
    let (roots_re, roots_im) = find_roots_reim(poly_complex.clone());
    let r = roots_re.len();
    let s = roots_im.len();
    let mut basis = vec![vec![Complex::new(0.0, 0.0); deg]; r + s];
    let mut lnmatrix = vec![];
    for i in 0..r + s {
        let root = if i < r {
            roots_re[i].into()
        } else {
            roots_im[i]
        };
        for j in 0..deg {
            let mut current = Complex::new(1.0, 0.0);
            for k in 0..deg {
                basis[i][j] += o.basis[j][k].to_f64().unwrap() * current;
                current *= root;
            }
        }
        eprintln!("root = {}, basis = {:?}", root, basis[i]);
    }
    for i in 0..unit_cand.len() {
        let num = &unit_cand[i];
        let mut lnvec = vec![];
        for j in 0..r + s {
            let mut val = Complex::new(0.0, 0.0);
            for k in 0..deg {
                val += basis[j][k] * num[k].to_f64().unwrap();
            }
            let ln = val.norm_sqr().ln() / 2.0;
            lnvec.push(ln);
        }
        eprintln!("num = {:?}, ln = {:?}", num, lnvec);
        lnmatrix.push(lnvec);
    }
    // randomly pick r + s - 1 elements
    let mut perm: Vec<usize> = (0..unit_cand.len()).collect();
    for i in 0..unit_cand.len() {
        let idx = rng.gen_range(0..i + 1);
        perm.swap(i, idx);
    }
    let mut matrix = vec![vec![0.0; r + s - 1]; r + s - 1];
    for i in 0..r + s - 1 {
        for j in 0..r + s - 1 {
            matrix[i][j] = lnmatrix[perm[i]][j];
        }
    }
    let reg = determinant_real(&matrix);
    eprintln!(
        "tentative Reg(K) = {}, Cl(K) = {}, prod = {}",
        reg,
        cl,
        reg * cl.to_f64().unwrap(),
    );
}
