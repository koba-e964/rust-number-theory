#![allow(clippy::needless_range_loop, clippy::many_single_char_names)]

use num::{bigint::Sign, BigInt, BigRational, Integer, One, ToPrimitive, Zero};
use number_theory_elementary::primes;
use number_theory_linear::determinant_real;
use number_theory_linear::hnf::{self, HNF};
use rand::Rng;
use rust_number_theory::{
    algebraic::Algebraic,
    class::roots_of_unity::find_muk,
    embeddings::CEmbeddings,
    ideal::Ideal,
    integral_basis::find_integral_basis,
    mult_table::MultTable,
    numerical_roots::find_roots_reim,
    order::{self, Order},
    poly_mod::find_linear_factors,
    polynomial::Polynomial,
};
use std::collections::{HashMap, HashSet};

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
    let deg = poly.deg();
    let pow_basis = Order::singly_gen(theta);
    let index = order::index(o, &pow_basis);
    if (index % p).is_zero() {
        // Ad-hoc factorization: if poly.deg() == 2 and poly = x^2 - a for a mod 4  == 1, we know (2) splits or is inert, depending on a mod 8.
        if p == &BigInt::from(2)
            && poly.coef_at(1) == 0.into()
            && poly.coef_at(2) == 1.into()
            && poly.coef_at(0).mod_floor(&BigInt::from(4)) == 3.into()
        {
            let a = -poly.coef_at(0);
            let rem = a.mod_floor(&BigInt::from(8));
            // (2)
            let mut pnum = vec![BigInt::zero(); deg];
            pnum[0] = BigInt::from(2);
            let ideal = Ideal::principal(&pnum, mult_table);
            if rem == 1.into() {
                // (2) splits. (2) = (2, eta)(2, eta+1) where eta = (1+theta) / 2
                let eta = Algebraic::with_expr(
                    theta.min_poly.clone(),
                    Polynomial::from_raw(vec![
                        BigRational::new(BigInt::one(), BigInt::from(2)),
                        BigRational::new(BigInt::one(), BigInt::from(2)),
                    ]),
                );
                let mut eta = o.to_z_basis_int(&eta);
                let ideal_0 = &ideal + &Ideal::principal(&eta, mult_table);
                eta[0] += 1;
                let ideal_1 = &ideal + &Ideal::principal(&eta, mult_table);
                return Some(vec![(ideal_0, 1), (ideal_1, 1)]);
            } else {
                // (2) is inert
                return Some(vec![(ideal, 2)]);
            }
        }
        return None;
    }
    let mut factors = find_linear_factors::<BigInt>(poly, p.clone());
    // ad-hoc factorization: if poly.deg() == 2 and there are no linear factors, we know that poly is irreducible.
    if factors.is_empty() && poly.deg() == 2 {
        // (p)
        let mut pnum = vec![BigInt::zero(); deg];
        pnum[0] = p.clone();
        let ideal = Ideal::principal(&pnum, mult_table);
        return Some(vec![(ideal, 2)]);
    }
    if factors.len() != poly.deg() {
        return None;
    }
    factors.sort();
    factors.dedup();
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
    let mut remaining = norm;
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
        if e == 0 {
            assert_eq!(
                fsum, 0,
                "fsum = {}, e = {}, rem = {}, dividing = {:?}",
                fsum, e, remaining, dividing
            );
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
        Some(factors)
    } else {
        None
    }
}

fn euler_prod<'mul>(primes: &[i32], map: &HashMap<BigInt, Vec<PrimeIdeal<'mul>>>) -> f64 {
    let mut ans = 1.0;
    for &p in primes {
        let big_p = BigInt::from(p);
        if !map.contains_key(&big_p) {
            eprintln!("error!! p = {}", p);
        }
        let on_p = &map[&big_p];
        let pinv = 1.0 / p as f64;
        for &(_, f) in on_p {
            let mut tmp = 1.0;
            for _ in 0..f {
                tmp *= pinv;
            }
            ans *= 1.0 - tmp;
        }
        ans /= 1.0 - pinv;
    }
    1.0 / ans
}

// TODOs:
// Use Minkowski bounds to enumerate primes
// Incrementally enumerate relations
// Factorize all primes
fn main() {
    let mut rng = rand::thread_rng();

    let poly_vec: Vec<BigInt> = vec![(-1141).into(), 1.into(), 1.into()];
    let poly = Polynomial::from_raw(poly_vec.clone());
    let poly_complex =
        Polynomial::from_raw(poly_vec.into_iter().map(|b| b.to_f64().unwrap()).collect());
    let deg = poly.deg();
    let theta = Algebraic::new(poly.clone());
    let o = find_integral_basis(&theta);
    eprintln!("o = {:?}", o);

    // Find a suitable bound
    let disc = o.discriminant(&theta);
    let disc_ln = disc.to_f64().unwrap().abs().ln();
    let coef = 1.0;
    let bound = coef * disc_ln * disc_ln;
    eprintln!("bound = {}", bound);

    // Find embeddings and roots of unity
    let (roots_re, roots_im) = find_roots_reim(poly_complex);
    let r = roots_re.len();
    let s = roots_im.len();
    let basis = CEmbeddings::new(&roots_re, &roots_im, &o);
    let muk = find_muk(&basis);

    let mult_table = o.get_mult_table(&theta);
    let primes: Vec<i32> = primes(bound.floor() as usize)
        .into_iter()
        .map(|x| x as i32)
        .collect();
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
    // First process rational primes so that every prime appears at least once
    for &p in &primes {
        let num: Vec<BigInt> = vec![p.into(), 0.into()];
        if let Some(factors) = factorize_with_known_primes(&num, &map, &mult_table) {
            eprintln!("prime p = {}", p);
            let mut row = vec![BigInt::zero(); w];
            for (p, idx) in factors {
                let offset = offsets[&p];
                row[offset + idx] += 1;
            }
            rows.push(row);
            nums.push(num);
        }
    }
    for a in 0..30 {
        for b in -10..10 {
            if b == 0 {
                continue;
            }
            let num: Vec<BigInt> = vec![a.into(), b.into()];
            if let Some(factors) = factorize_with_known_primes(&num, &map, &mult_table) {
                let mut row = vec![BigInt::zero(); w];
                for (p, idx) in factors {
                    let offset = offsets[&p];
                    row[offset + idx] += 1;
                }
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
    let mut unseen: HashSet<usize> = (0..w).collect();
    for p in &principal.0 {
        for i in 0..w {
            if !p[i].is_zero() {
                unseen.remove(&i);
            }
        }
    }
    if cl.is_zero() {
        for idx in unseen {
            let mut pid = None;
            for (p, &o) in &offsets {
                if o <= idx && idx < map[p].len() + o {
                    pid = Some(map[p][idx - o].clone());
                }
            }
            eprintln!("idx = {}, p = {:?}", idx, pid);
        }
        return;
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
        unit_cand.push(res);
    }
    let mut lnmatrix = vec![];
    for i in 0..unit_cand.len() {
        let num = &unit_cand[i];
        let mut lnvec = vec![];
        for j in 0..r + s {
            let val = basis.compute(j, num);
            let ln = val.norm_sqr().ln() / 2.0;
            lnvec.push(ln);
        }
        lnmatrix.push(lnvec);
    }
    loop {
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
        let reg = determinant_real(&matrix).abs();
        let product = reg * cl.to_f64().unwrap();
        eprintln!(
            "tentative Reg(K) = {}, Cl(K) = {}, prod = {}",
            reg, cl, product,
        );
        let euler_prod = euler_prod(&primes, &map);
        // https://www.isibang.ac.in/~sury/algoiisc.pdf
        let mut a = euler_prod;
        a *= muk as f64;
        a *= o.discriminant(&theta).to_f64().unwrap().abs().sqrt();
        a /= 2.0f64.powf(r as f64);
        a /= (2.0 * std::f64::consts::PI).powf(s as f64);
        eprintln!("a = {}", a);
        if product > 0.707 * a && product < 1.414 * a {
            eprintln!("We found the correct unit group and the class group. Stopping.");
            eprintln!("generators:");
            for i in 0..r + s - 1 {
                eprintln!("{:?}", unit_cand[perm[i]]);
            }
            break;
        }
    }
}
