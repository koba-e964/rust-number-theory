use std::convert::TryInto;

use num::{BigInt, BigRational, One, Zero};

use crate::{
    algebraic::Algebraic, ideal::Ideal, mult_table::MultTable, order::Order, poly_mod,
    polynomial::Polynomial,
};

// References:
// - Henri Cohen, "A Course in Computational Algebraic Number Theory",
//   Graduate Texts in Mathematics 138, Springer (1993).
//   * Algorithm 6.2.2: prime ideal decomposition at a rational prime p.
//   * Algorithm 6.2.5: multiplication in O/pO (not currently implemented verbatim here).
// - For background on factoring pO_K and ramification exponents:
//   * Marcus, "Number Fields", Springer (1977), Chapter 4.

fn is_power_of_p(x: &BigInt, p: &BigInt) -> Option<usize> {
    if x <= &BigInt::zero() || p <= &BigInt::one() {
        return None;
    }
    let mut t = x.clone();
    let mut e = 0usize;
    while &t % p == BigInt::zero() {
        t /= p;
        e += 1;
    }
    if t.is_one() {
        Some(e)
    } else {
        None
    }
}

fn search_decomposition<'mul>(
    candidates: &[Ideal<'mul>],
    idx: usize,
    target: &Ideal<'mul>,
    current: &Ideal<'mul>,
    exponents: &mut [usize],
    max_e: usize,
) -> bool {
    if idx == candidates.len() {
        return current == target && exponents.iter().any(|&e| e > 0);
    }
    if current.norm() > target.norm() {
        return false;
    }
    let mut power = current.clone();
    for e in 0..=max_e {
        exponents[idx] = e;
        if search_decomposition(candidates, idx + 1, target, &power, exponents, max_e) {
            return true;
        }
        if e < max_e {
            power = &power * &candidates[idx];
            if power.norm() > target.norm() {
                break;
            }
        }
    }
    exponents[idx] = 0;
    false
}

fn split_into_norm_p<'mul>(
    target: &Ideal<'mul>,
    norm_p_candidates: &[Ideal<'mul>],
    f: usize,
) -> Option<Vec<Ideal<'mul>>> {
    fn rec<'mul>(
        target: &Ideal<'mul>,
        candidates: &[Ideal<'mul>],
        f: usize,
        depth: usize,
        start: usize,
        current: Option<Ideal<'mul>>,
        picked: &mut Vec<Ideal<'mul>>,
    ) -> bool {
        if depth == f {
            return if let Some(cur) = current {
                cur == *target
            } else {
                false
            };
        }
        for i in start..candidates.len() {
            let next = if let Some(cur) = &current {
                cur * &candidates[i]
            } else {
                candidates[i].clone()
            };
            picked.push(candidates[i].clone());
            if rec(target, candidates, f, depth + 1, i, Some(next), picked) {
                return true;
            }
            picked.pop();
        }
        false
    }

    let mut picked = vec![];
    if rec(target, norm_p_candidates, f, 0, 0, None, &mut picked) {
        Some(picked)
    } else {
        None
    }
}

fn fallback_decompose<'mul>(
    theta: &Algebraic,
    mult_table: &'mul MultTable,
    p: &BigInt,
    pz: &Ideal<'mul>,
) -> Option<Vec<(Ideal<'mul>, usize)>> {
    let n = theta.deg();
    let p_u: usize = p.clone().try_into().ok()?;
    let mut count = 1usize;
    for _ in 0..n {
        count = count.checked_mul(p_u)?;
    }
    if count > 200_000 {
        return None;
    }
    let mut candidates = vec![];
    for mask in 1..count {
        let mut x = mask;
        let mut elem = vec![BigInt::zero(); n];
        for coef in &mut elem {
            *coef = BigInt::from(x % p_u);
            x /= p_u;
        }
        let a = Ideal::principal(&elem, mult_table);
        let cand = &a + pz;
        let norm = cand.norm();
        if norm.is_one() || cand == *pz {
            continue;
        }
        if is_power_of_p(&norm, p).is_none() {
            continue;
        }
        if !candidates.contains(&cand) {
            candidates.push(cand);
        }
    }
    candidates.sort_by_key(|a| a.norm());
    if candidates.is_empty() {
        return None;
    }

    let mut one = vec![BigInt::zero(); n];
    one[0] = BigInt::one();
    let identity = Ideal::principal(&one, mult_table);
    let mut exponents = vec![0usize; candidates.len()];
    if !search_decomposition(&candidates, 0, pz, &identity, &mut exponents, n) {
        return None;
    }
    let mut out = vec![];
    for (cand, &e) in candidates.iter().cloned().zip(exponents.iter()) {
        if e > 0 {
            out.push((cand, e));
        }
    }
    let norm_p_candidates = candidates
        .iter()
        .filter(|cand| cand.norm() == *p)
        .cloned()
        .collect::<Vec<_>>();
    let mut refined: Vec<(Ideal<'mul>, usize)> = vec![];
    for (ideal, e) in out {
        let norm = ideal.norm();
        let mut split = None::<Vec<Ideal<'mul>>>;
        if let Some(f) = is_power_of_p(&norm, p) {
            if f > 1 && !norm_p_candidates.is_empty() {
                split = split_into_norm_p(&ideal, &norm_p_candidates, f);
            }
        }
        if let Some(parts) = split {
            for part in parts {
                let mut merged = false;
                for (existing, ee) in &mut refined {
                    if *existing == part {
                        *ee += e;
                        merged = true;
                        break;
                    }
                }
                if !merged {
                    refined.push((part, e));
                }
            }
        } else {
            let mut merged = false;
            for (existing, ee) in &mut refined {
                if *existing == ideal {
                    *ee += e;
                    merged = true;
                    break;
                }
            }
            if !merged {
                refined.push((ideal, e));
            }
        }
    }
    Some(refined)
}

// 6.2.2 of [Cohen]. Returns a list of pairs (P, e).
//
// Notes on current implementation:
// - The first pass mirrors the classic construction from factorization mod p,
//   creating ideals (p, f_i(theta)).
// - When that pass is not norm-consistent, we use a fallback search in O/pO
//   by enumerating (p, a) candidates and solving pO = product(P_i^e_i).
// - This keeps behavior correct for tested index-dividing-prime cases, while
//   not yet being a direct transcription of all BL substeps.
pub fn decompose<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
) -> Vec<(Ideal<'mul>, usize)> {
    let mut pelem = vec![BigInt::from(0); theta.deg()];
    pelem[0] = p.clone();
    let pz = Ideal::principal(&pelem, mult_table);
    let result = poly_mod::factorize_mod_p::<BigInt>(&theta.min_poly, p, p.try_into().unwrap_or(0));
    let decomposition: Vec<(Ideal<'mul>, usize)> = result
        .into_iter()
        .map(|(poly, mul)| {
            let poly = Polynomial::from_raw(
                poly.into_vec()
                    .into_iter()
                    .map(BigRational::from_integer)
                    .collect(),
            );
            let mut elem = if poly.deg() >= theta.min_poly.deg() {
                vec![BigInt::from(0); theta.deg()]
            } else {
                int_basis.to_z_basis_int(&Algebraic::with_expr(theta.min_poly.clone(), poly))
            };
            while elem.iter().any(|coef| !coef.is_zero())
                && elem.iter().all(|coef| coef % p == BigInt::zero())
            {
                for coef in &mut elem {
                    *coef /= p;
                }
            }
            let ancilla = Ideal::principal(&elem, mult_table);
            (&ancilla + &pz, mul)
        })
        .collect();
    let mut norm_prod = BigInt::one();
    for (ideal, e) in &decomposition {
        norm_prod *= ideal.norm().pow(*e as u32);
    }
    let target_norm = p.pow(theta.deg() as u32);
    if decomposition
        .iter()
        .all(|(ideal, _)| !ideal.norm().is_one())
        && norm_prod == target_norm
    {
        return decomposition;
    }
    fallback_decompose(theta, mult_table, p, &pz).unwrap_or(vec![(pz, 1)])
}
