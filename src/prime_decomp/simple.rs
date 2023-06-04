use std::convert::TryInto;

use num::{BigInt, BigRational, Zero};

use crate::{
    algebraic::Algebraic,
    ideal::Ideal,
    mult_table::MultTable,
    order::{self, Order},
    poly_mod,
    polynomial::Polynomial,
};

pub fn decompose<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
) -> Vec<(Ideal<'mul>, usize)> {
    let z_theta = order::trivial_order_monic(theta);
    let index = order::index(int_basis, &z_theta);
    if index % p == BigInt::zero() {
        panic!("not (p | (Z_K : Z[theta])) must hold");
    }
    let result = poly_mod::factorize_mod_p::<BigInt>(&theta.min_poly, p, p.try_into().unwrap_or(0));
    result
        .into_iter()
        .map(|(poly, mul)| {
            let poly = Polynomial::from_raw(
                poly.into_vec()
                    .into_iter()
                    .map(BigRational::from_integer)
                    .collect(),
            );
            let elem = if poly.deg() >= theta.min_poly.deg() {
                vec![BigInt::from(0); theta.deg()]
            } else {
                int_basis.to_z_basis_int(&Algebraic::with_expr(theta.min_poly.clone(), poly))
            };
            let ancilla = Ideal::principal(&elem, mult_table);
            let mut pelem = vec![BigInt::from(0); theta.deg()];
            pelem[0] = p.clone();
            let pz = Ideal::principal(&pelem, mult_table);
            (&ancilla + &pz, mul)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::integral_basis::find_integral_basis;

    use super::*;

    #[test]
    fn decompose_works_0() {
        let p = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        let mult_table = order.get_mult_table(&theta);
        let result = decompose(&theta, &order, &mult_table, &3.into());
        // (3) is a prime ideal in Z[i]
        assert_eq!(result.len(), 1);
        let (pid, mul) = result[0].clone();
        assert_eq!(pid.norm(), 9.into());
        assert_eq!(mul, 1);
    }

    #[test]
    fn decompose_works_1() {
        let p = Polynomial::from_raw(vec![1.into(), 0.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        let mult_table = order.get_mult_table(&theta);
        let result = decompose(&theta, &order, &mult_table, &5.into());
        // (5) = P1 P2 in Z[i] where P1 != P2
        assert_eq!(result.len(), 2);
        for (pid, mul) in result {
            assert_eq!(pid.norm(), 5.into());
            assert_eq!(mul, 1);
        }
    }

    #[test]
    fn decompose_works_2() {
        let p = Polynomial::from_raw(vec![2.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        let mult_table = order.get_mult_table(&theta);
        let result = decompose(&theta, &order, &mult_table, &5.into());
        // Every prime gives a prime ideal in Z
        assert_eq!(result.len(), 1);
        let (pid, mul) = result[0].clone();
        assert_eq!(pid.norm(), 5.into());
        assert_eq!(mul, 1);
    }
}
