use num::BigInt;

use num::Zero;

use crate::{
    algebraic::Algebraic,
    ideal::Ideal,
    mult_table::MultTable,
    order::{self, Order},
};

/// Buchmann-Lenstra's algorithm for prime decomposition.
mod bl;
mod simple;

// Decompose a prime into prime ideals in Z_Q(theta).
pub fn decompose<'mul>(
    theta: &Algebraic,
    int_basis: &Order,
    mult_table: &'mul MultTable,
    p: &BigInt,
) -> Vec<(Ideal<'mul>, usize)> {
    let z_theta = order::trivial_order_monic(theta);
    let index = order::index(int_basis, &z_theta);
    if index % p == BigInt::zero() {
        bl::decompose(theta, int_basis, mult_table, p)
    } else {
        simple::decompose(theta, int_basis, mult_table, p)
    }
}

#[cfg(test)]
mod tests {
    use num::traits::Pow;

    use crate::{integral_basis, polynomial::Polynomial};

    use super::*;

    fn p_ideal<'mul>(p: &BigInt, mult_table: &'mul MultTable) -> Ideal<'mul> {
        let mut pelem = vec![BigInt::from(0); mult_table.deg()];
        pelem[0] = p.clone();
        Ideal::principal(&pelem, mult_table)
    }

    fn is_power_of_p(mut x: BigInt, p: &BigInt) -> bool {
        if x <= BigInt::from(0) {
            return false;
        }
        while &x % p == BigInt::from(0) {
            x /= p;
        }
        x == BigInt::from(1)
    }

    fn product_from_decomposition<'mul>(factors: &[(Ideal<'mul>, usize)]) -> Ideal<'mul> {
        assert!(!factors.is_empty());
        let mut result = factors[0].0.clone();
        for _ in 1..factors[0].1 {
            result = &result * &factors[0].0;
        }
        for (ideal, e) in factors.iter().skip(1) {
            for _ in 0..*e {
                result = &result * ideal;
            }
        }
        result
    }

    #[test]
    fn decompose_dispatches_to_bl_for_index_dividing_prime() {
        // Q(sqrt(5)): Z_K = Z[(1 + sqrt(5)) / 2], so (Z_K : Z[sqrt(5)]) = 2.
        let theta = Algebraic::new(Polynomial::from_raw(vec![(-5).into(), 0.into(), 1.into()]));
        let int_basis = integral_basis::find_integral_basis(&theta);
        let mult_table = int_basis.get_mult_table(&theta);
        let p: BigInt = 2.into();
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        assert!(!result.is_empty());
        let mut total: BigInt = 1.into();
        for (ideal, e) in &result {
            assert!(*e > 0);
            assert!(is_power_of_p(ideal.norm(), &p));
            total *= ideal.norm().pow(*e);
        }
        assert_eq!(total, p.clone().pow(2usize));
        let lhs = product_from_decomposition(&result);
        let rhs = p_ideal(&p, &mult_table);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn decompose_dispatches_to_simple_for_non_dividing_prime() {
        let theta = Algebraic::new(Polynomial::from_raw(vec![(-5).into(), 0.into(), 1.into()]));
        let int_basis = integral_basis::find_integral_basis(&theta);
        let mult_table = int_basis.get_mult_table(&theta);
        let p: BigInt = 3.into();
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        assert_eq!(result.len(), 1);
        let (ideal, e) = &result[0];
        assert_eq!(*e, 1);
        assert_eq!(ideal.norm(), p.clone().pow(theta.deg()));
        assert!(is_power_of_p(ideal.norm(), &p));
        let lhs = product_from_decomposition(&result);
        let rhs = p_ideal(&p, &mult_table);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn decompose_cubic_x3_plus_9x_plus_1_at_3() {
        let theta = Algebraic::new(Polynomial::from_raw(vec![
            1.into(),
            9.into(),
            0.into(),
            1.into(),
        ]));
        let int_basis = integral_basis::find_integral_basis(&theta);
        let mult_table = int_basis.get_mult_table(&theta);
        let p: BigInt = 3.into();
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        assert_eq!(result.len(), 2);
        for (ideal, _) in &result {
            assert_eq!(ideal.norm(), 3.into());
            assert!(is_power_of_p(ideal.norm(), &p));
        }
        let mut exponents = result.into_iter().map(|(_, e)| e).collect::<Vec<_>>();
        exponents.sort();
        assert_eq!(exponents, vec![1, 2]);
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        let lhs = product_from_decomposition(&result);
        let rhs = p_ideal(&p, &mult_table);
        assert_eq!(lhs, rhs);
    }
}
