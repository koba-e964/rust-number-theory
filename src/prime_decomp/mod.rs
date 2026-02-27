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
    use num::{traits::Pow, One};

    use crate::{integral_basis, polynomial::Polynomial};

    use super::*;

    #[test]
    fn decompose_dispatches_to_bl_for_index_dividing_prime() {
        // Q(sqrt(5)): Z_K = Z[(1 + sqrt(5)) / 2], so (Z_K : Z[sqrt(5)]) = 2.
        let theta = Algebraic::new(Polynomial::from_raw(vec![(-5).into(), 0.into(), 1.into()]));
        let int_basis = integral_basis::find_integral_basis(&theta);
        let mult_table = int_basis.get_mult_table(&theta);
        let p: BigInt = 2.into();
        let result = decompose(&theta, &int_basis, &mult_table, &p);
        assert!(!result.is_empty());
        let mut total = BigInt::one();
        for (ideal, e) in result {
            assert!(e > 0);
            assert!(ideal.norm() >= BigInt::one());
            total *= ideal.norm().pow(e);
        }
        assert!(total >= BigInt::one());
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
        assert_eq!(ideal.norm(), p.pow(theta.deg()));
    }
}
