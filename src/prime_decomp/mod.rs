use num::BigInt;

use crate::{algebraic::Algebraic, ideal::Ideal, mult_table::MultTable, order::Order};

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
    // TODO: support if (p | (Z_K : Z[theta]))
    simple::decompose(theta, int_basis, mult_table, p)
}
