use core::ops::Neg;

use num::{
    traits::{NumAssign, NumOps},
    Integer, Signed,
};
use rand::distributions::uniform::SampleUniform;

use crate::polynomial::Polynomial;

pub fn factorize<
    Int: Clone + Integer + NumAssign + Signed + Neg<Output = Int> + From<i32> + SampleUniform,
>(
    a: &Polynomial<Int>,
) -> Vec<(Polynomial<Int>, usize)>
where
    for<'a> &'a Int: NumOps<&'a Int, Int>,
{
    // Theorem 3.5.1 in [Cohen]
    let n = a.deg();
    let mut sum = a.coef_at(n).abs();
    for i in 0..n + 1 {
        sum += a.coef_at(i).abs();
    }
    let two = Int::one() + Int::one();
    let mut bound = sum;
    for _ in 0..n - 1 {
        bound = &bound * &two;
    }
    bound = &(&bound * &two) * &a.coef_at(n);

    let p = two; // TODO: find optimal primes
    let mut e = 0;
    {
        let mut pe = Int::one();
        while pe <= bound {
            pe = &pe * &p;
            e += 1;
        }
    }
    eprintln!("e = {}", e);
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorize_works_0() {
        let a = Polynomial::<i32>::from_raw(vec![3, 2, 1]);
        let result = factorize::<i32>(&a);
        eprintln!("result = {:?}", result);
    }
}
