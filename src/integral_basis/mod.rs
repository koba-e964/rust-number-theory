use num::Signed;

mod round2;

use crate::algebraic::Algebraic;
use crate::factorize;
use crate::order::{non_monic_initial_order, Order};

/// Finds an integral basis of Q(theta).
pub fn find_integral_basis(theta: &Algebraic) -> Order {
    let mut o = non_monic_initial_order(theta).hnf_reduce();
    let disc = o.discriminant(theta);
    let disc_fac = factorize::factorize(&disc.abs());
    for &(ref p, mut e) in &disc_fac {
        while e >= 2 {
            let (new_o, howmany) = round2::one_step(theta, &o, p);
            e -= 2 * howmany;
            o = new_o;
            if howmany == 0 {
                break;
            }
        }
    }
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::polynomial::Polynomial;

    #[test]
    fn find_integral_basis_works_0() {
        let p = Polynomial::from_raw(vec![
            5.into(),
            6.into(),
            (-7).into(),
            6.into(),
            (-7).into(),
            6.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), 7601837.into());
    }

    #[test]
    fn find_integral_basis_works_1() {
        // theta = -1+6i
        let p = Polynomial::from_raw(vec![37.into(), 2.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        // Z_{Q(theta)} = Z[i]
        assert_eq!(order.discriminant(&theta), (-4).into());
    }

    #[test]
    fn find_integral_basis_works_seq3() {
        let p = Polynomial::from_raw(vec![4.into(), 3.into(), 2.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), (-200).into());
    }

    #[test]
    fn find_integral_basis_works_seq4() {
        let p = Polynomial::from_raw(vec![5.into(), 4.into(), 3.into(), 2.into(), 1.into()]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), 10800.into());
    }

    #[test]
    fn find_integral_basis_works_seq5() {
        let p = Polynomial::from_raw(vec![
            6.into(),
            5.into(),
            4.into(),
            3.into(),
            2.into(),
            1.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), 1037232.into());
    }

    #[test]
    fn find_integral_basis_works_seq6() {
        let p = Polynomial::from_raw(vec![
            7.into(),
            6.into(),
            5.into(),
            4.into(),
            3.into(),
            2.into(),
            1.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), (-9834496).into());
    }

    #[test]
    fn find_integral_basis_works_seq7() {
        let p = Polynomial::from_raw(vec![
            8.into(),
            7.into(),
            6.into(),
            5.into(),
            4.into(),
            3.into(),
            2.into(),
            1.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), (-241864704).into());
    }

    // Too slow. Takes around 6 sec in debug build.
    #[allow(unused)]
    fn find_integral_basis_works_seq8() {
        let p = Polynomial::from_raw(vec![
            9.into(),
            8.into(),
            7.into(),
            6.into(),
            5.into(),
            4.into(),
            3.into(),
            2.into(),
            1.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        eprintln!("disc={}", order.discriminant(&theta));
        assert_eq!(order.discriminant(&theta), 1180980000000u64.into());
    }

    // Too slow. Needs optimization. Takes around 3s in *release* build.
    #[allow(unused)]
    fn find_integral_basis_works_seq10() {
        let p = Polynomial::from_raw(vec![
            11.into(),
            10.into(),
            9.into(),
            8.into(),
            7.into(),
            6.into(),
            5.into(),
            4.into(),
            3.into(),
            2.into(),
            1.into(),
        ]);
        let theta = Algebraic::new(p);
        let order = find_integral_basis(&theta);
        assert_eq!(order.discriminant(&theta), (-8640974550472704i64).into());
    }
}
