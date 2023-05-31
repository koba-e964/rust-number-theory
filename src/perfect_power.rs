use num::{bigint::Sign, BigInt, One};

pub fn perfect_power(n: &BigInt) -> (BigInt, u32) {
    if n.sign() == Sign::Minus {
        panic!("n should be >= 0, but got: {}", n);
    }
    if *n <= BigInt::one() {
        return (n.clone(), 1);
    }
    let numbits = n.bits() as u32;
    for k in (2..=numbits).rev() {
        if let Some(b) = is_perfect_power(n, k) {
            return (b, k);
        }
    }
    (n.clone(), 1)
}

pub fn is_perfect_power(n: &BigInt, k: u32) -> Option<BigInt> {
    let x = n.nth_root(k);
    if x.pow(k) == *n {
        Some(x)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_power_works_0() {
        assert_eq!(perfect_power(&27.into()), (3.into(), 3));
        assert_eq!(perfect_power(&64.into()), (2.into(), 6)); // not 4^3 or 8^2
        assert_eq!(perfect_power(&26.into()), (26.into(), 1));
    }
}
