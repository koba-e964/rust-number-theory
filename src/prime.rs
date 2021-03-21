use num::bigint::RandBigInt;
use num::{BigInt, One};

#[allow(clippy::many_single_char_names)]
pub fn is_prime(n: &BigInt) -> bool {
    if n <= &BigInt::one() {
        return false;
    }
    let two = BigInt::from(2u8);
    if n == &two {
        return true;
    }
    // if n % 2 == 0
    if !n.bit(0) {
        return false;
    }

    let mut c = 0u64;
    let mut d: BigInt = n - 1;

    while !d.bit(0) {
        d >>= 1;
        c += 1;
    }

    let mut rng = rand::thread_rng();

    let k = 20;
    for _ in 0..k {
        let r = rng.gen_bigint_range(&1.into(), n);
        let mut tmp = r.modpow(&d, n);
        if tmp == BigInt::one() {
            continue;
        }
        let mut aborted = false;
        for _ in 0..c {
            if tmp == n - 1 {
                aborted = true;
                break;
            }
            tmp = tmp.modpow(&2.into(), n);
            if tmp == BigInt::one() {
                // witness found
                return false;
            }
        }
        if !aborted && tmp != BigInt::one() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_prime_works_0() {
        let expected = [
            false, false, true, true, false, true, false, true, false, false,
        ];
        #[allow(clippy::needless_range_loop)]
        for i in 0..10 {
            assert_eq!(is_prime(&BigInt::from(i)), expected[i], "i = {}", i);
        }
    }

    #[test]
    fn is_prime_works_1() {
        // These two integers are known to be primes.
        let large1 = BigInt::from(50920368630756370485097u128);
        let large2 = BigInt::from(1106020951959326672306393u128);
        assert!(is_prime(&large1));
        assert!(is_prime(&large2));
        assert!(!is_prime(&(large1 * large2)));
    }
}
