pub fn primes(bound: usize) -> Vec<usize> {
    let mut is_prime = vec![true; bound + 1];
    is_prime[0] = false;
    if bound >= 1 {
        is_prime[1] = false;
    }
    for i in 2..=bound {
        if !is_prime[i] {
            continue;
        }
        for j in 2..=bound / i {
            is_prime[i * j] = false;
        }
    }
    (2..=bound).filter(|&i| is_prime[i]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primes_works_1() {
        let bound = 50;
        let expected = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        assert_eq!(primes(bound), expected);
    }

    #[test]
    fn primes_works_2() {
        let bound = 53;
        let expected = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
        assert_eq!(primes(bound), expected);
    }

    #[test]
    fn primes_works_3() {
        let bound = 49;
        let expected = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        assert_eq!(primes(bound), expected);
    }
}
