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

fn is_prime(a: usize) -> bool {
    if a <= 1 {
        return false;
    }
    let mut d = 2;
    while d * d <= a {
        if a % d == 0 {
            return false;
        }
        d += 1;
    }
    true
}

/// An iterator that returns primes in the increasing order.
///
/// Complexity: O(n^1.5)
pub struct Primes {
    now: usize,
}

impl Primes {
    pub fn new() -> Self {
        Self { now: 2 }
    }
}

impl Iterator for Primes {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        let mut now = self.now;
        while !is_prime(now) {
            now += 1;
        }
        self.now = now + 1;
        Some(now)
    }
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

    #[test]
    fn primes_iterator_works_0() {
        let iter = Primes::new();
        let primes: Vec<usize> = iter.take(15).collect();
        let expected = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        assert_eq!(primes, expected);
    }
}
