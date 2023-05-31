/// Computes a Kronecker-Jacobi symbol.
///
/// This algorithm is Algorithm 1.4.10 in [Cohen].
pub fn kronecker_symbol_i64(mut a: i64, mut b: i64) -> i32 {
    if b == 0 {
        return if a == 1 || a == -1 { 1 } else { 0 };
    }
    if ((a | b) & 1) == 0 {
        return 0;
    }
    let mut k = 1;
    let recip_table = [0, 1, 0, -1, 0, -1, 0, 1];
    // 2.
    {
        let mut v = 0;
        while (b & 1) == 0 {
            v += 1;
            b /= 2;
        }
        if v % 2 == 1 {
            k = recip_table[(a & 7) as usize];
        }
        if b < 0 {
            b = -b;
        }
        if a < 0 {
            k = -k;
        }
    }
    // 3. and 4.
    while a != 0 {
        let mut v = 0;
        while (a & 1) == 0 {
            v += 1;
            a /= 2;
        }
        if v % 2 == 1 {
            k *= recip_table[(b & 7) as usize];
        }
        // 4.
        if (a & b & 2) != 0 {
            k = -k;
        }
        let r = a.abs();
        a = b % r;
        b = r;
    }
    if b == 1 {
        k
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kronecker_symbol_i64_works_0() {
        let n = 13;
        let v: Vec<i32> = (0..n + 1).map(|x| kronecker_symbol_i64(x, n)).collect();
        assert_eq!(v, [0, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, 0]);
    }

    #[test]
    fn kronecker_symbol_i64_works_1() {
        let n = 17;
        let v: Vec<i32> = (0..n + 1).map(|x| kronecker_symbol_i64(x, n)).collect();
        assert_eq!(
            v,
            [0, 1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 0],
        );
    }
}
