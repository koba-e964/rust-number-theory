pub fn determinant_real(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let mut a = a.to_vec();
    let mut result = 1.0;
    for i in 0..n {
        let mut idx = None;
        #[allow(clippy::needless_range_loop)]
        for j in i..n {
            if a[j][i] != 0.0 {
                idx = Some(j);
                break;
            }
        }
        let idx = match idx {
            None => return 0.0,
            Some(idx) => idx,
        };
        a.swap(i, idx);
        if i != idx {
            result = -result;
        }
        for j in i + 1..n {
            let factor = a[j][i] / a[i][i];
            for k in i..n {
                let tmp = factor * a[i][k];
                a[j][k] -= tmp;
            }
        }
        result *= a[i][i];
    }
    result
}
