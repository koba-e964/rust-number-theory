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
            let (left, right) = a.split_at_mut(j);
            let row_i = &left[i];
            let row_j = &mut right[0];
            for (a_jk, a_ik) in row_j.iter_mut().skip(i).zip(row_i.iter().skip(i)) {
                let tmp = factor * a_ik;
                *a_jk -= tmp;
            }
        }
        result *= a[i][i];
    }
    result
}
