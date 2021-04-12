#[derive(Debug)]
pub struct Cholesky {
    q: Vec<Vec<f64>>,
}

impl Cholesky {
    pub fn find(matrix: &[Vec<f64>]) -> Self {
        let n = matrix.len();
        let mut q = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i..n {
                q[i][j] = matrix[i][j];
            }
        }
        for i in 0..n {
            for j in i + 1..n {
                q[j][i] = q[i][j];
                q[i][j] /= q[i][i];
            }
            for k in i + 1..n {
                for l in k..n {
                    q[k][l] -= q[k][i] * q[i][l];
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                q[i][j] = 0.0;
            }
        }
        Self { q }
    }

    pub fn find_value(&self, x: &[f64]) -> f64 {
        let n = self.q.len();
        assert_eq!(x.len(), n);
        let mut val = 0.0;
        for i in 0..n {
            let mut tmp = x[i];
            for j in i + 1..n {
                tmp += self.q[i][j] * x[j];
            }
            val += self.q[i][i] * tmp * tmp;
        }
        val
    }

    /// Algorithm 2.7.5 in [Cohen]
    /// Returns all non-zero vectors x that satisfy (x, Qx) <= c
    /// Only one of the two vectors x, -x is returned.
    /// [Cohen]: Cohen, Henri. A course in computational algebraic number theory. Vol. 138. Springer Science & Business Media, 2013.
    pub fn find_short_vectors(&self, c: f64) -> Vec<(f64, Vec<i64>)> {
        let mut result = vec![];
        Self::dfs(
            &self.q,
            c,
            c,
            self.q.len(),
            &mut vec![0; self.q.len()],
            &mut result,
        )
        .unwrap_err();
        result
    }
    fn dfs(
        q: &[Vec<f64>],
        c: f64,
        rem: f64,
        idx: usize,
        x: &mut [i64],
        result: &mut Vec<(f64, Vec<i64>)>,
    ) -> Result<(), ()> {
        let n = q.len();
        if rem < 0.0 {
            return Ok(());
        }
        if idx == 0 {
            if x.iter().all(|&x| x == 0) {
                return Err(());
            }
            result.push((c - rem, x.to_vec()));
            return Ok(());
        }
        let idx = idx - 1;
        let mut u = 0.0;
        for i in idx + 1..n {
            u += q[idx][i] * x[i] as f64;
        }
        let z = (rem / q[idx][idx]).sqrt();
        let up = (z - u).floor() as i64;
        let down = (-z - u).ceil() as i64;
        for xidx in down..=up {
            let new_rem = rem - q[idx][idx] * (xidx as f64 + u) * (xidx as f64 + u);
            x[idx] = xidx;
            Self::dfs(q, c, new_rem, idx, x, result)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_works() {
        let q = vec![vec![2.0, 1.0], vec![1.0, 1.0]];
        let cho = Cholesky::find(&q);
        for x in -10..10 {
            for y in -10..10 {
                let x = [x as f64, y as f64];
                let val_cho = cho.find_value(&x);
                let mut val = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        val += q[i][j] * x[i] * x[j];
                    }
                }
                assert!((val_cho - val).abs() <= 1e-9);
            }
        }
    }

    #[test]
    fn find_short_vectors_works() {
        let q = vec![vec![2.0, 1.0], vec![1.0, 1.0]];
        let cho = Cholesky::find(&q);
        let c = 3.0;
        let vecs = cho.find_short_vectors(c);
        assert_eq!(vecs.len(), 4);
        for (ret_val, x) in vecs {
            let mut val = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    val += q[i][j] * (x[i] * x[j]) as f64;
                }
            }
            assert!((ret_val - val).abs() <= 1.0e-9);
            assert!(ret_val <= c);
        }
    }
}
