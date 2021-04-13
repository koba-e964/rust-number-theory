use crate::order::Order;
use num::{BigInt, Complex, ToPrimitive};

/// Embeddings into R or C.
#[derive(Clone, Debug)]
pub struct CEmbeddings {
    // #real embeddings
    r: usize,
    // #complex embeddings
    s: usize,
    basis: Vec<Vec<Complex<f64>>>,
}

impl CEmbeddings {
    #[allow(clippy::needless_range_loop)]
    pub fn new(roots_re: &[f64], roots_im: &[Complex<f64>], o: &Order) -> Self {
        let r = roots_re.len();
        let s = roots_im.len();
        let deg = r + 2 * s;
        let mut basis = vec![vec![Complex::new(0.0, 0.0); deg]; r + s];
        for i in 0..r + s {
            let root = if i < r {
                roots_re[i].into()
            } else {
                roots_im[i - r]
            };
            for j in 0..deg {
                let mut current = Complex::new(1.0, 0.0);
                for k in 0..deg {
                    basis[i][j] += o.basis[j][k].to_f64().unwrap() * current;
                    current *= root;
                }
            }
        }
        Self { r, s, basis }
    }
    pub fn deg(&self) -> usize {
        self.basis[0].len()
    }
    /// Computes the idx-th embedding of num[0]*w[0] + ... num[n-1]*w[n-1].
    #[allow(clippy::needless_range_loop)]
    pub fn compute(&self, idx: usize, num: &[BigInt]) -> Complex<f64> {
        let mut val = Complex::new(0.0, 0.0);
        let deg = self.deg();
        for k in 0..deg {
            val += self.basis[idx][k] * num[k].to_f64().unwrap();
        }
        val
    }
    /// Returns the idx-th embedding of w[k].
    pub fn get(&self, idx: usize, k: usize) -> Complex<f64> {
        self.basis[idx][k]
    }
    /// Returns the number of real embeddings.
    pub fn real(&self) -> usize {
        self.r
    }
    /// Returns the number of complex embeddings / 2.
    pub fn complex(&self) -> usize {
        self.s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebraic::Algebraic;
    use crate::integral_basis::find_integral_basis;
    use crate::numerical_roots::find_roots_reim;
    use crate::polynomial::Polynomial;
    use num::{BigInt, ToPrimitive};
    use num::{BigRational, One, Zero};

    #[test]
    fn cembeddings_works_1() {
        // 1 \pm 6i
        let roots_re = vec![];
        let roots_im = vec![Complex::new(1.0, 6.0)];
        let order = Order {
            basis: vec![
                vec![BigRational::one(), BigRational::zero()],
                vec![BigRational::zero(), BigRational::one()],
            ],
        };
        let embeddings = CEmbeddings::new(&roots_re, &roots_im, &order);
        assert_eq!(
            embeddings.compute(0, &[BigInt::from(0), BigInt::from(1)]),
            Complex::new(1.0, 6.0)
        );
        // 2 + 3 * (1 + 6i) = 5 + 18i
        assert_eq!(
            embeddings.compute(0, &[BigInt::from(2), BigInt::from(3)]),
            Complex::new(5.0, 18.0)
        );
    }
    #[test]
    fn cembeddings_works_2() {
        // One real root, four imaginary roots.
        let poly_vec: Vec<BigInt> = vec![
            31.into(),
            36.into(),
            27.into(),
            (-4).into(),
            9.into(),
            1.into(),
        ];
        let poly = Polynomial::from_raw(poly_vec.to_vec());
        let poly_complex =
            Polynomial::from_raw(poly_vec.iter().map(|b| b.to_f64().unwrap()).collect());
        let theta = Algebraic::new(poly);
        let o = find_integral_basis(&theta);

        // Find embeddings and roots of unity
        let (roots_re, roots_im) = find_roots_reim(poly_complex);
        let basis = CEmbeddings::new(&roots_re, &roots_im, &o);
        assert_eq!(basis.real(), 1);
        assert_eq!(basis.complex(), 2);
    }
}
