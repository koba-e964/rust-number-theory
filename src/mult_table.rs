use num::BigInt;

/// Multiplication table of a ring of integers (or orders).
struct MultTable {
    table: Vec<Vec<Vec<BigInt>>>,
}

impl MultTable {
    pub fn new(table: Vec<Vec<Vec<BigInt>>>) -> Self {
        MultTable { table }
    }
    pub fn deg(&self) -> usize {
        self.table.len()
    }
    pub fn mul(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        debug_assert_eq!(a.len(), b.len());
        let n = a.len();
        debug_assert_eq!(n, self.deg());
        let mut result = vec![0.into(); n];
        for i in 0..n {
            for j in 0..n {
                let prod = &a[i] * &b[j];
                for k in 0..n {
                    result[k] += &prod * &self.table[i][j][k];
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mult_table_works() {
        // Multiplication table for Z[i]
        let table = vec![
            vec![vec![1.into(), 0.into()], vec![0.into(), 1.into()]],
            vec![vec![0.into(), 1.into()], vec![(-1).into(), 0.into()]],
        ];
        let table = MultTable::new(table);
        let a = vec![2.into(), 3.into()];
        let b = vec![4.into(), 1.into()];
        let prod = table.mul(&a, &b);
        assert_eq!(prod, vec![5.into(), 14.into()]); // (2+3i) * (4+i) = 5 + 14i
    }
}
