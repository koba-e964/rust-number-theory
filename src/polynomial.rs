extern crate num;
use num::BigInt;

// The leading coefficient (a[a.len() - 1]) must not be 0.
#[derive(Clone)]
pub struct Polynomial {
    pub dat: Vec<BigInt>
}

impl Polynomial {
    pub fn is_zero(&self) -> bool {
        self.dat.len() == 0
    }
    // If 0, returns usize::max_value().
    pub fn deg(&self) -> usize {
        if self.dat.len() == 0 {
            usize::max_value()
        } else {
            self.dat.len() - 1
        }
    }
    pub fn from_raw(mut raw: Vec<BigInt>) -> Self {
        let mut ma = 0;
        for i in (0 .. raw.len()).rev() {
            if raw[i] != 0.into() {
                ma = i + 1;
                break;
            }
        }
        raw.drain(ma..);
        Polynomial { dat: raw }
    }
    pub fn from_int(v: impl Into<BigInt>) -> Self {
        Polynomial::from_raw(vec![v.into()])
    }
}

impl std::fmt::Debug for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        (self as &std::fmt::Display).fmt(f)
    }
}
impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_zero() {
            return 0.fmt(f);
        }
        let d = self.deg();
        let mut term_appear = false;
        for i in (0 .. d + 1).rev() {
            if self.dat[i] == 0.into() { continue; }
            if term_appear {
                write!(f, " + ")?;
            }
            write!(f, "({})", self.dat[i])?;
            if i >= 2 {
                write!(f, "X^{}", i)?;
            } else if i == 1 {
                write!(f, "X")?;
            }
            term_appear = true;
        }
        Ok(())
    }
}


// b should be monic, or it may panic.
pub fn div_rem(a: &Polynomial, b: &Polynomial) -> (Polynomial, Polynomial) {
    let a_deg = a.deg();
    let b_deg = b.deg();
    if b_deg == usize::max_value() || a_deg == usize::max_value()
        || a_deg < b_deg {
            return (Polynomial::from_int(0), a.clone());
        }
    assert_eq!(b.dat[b_deg], 1.into());
    let mut tmp = a.dat.clone();
    let mut quo = vec![0.into(); a_deg - b_deg + 1];
    // Naive division
    for i in (0 .. a_deg - b_deg + 1).rev() {
        let mut coef = tmp[i + b_deg].clone();
        for j in 0 .. b_deg + 1 {
            tmp[i + j] -= &coef * &b.dat[j];
        }
        quo[i] = coef;
    }
    (Polynomial::from_raw(quo), Polynomial::from_raw(tmp))
}
