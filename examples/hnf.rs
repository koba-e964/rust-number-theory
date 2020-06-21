extern crate num;
extern crate rust_number_theory;

use num::BigInt;

use rust_number_theory::hnf;

fn main() {
    // 3 1
    // 1 1
    let a: Vec<Vec<BigInt>> = vec![vec![3.into(), 1.into()], vec![1.into(), 1.into()]];
    let hnf = hnf::hnf(&a);
    eprintln!("hnf =\n{}", hnf);

    // 0 3 1
    // 0 1 1
    let a: Vec<Vec<BigInt>> = vec![
        vec![0.into(), 3.into(), 1.into()],
        vec![0.into(), 1.into(), 1.into()],
    ];
    let hnf = hnf::hnf(&a);
    eprintln!("hnf =\n{}", hnf);

    // 1 3 1
    // 0 1 1
    let a: Vec<Vec<BigInt>> = vec![
        vec![1.into(), 3.into(), 1.into()],
        vec![0.into(), 1.into(), 1.into()],
    ];
    let hnf = hnf::hnf(&a);
    eprintln!("hnf =\n{}", hnf);
}
