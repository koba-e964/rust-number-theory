extern crate num;
extern crate rust_number_theory;

use num::BigInt;

use number_theory_linear::hnf::HNF;

fn main() {
    // 3 1
    // 1 1
    let a: Vec<Vec<BigInt>> = vec![vec![3.into(), 1.into()], vec![1.into(), 1.into()]];
    let hnf = HNF::new(&a);
    let kern = HNF::kernel(&a);
    let kern = HNF::new(&kern);
    eprintln!("hnf =\n{}", hnf);
    eprintln!("kernel =\n{}", kern);

    // 0 0
    // 3 1
    // 1 1
    let a: Vec<Vec<BigInt>> = vec![
        vec![0.into(), 0.into()],
        vec![3.into(), 1.into()],
        vec![1.into(), 1.into()],
    ];
    let hnf = HNF::new(&a);
    let kern = HNF::kernel(&a);
    let kern = HNF::new(&kern);
    eprintln!("hnf =\n{}", hnf);
    eprintln!("kernel =\n{}", kern);

    // 1 0
    // 3 1
    // 1 1
    let a: Vec<Vec<BigInt>> = vec![
        vec![1.into(), 0.into()],
        vec![3.into(), 1.into()],
        vec![1.into(), 1.into()],
    ];
    let hnf = HNF::new(&a);
    let kern = HNF::kernel(&a);
    let kern = HNF::new(&kern);
    eprintln!("hnf =\n{}", hnf);
    eprintln!("kernel =\n{}", kern);

    // 5 0
    // 7 0
    // 2 0
    let a: Vec<Vec<BigInt>> = vec![
        vec![5.into(), 0.into()],
        vec![7.into(), 0.into()],
        vec![2.into(), 0.into()],
    ];
    let hnf = HNF::new(&a);
    let kern = HNF::kernel(&a);
    let kern = HNF::new(&kern);
    eprintln!("hnf =\n{}", hnf);
    eprintln!("kernel =\n{}", kern);

    // 1516
    // -154
    // -336
    // -1423
    let a: Vec<Vec<BigInt>> = vec![
        vec![1516.into()],
        vec![(-154).into()],
        vec![(-336).into()],
        vec![(-1423).into()],
    ];
    let hnf = HNF::new(&a);
    let kern = HNF::kernel(&a);
    let kern = HNF::new(&kern);
    eprintln!("hnf =\n{}", hnf);
    eprintln!("kernel =\n{}", kern);
}
