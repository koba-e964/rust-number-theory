use rust_number_theory::algebraic::Algebraic;
use rust_number_theory::order::Order;
use rust_number_theory::polynomial::Polynomial;
use rust_number_theory::{hnf::HNF, ideal::Ideal};

fn main() {
    // Z[sqrt(-5)], (2, 1 + sqrt(-5))
    let p = Polynomial::from_raw(vec![5.into(), 0.into(), 1.into()]);
    let theta = Algebraic::new(p);
    let hnf = HNF::hnf(&[
        vec![1.into(), 1.into()],
        vec![5.into(), 1.into()],
        vec![2.into(), 0.into()],
        vec![0.into(), 2.into()],
    ]);
    let o = Order::singly_gen(&theta);
    let mult_table = o.get_mult_table(&theta);
    let x = Ideal::new(hnf, &mult_table);
    assert_eq!(x.norm(), 2.into());
    eprintln!("x = {:?}", x);
    let two = &x * &x; // (2, 1 + sqrt(-5))^2 = (2)
    assert_eq!(two.norm(), 4.into());
}
