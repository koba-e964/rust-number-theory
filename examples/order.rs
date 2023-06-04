extern crate num;
extern crate rust_number_theory;

use num::{BigInt, BigRational, One, Zero};

use rust_number_theory::algebraic::Algebraic;
use rust_number_theory::order::{
    index, non_monic_initial_order, trivial_order_monic, union, Order,
};
use rust_number_theory::polynomial::Polynomial;

fn main() {
    // O1 = (1, (1 + x) / 2), O2 = (1, (2 + x) / 3)
    // (O1 + O2 : Z[x]) = 6
    let inv2 = BigRational::new(1.into(), 2.into());
    let inv3 = BigRational::new(1.into(), 3.into());
    // (x-1)^2 + 36, the minimum polynomial of 1 + 6i
    let theta = Algebraic::new(Polynomial::from_raw(vec![37.into(), (-2).into(), 1.into()]));
    let o = trivial_order_monic(&theta);
    let o1 = Order::from_basis(&[
        vec![BigRational::one(), BigRational::zero()],
        vec![inv2.clone(), inv2],
    ]);
    let o2 = Order::from_basis(&[
        vec![BigRational::one(), BigRational::zero()],
        vec![inv3.clone() * BigInt::from(2), inv3],
    ]);
    eprintln!("(O1: Z[x]) = {}", index(&o1, &o));
    eprintln!("(O2: Z[x]) = {}", index(&o2, &o));
    let union = union(&o1, &o2);
    eprintln!("(O1 + O2: Z[x]) = {}", index(&union, &o));

    let p = Polynomial::from_raw(vec![
        5.into(),
        6.into(),
        (-7).into(),
        6.into(),
        (-7).into(),
        6.into(),
    ]);
    let theta = Algebraic::new(p.clone());
    let o = non_monic_initial_order(&theta);
    let o1_old: Vec<_> = vec![
        vec![1, 0, 0, 0, 0],
        vec![0, 6, 0, 0, 0],
        vec![0, 5, 6, 0, 0],
    ]
    .into_iter()
    .map(|row| {
        row.into_iter()
            .map(|x| BigRational::new(x.into(), 1.into()))
            .collect()
    })
    .collect();
    let o1_new: Vec<_> = vec![vec![1, 9, 11, 6, 0], vec![1, 2, 10, 5, 6]]
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|x| BigRational::new(x.into(), 2.into()))
                .collect()
        })
        .collect();
    let o1 = Order::from_basis(&{
        let mut v = o1_old;
        let mut o1_new = o1_new;
        v.append(&mut o1_new);
        v
    });
    eprintln!("(O1 : O) = {}", index(&o1, &o));

    let disc = o.discriminant(&theta);
    eprintln!("D(O) = {}", disc);

    let obig = Order::singly_gen(&(&theta * &Algebraic::from_int(p, 6)));
    eprintln!("D(Obig) = {}", obig.discriminant(&theta));
    eprintln!("(O : Obig) = {}", index(&o, &obig));
}
