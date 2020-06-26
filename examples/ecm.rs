extern crate num;
extern crate rust_number_theory;

use num::BigInt;

use rust_number_theory::ecm;

fn main() {
    let large1 = BigInt::from(1_000_000_007);
    let large2 = BigInt::from(1_000_000_009);
    let ans = ecm::ecm(
        &(large1 * large2),
        ecm::ECMConfig {
            b1: 100,
            b2: 10000,
            verbose: true,
        },
    );
    eprintln!("{:?}", ans);
}
