use num::BigInt;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

use rust_number_theory::discriminant;
use rust_number_theory::polynomial::Polynomial;
use rust_number_theory::resultant::resultant;

#[derive(Serialize, Deserialize, Debug)]
struct InputConfig {
    input: Input,
    to_find: Vec<String>,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
enum Input {
    Polynomials(Vec<Vec<BigIntBridge>>),
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq)]
#[serde(transparent)]
struct BigIntBridge(String);

impl From<BigIntBridge> for BigInt {
    fn from(BigIntBridge(x): BigIntBridge) -> Self {
        BigInt::from_str(&x).expect("invalid number")
    }
}

impl From<BigInt> for BigIntBridge {
    fn from(x: BigInt) -> Self {
        BigIntBridge(x.to_string())
    }
}

fn polynomial_bridge(x: Polynomial<BigInt>) -> Vec<BigIntBridge> {
    let mut dat = vec![];
    for elem in x.dat {
        dat.push(elem.into());
    }
    dat
}

fn polynomial_unbridge(x: Vec<BigIntBridge>) -> Polynomial<BigInt> {
    let mut dat = vec![];
    for elem in x {
        dat.push(elem.into());
    }
    Polynomial { dat }
}

fn main() {
    // 9x^5 + 6x^4 + 2x^2 + 5
    let p: Polynomial<BigInt> = Polynomial::from_raw(vec![
        5.into(),
        0.into(),
        2.into(),
        0.into(),
        6.into(),
        9.into(),
    ]);
    println!(
        "{}",
        serde_json::to_string(&polynomial_bridge(p.clone())).unwrap()
    );
    // 7x^4 + x^3 + 6x^2 + 6x + 6
    let q: Polynomial<BigInt> =
        Polynomial::from_raw(vec![6.into(), 6.into(), 6.into(), 1.into(), 7.into()]);
    let res = resultant(&p, &q);
    eprintln!("resultant = {}", res);
    assert_eq!(res, 335159672.into());
    // 2x^3 + x^2 - 2x + 3
    let p: Polynomial<BigInt> =
        Polynomial::from_raw(vec![3.into(), (-2).into(), 1.into(), 2.into()]);
    eprintln!("det = {}", discriminant::discriminant(&p));

    let filename = "data/input-resultant.yml";
    let file = std::fs::File::open(&filename).unwrap();
    let input_config: InputConfig = serde_yaml::from_reader(&file).unwrap();
    eprintln!("input: {:?}", input_config);
}
