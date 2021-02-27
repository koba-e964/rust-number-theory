use num::BigInt;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::str::FromStr;

use rust_number_theory::discriminant;
use rust_number_theory::polynomial::Polynomial;
use rust_number_theory::resultant::resultant;

#[derive(Serialize, Deserialize, Debug, Clone)]
struct InputConfig {
    input: Input,
    to_find: Vec<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "lowercase")]
enum Input {
    Polynomials(Vec<Vec<BigIntBridge>>),
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
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

#[allow(unused)]
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
    let args: Vec<String> = std::env::args().collect();
    if args.len() <= 1 {
        eprintln!("This executable needs a configuration file.");
        std::process::exit(1);
    }
    let filename = args[1].clone();
    let file = std::fs::File::open(&filename).unwrap();
    let input_config: InputConfig = serde_yaml::from_reader(&file).unwrap();
    for to_find in input_config.to_find {
        if to_find == "resultant" {
            let polys = match input_config.input {
                Input::Polynomials(ref polys) => polys.clone(),
            };
            let p = polynomial_unbridge(polys[0].clone());
            let q = polynomial_unbridge(polys[1].clone());
            let res = resultant(&p, &q);
            let mut result = HashMap::new();
            result.insert("resultant", res.to_string());
            println!("{}", serde_json::to_string(&result).unwrap());
        }
        if to_find == "discriminant" {
            let polys = match input_config.input {
                Input::Polynomials(ref polys) => polys.clone(),
            };
            let p = polynomial_unbridge(polys[0].clone());
            let disc = discriminant::discriminant(&p);
            let mut result = HashMap::new();
            result.insert("discriminant", disc.to_string());
            println!("{}", serde_json::to_string(&result).unwrap());
        }
    }
}
