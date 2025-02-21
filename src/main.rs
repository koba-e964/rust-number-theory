use num::bigint::Sign;
use num::BigInt;
use num::ToPrimitive;
use rust_number_theory::poly_mod::factorize_mod_p;
use rust_number_theory::poly_z;
use rust_number_theory::prime_decomp::decompose;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::str::FromStr;

use rust_number_theory::algebraic::Algebraic;
use rust_number_theory::discriminant;
use rust_number_theory::ecm;
use rust_number_theory::integral_basis;
use rust_number_theory::order;
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
    Integer(BigIntBridge),
    #[serde(rename = "polynomial_and_primes")]
    PolynomialAndPrimes {
        polynomial: Vec<BigIntBridge>,
        primes: Vec<BigIntBridge>,
    },
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

#[derive(Debug, Clone)]
struct Opts {
    config: PathBuf,
}

fn opts() -> Opts {
    use bpaf::*;
    let config = positional("YAML|TOML").help("Configuration file");
    construct!(Opts { config })
        .to_options()
        .fallback_to_usage()
        .run()
}

fn main() {
    let opts = opts();
    let content = fs::read_to_string(&opts.config).unwrap();
    let input_config: InputConfig = if let Ok(conf) = serde_yaml::from_str(&content) {
        conf
    } else {
        // A "nasty hack" suggested in https://github.com/toml-rs/toml-rs/issues/390#issuecomment-1095625417
        let config_json = toml::from_str::<serde_json::Value>(&content).unwrap();
        serde_json::from_value(config_json).unwrap()
    };
    for to_find in input_config.to_find {
        if to_find == "resultant" {
            let polys = match input_config.input {
                Input::Polynomials(ref polys) => polys.clone(),
                _ => {
                    eprintln!("resultant accepts polynomials only");
                    continue;
                }
            };
            let p = polynomial_unbridge(polys[0].clone());
            let q = polynomial_unbridge(polys[1].clone());
            let res = resultant(&p, &q);
            let mut result = HashMap::new();
            result.insert("resultant", res.to_string());
            println!("{}", serde_json::to_string_pretty(&result).unwrap());
            continue;
        }
        if to_find == "discriminant" {
            let polys = match input_config.input {
                Input::Polynomials(ref polys) => polys.clone(),
                _ => {
                    eprintln!("discriminant accepts polynomials only");
                    continue;
                }
            };
            let p = polynomial_unbridge(polys[0].clone());
            let disc = discriminant::discriminant(&p);
            let mut result = HashMap::new();
            result.insert("discriminant", disc.to_string());
            println!("{}", serde_json::to_string_pretty(&result).unwrap());
            continue;
        }
        if to_find == "integral_basis" {
            let polys = match input_config.input {
                Input::Polynomials(ref polys) => polys.clone(),
                _ => {
                    eprintln!("integral_basis accepts polynomials only");
                    continue;
                }
            };
            let p = polynomial_unbridge(polys[0].clone());
            let theta = Algebraic::new(p);
            let o = integral_basis::find_integral_basis(&theta);
            let index = order::index(&o, &order::non_monic_initial_order(&theta));
            eprintln!("Z_K = {}", o);
            let mut result = HashMap::new();
            result.insert("reduced_index", index.to_string());
            result.insert("discriminant", o.discriminant(&theta).to_string());
            println!("{}", serde_json::to_string_pretty(&result).unwrap());
            continue;
        }
        if to_find == "factorization" {
            // TODO: employ faster algorithms
            match input_config.input {
                Input::Integer(ref value) => {
                    let value: BigInt = value.clone().into();
                    let result = ecm::factorize(&value);
                    let mut map = serde_json::Map::new();
                    for (p, e) in result {
                        map.insert(p.to_string(), serde_json::Value::Number(e.into()));
                    }
                    println!(
                        "{}",
                        serde_json::to_string_pretty(&serde_json::Value::Object(map)).unwrap(),
                    );
                }
                Input::Polynomials(ref polys) => {
                    let poly = polynomial_unbridge(polys[0].clone());
                    let (content, result) = poly_z::factorize(&poly);
                    #[derive(Serialize)]
                    struct Factor {
                        factor_vec: Vec<BigIntBridge>,
                        factor_str: String,
                        e: usize,
                    }
                    #[derive(Serialize)]
                    struct FactorizationData {
                        content: BigIntBridge,
                        factors: Vec<Factor>,
                    }
                    let mut factors = vec![];
                    for (f, e) in result {
                        let factor_vec: Vec<BigIntBridge> = polynomial_bridge(f.clone());
                        let factor_str = format!("{:?}", f);
                        factors.push(Factor {
                            factor_vec,
                            factor_str,
                            e,
                        });
                    }
                    println!(
                        "{}",
                        serde_json::to_string_pretty(&FactorizationData {
                            content: content.into(),
                            factors,
                        })
                        .unwrap(),
                    );
                }
                _ => {
                    eprintln!("factorization accepts integer or polynomials only");
                    continue;
                }
            };
            continue;
        }
        if to_find == "factorization-mod-p" {
            let (polynomial, primes) = match input_config.input {
                Input::PolynomialAndPrimes {
                    ref polynomial,
                    ref primes,
                } => (
                    polynomial_unbridge(polynomial.clone()),
                    primes
                        .iter()
                        .map(|p| p.clone().into())
                        .collect::<Vec<BigInt>>(),
                ),
                _ => {
                    eprintln!("factorization-mod-p accepts (polynomial, primes) only");
                    continue;
                }
            };
            #[derive(Serialize)]
            struct Factor {
                factor_vec: Vec<BigIntBridge>,
                factor_str: String,
                e: usize,
            }
            #[derive(Serialize)]
            struct FactorizationData {
                modulus: BigIntBridge,
                factors: Vec<Factor>,
            }
            let mut data = vec![];
            for p in primes {
                let mut factors = vec![];
                fn as_usize(a: &BigInt) -> usize {
                    let (sign, digits) = a.to_u64_digits();
                    match sign {
                        Sign::Plus => {}
                        _ => return 0,
                    }
                    if digits.len() >= 2 {
                        return 0;
                    }
                    digits[0].to_usize().unwrap_or(0)
                }
                let result = factorize_mod_p::<BigInt>(&polynomial, &p, as_usize(&p));
                for (f, e) in result {
                    let factor_vec: Vec<BigIntBridge> = polynomial_bridge(f.clone());
                    let factor_str = format!("{:?}", f);
                    factors.push(Factor {
                        factor_vec,
                        factor_str,
                        e,
                    });
                }
                data.push(FactorizationData {
                    modulus: p.into(),
                    factors,
                });
            }
            println!("{}", serde_json::to_string_pretty(&data).unwrap());
            continue;
        }
        if to_find == "prime-decomposition" {
            let (polynomial, primes) = match input_config.input {
                Input::PolynomialAndPrimes {
                    ref polynomial,
                    ref primes,
                } => (
                    polynomial_unbridge(polynomial.clone()),
                    primes
                        .iter()
                        .map(|p| p.clone().into())
                        .collect::<Vec<BigInt>>(),
                ),
                _ => {
                    eprintln!("factorization-mod-p accepts (polynomial, primes) only");
                    continue;
                }
            };
            #[derive(Serialize)]
            struct Factor {
                norm: BigIntBridge,
                e: usize,
            }
            #[derive(Serialize)]
            struct FactorizationData {
                modulus: BigIntBridge,
                factors: Vec<Factor>,
            }
            let mut data = vec![];
            for p in primes {
                let mut factors = vec![];
                let theta = Algebraic::new(polynomial.clone());
                let int_basis = integral_basis::find_integral_basis(&theta);
                let mult_table = int_basis.get_mult_table(&theta);
                let result = decompose(&theta, &int_basis, &mult_table, &p);
                for (id, e) in result {
                    factors.push(Factor {
                        norm: id.norm().into(),
                        e,
                    });
                }
                data.push(FactorizationData {
                    modulus: p.into(),
                    factors,
                });
            }
            println!("{}", serde_json::to_string_pretty(&data).unwrap());
            continue;
        }
        eprintln!("Unrecognized command: {}", to_find);
    }
}
