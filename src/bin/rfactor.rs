use std::io::Write;
use std::str::FromStr;
use std::{io, time::Instant};

use clap::Parser;
use num::BigInt;
use rust_number_theory::ecm::EcmStats;
use rust_number_theory::ecm_parallel;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Optional integer argument to factorize
    integer: Option<String>,
    #[arg(short, long)]
    verbose: bool,
    #[arg(long)]
    json: bool,
}

fn main() {
    let cli = Cli::parse();

    // You can check the value provided by positional arguments, or option arguments
    let value = if let Some(integer) = cli.integer.as_deref() {
        integer.to_string()
    } else {
        // Reads from stdin
        print!("> ");
        io::stdout().flush().ok().unwrap();
        let mut s = "".to_string();
        match io::stdin().read_line(&mut s) {
            Ok(_) => {}
            Err(err) => {
                panic!("{err}");
            }
        }
        s = s.trim().to_string();
        s
    };

    let value = BigInt::from_str(&value).unwrap();
    let start = Instant::now();
    let (result, ecm_stats) = ecm_parallel::factorize_verbose(&value, cli.verbose);
    let elapsed = start.elapsed();
    present(cli, result, ecm_stats, elapsed);
}

fn present(
    cli: Cli,
    result: Vec<(BigInt, u64)>,
    ecm_stats: EcmStats,
    elapsed: std::time::Duration,
) {
    if cli.json {
        #[derive(serde::Serialize)]
        struct Entry {
            p: String,
            e: u64,
            #[serde(skip_serializing_if = "std::ops::Not::not")]
            is_composite: bool,
        }
        let mut entries = Vec::new();
        for (p, e) in result {
            entries.push(Entry {
                p: p.to_string(),
                e,
                is_composite: false,
            });
        }
        let mut object = serde_json::Map::new();
        object.insert(
            "entries".to_string(),
            serde_json::Value::Array(
                entries
                    .into_iter()
                    .map(|x| serde_json::to_value(x).unwrap())
                    .collect(),
            ),
        );
        if cli.verbose {
            #[derive(serde::Serialize)]
            struct Data {
                curve_count: u64,
                time: f64,
                average: f64,
            }
            let data = Data {
                curve_count: ecm_stats.curve_count,
                time: elapsed.as_secs_f64(),
                average: elapsed.as_secs_f64() / ecm_stats.curve_count as f64,
            };
            object.insert("stats".to_string(), serde_json::to_value(data).unwrap());
        }
        println!("{}", serde_json::to_string_pretty(&object).unwrap(),);
    } else {
        let mut first = true;
        for (p, e) in result {
            if !first {
                print!(" ");
            }
            first = false;
            for i in 0..e {
                print!("{p}");
                if i + 1 < e {
                    print!(" ");
                }
            }
        }
    }
    println!();
}
