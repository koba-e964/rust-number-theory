use std::io::Write;
use std::str::FromStr;
use std::{io, time::Instant};

use num::BigInt;
use rust_number_theory::ecm::EcmStats;
use rust_number_theory::ecm_parallel;

struct Cli {
    verbose: bool,
    json: bool,
    integer: Option<String>,
}

fn main() {
    let cli = match parse_cli() {
        Ok(cli) => cli,
        Err(err) => {
            eprintln!("{err}");
            eprintln!();
            eprintln!("Use --help for usage.");
            std::process::exit(2);
        }
    };

    let value = if let Some(integer) = cli.integer.as_deref() {
        integer.to_string()
    } else {
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

fn parse_cli() -> Result<Cli, String> {
    let mut args = pico_args::Arguments::from_env();
    if args.contains(["-h", "--help"]) {
        print_help();
        std::process::exit(0);
    }
    let verbose = args.contains(["-v", "--verbose"]);
    let json = args.contains("--json");
    let mut free = Vec::new();
    for arg in args.finish() {
        let arg = arg
            .into_string()
            .map_err(|arg| format!("non-UTF-8 argument: {:?}", arg))?;
        free.push(arg);
    }
    if let Some(arg) = free.iter().find(|arg| looks_like_option(arg)) {
        return Err(format!("unknown argument: {arg}"));
    }
    if free.len() > 1 {
        return Err(format!(
            "unexpected extra positional argument: {}",
            free[1]
        ));
    }
    let integer = free.pop();
    Ok(Cli {
        integer,
        verbose,
        json,
    })
}

fn looks_like_option(arg: &str) -> bool {
    if !arg.starts_with('-') || arg == "-" {
        return false;
    }
    let bytes = arg.as_bytes();
    if bytes.len() >= 2 && bytes[1].is_ascii_digit() {
        return false;
    }
    true
}

fn print_help() {
    println!("rfactor {}", env!("CARGO_PKG_VERSION"));
    println!("Factorize an integer.");
    println!();
    println!("Usage:");
    println!("  rfactor [OPTIONS] [integer]");
    println!();
    println!("Options:");
    println!("  -v, --verbose   Show verbose ECM statistics");
    println!("      --json      Output in JSON format");
    println!("  -h, --help      Print help");
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
