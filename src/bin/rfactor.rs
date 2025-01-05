use std::io;
use std::io::Write;
use std::str::FromStr;

use clap::Parser;
use num::BigInt;
use rust_number_theory::ecm;

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
    let result = ecm::factorize(&value);
    present(cli, result);
}

fn present(cli: Cli, result: Vec<(BigInt, u64)>) {
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
        println!("{}", serde_json::to_string_pretty(&entries).unwrap(),);
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
