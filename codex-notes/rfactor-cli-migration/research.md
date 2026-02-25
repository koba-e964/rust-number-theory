# Research: rfactor CLI library migration

## Scope
- Task: migrate `rfactor` binary from current CLI parsing library to a smaller one while preserving behavior.
- Primary target binary: `rfactor` (`src/bin/rfactor.rs`).

## Relevant files and modules
- `src/bin/rfactor.rs`
- `Cargo.toml`
- `README.md`
- `.github/workflows/release-rfactor.yml`
- `target/release/rfactor.map` (local analysis artifact)

## Current execution flow and call graph
1. Process starts in `main()` in `src/bin/rfactor.rs`.
2. CLI is parsed via `clap::Parser` derive on `Cli`:
   - Positional optional `integer: Option<String>`
   - Flag `--verbose` / `-v`
   - Flag `--json`
3. Input selection:
   - If positional integer exists, use it.
   - Else prompt `> ` and read one line from stdin, trim it.
4. Parse selected input as `BigInt` via `BigInt::from_str(...).unwrap()`.
5. Run factorization: `ecm_parallel::factorize_verbose(&value, cli.verbose)`.
6. Format output via `present(cli, result, ecm_stats, elapsed)`:
   - JSON mode:
     - prints `{ entries: [...], stats?: ... }` pretty-printed.
     - `stats` included only if `verbose`.
   - Non-JSON mode:
     - prints factors space-separated with multiplicity expansion.
   - trailing newline is always printed.

## Data structures and invariants
- `Cli` carries user-facing parse state (`integer`, `verbose`, `json`).
- `result: Vec<(BigInt, u64)>` from factorization:
  - each tuple is `(factor, exponent)`.
- JSON output shape:
  - `entries`: array of `{ p: String, e: u64, is_composite?: bool }`.
  - `is_composite` omitted when false (current code always false).
  - optional `stats` object only when verbose.
- Non-JSON output invariant:
  - factors are printed repeated `e` times, separated by single spaces.

## Existing architectural and coding patterns
- Direct `unwrap()`/`panic!` used for parse and I/O errors in this binary.
- CLI parsing currently depends on derive macros (`clap` derive).
- Output formatting is split into a helper (`present`).
- Minimal abstraction around CLI parsing (single `Cli::parse()` call).

## Naming conventions
- CLI struct named `Cli`.
- Booleans are positive flags (`verbose`, `json`).
- Presentation helper named `present`.

## Error handling patterns
- Invalid integer parse currently panics via `unwrap()`.
- stdin read error currently panics.
- No custom exit codes/messages for parse errors besides parser defaults.

## Release and distribution constraints
- `rfactor` is built/released in CI workflow:
  - `cargo build --release --target=... --bin rfactor`
  - runnable check executes `rfactor 12345` on most targets.
- README documents `cargo install ... --bin rfactor` and Homebrew tap install.
- Any migration should keep basic CLI usage stable for existing users.

## Baseline size observations (already measured)
- Binary file size: about `1.2 MiB` (`target/release/rfactor`, arm64 macOS local build).
- `cargo bloat --release --bin rfactor --crates` indicates `.text` dominated by:
  - `std`: 272.6 KiB (38.4% of `.text`)
  - `clap_builder`: 269.2 KiB (37.9% of `.text`)
  - `num_bigint`: 66.4 KiB (9.3% of `.text`)
  - `rust_number_theory`: 62.7 KiB (8.8% of `.text`)
- Top large functions include `clap_builder` parser/build paths and ECM functions.

## Potential migration pitfalls
- Behavioral drift in argument parsing:
  - preserving optional positional + two flags (`--json`, `--verbose`) and `-v` alias.
  - handling unknown flags/options may change user-visible errors/exit codes.
- Maintaining stdin fallback behavior when positional argument absent.
- Avoiding accidental changes to JSON shape or factor print formatting.
- Ensuring parser library replacement does not break MSRV (`rust-version = 1.74`).

## Typing conventions and compatibility constraints
- Project uses Rust 2021 edition and MSRV 1.74.
- Current CLI stores raw string then parses to `BigInt` manually.
- Migration can preserve this typed flow and minimize downstream changes.

## Unknowns
- Whether maintainers want strict compatibility for clap-style help/usage text.
- Whether `-v` short flag compatibility is required (currently provided).
- Whether there are external scripts depending on exact error messages.
