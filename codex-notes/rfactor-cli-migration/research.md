# Research: rfactor CLI parser size migration (bpaf -> pico-args)

## Scope
- Evaluate current `rfactor` CLI parsing on `origin/master` and migrate to a smaller parser while preserving behavior.
- Focus binary: `src/bin/rfactor.rs` (`[[bin]] name = "rfactor"`).

## Relevant files
- `src/bin/rfactor.rs`
- `Cargo.toml`
- `Cargo.lock`
- `.github/workflows/release-rfactor.yml`
- `README.md`

## Current behavior on origin/master
- CLI parser: `bpaf` derive (`#[derive(bpaf::Bpaf)]`, `#[bpaf(options)]`).
- Supported args:
  - Optional positional integer.
  - `-v` / `--verbose`.
  - `--json`.
- If no positional integer is provided, it prompts and reads one line from stdin.
- Factorization path:
  - Parse string to `BigInt`.
  - Call `ecm_parallel::factorize_verbose(&value, cli.verbose)`.
- Output path:
  - JSON mode prints `{ entries, stats? }` (pretty JSON).
  - Text mode prints factors with multiplicity expanded and space separated.

## Constraints and invariants
- Preserve existing invocation semantics for `rfactor` users.
- Preserve stdin fallback and output format.
- Keep repository MSRV (`rust-version = "1.74"`).
- Keep release build compatibility across CI targets.

## Size baseline (origin/master, before migration)
From `cargo bloat --release --bin rfactor --crates -n 20`:
- `.text`: `514.4 KiB`
- file size: `948.6 KiB`
- notable contributors:
  - `std`: `262.6 KiB`
  - `bpaf`: `93.9 KiB`
  - `num_bigint`: `66.4 KiB`
  - `rust_number_theory`: `62.7 KiB`

## Risk areas
- CLI behavior drift for unknown/extra arguments.
- Help output changes from derive parser defaults.
- Accidental output formatting changes.

## Conclusion
- `bpaf` is a measurable contributor in this binary and is a valid migration target for size reduction.
- A lightweight manual parser (`pico-args`) can preserve behavior with lower footprint.
