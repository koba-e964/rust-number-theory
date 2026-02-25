# Plan: Replace bpaf with pico-args in rfactor

## Overview
Migrate `rfactor` from `bpaf` derive parsing to `pico-args` manual parsing, preserving behavior and validating size reduction with `cargo bloat`.

## Files to change
- `Cargo.toml`
- `Cargo.lock`
- `src/bin/rfactor.rs`
- `codex-notes/rfactor-cli-migration/research.md`
- `codex-notes/rfactor-cli-migration/plan.md`

## Implementation steps
1. Dependency swap
- Remove `bpaf`.
- Add `pico-args`.

2. Parser migration in `src/bin/rfactor.rs`
- Remove bpaf derive attributes.
- Add `parse_cli() -> Result<Cli, String>` with:
  - `-v`/`--verbose`
  - `--json`
  - optional positional integer
  - `--help`
  - unknown/extra argument errors

3. Preserve behavior
- Keep stdin fallback when integer is omitted.
- Keep JSON/text output formatting unchanged.

4. Validate
- Functional runs for positional, JSON, stdin, help, and unknown flag cases.
- Re-run `cargo bloat --release --bin rfactor --crates` and compare to baseline.

## Risks
- Slightly different help/error phrasing vs bpaf defaults.
- Need careful unknown-argument detection for hyphen-prefixed values.

## Test strategy
- `cargo run --bin rfactor -- 12345`
- `cargo run --bin rfactor -- --json 12345`
- `printf '12345\n' | cargo run --bin rfactor`
- `cargo run --bin rfactor -- --help`
- `cargo run --bin rfactor -- --unknown 12345`
- `cargo bloat --release --bin rfactor --crates -n 20`
- `cargo build --release --bin rfactor`

## Implementation Checklist
- [x] Swap dependencies (`bpaf` -> `pico-args`).
- [x] Replace parser implementation in `src/bin/rfactor.rs`.
- [x] Preserve stdin fallback and output formatting.
- [x] Run functional CLI checks.
- [x] Run release size checks and confirm size improvement.
