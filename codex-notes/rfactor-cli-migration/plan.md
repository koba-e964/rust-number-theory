# Plan: Migrate `rfactor` to a smaller CLI parser

## Overview
Replace `clap` usage in `src/bin/rfactor.rs` with a lightweight parser to reduce binary size while preserving current CLI behavior (`[integer]`, `-v/--verbose`, `--json`, stdin fallback). Validate with `cargo bloat` and binary size checks before/after.

Recommended target parser: `pico-args`.
Reason: very small runtime footprint, straightforward flag/positional parsing, no derive macro expansion overhead.

## Files to change
- `Cargo.toml`
- `src/bin/rfactor.rs`
- `codex-notes/rfactor-cli-migration/plan.md` (checklist progress updates during implementation)

## Detailed implementation steps
1. Update dependencies in `Cargo.toml`.
- Remove `clap` dependency entry.
- Add `pico-args` dependency (minimal feature usage).

2. Replace CLI parsing in `src/bin/rfactor.rs`.
- Remove `use clap::Parser` and derive-based `Cli` annotations.
- Keep `Cli` struct but make it plain Rust data.
- Add a dedicated parse function (e.g., `parse_cli() -> Result<Cli, String>`):
  - Parse `-v` and `--verbose` as equivalent.
  - Parse `--json`.
  - Accept at most one free positional integer argument.
  - Reject unknown arguments with clear error text.
  - Preserve default behavior: when no positional integer, read from stdin.
- In `main`, map parse errors to stderr + exit code 2.

3. Preserve output and factorization behavior.
- Keep `BigInt` parse + factorization flow unchanged.
- Keep `present(...)` output format unchanged for both JSON and non-JSON modes.

4. Validate behavior and size.
- Functional checks:
  - `cargo run --bin rfactor -- 12345`
  - `cargo run --bin rfactor -- --json 12345`
  - `printf '12345\n' | cargo run --bin rfactor`
  - Unknown flag error path check.
- Size checks:
  - `cargo bloat --release --bin rfactor --crates -n 20`
  - `ls -lh target/release/rfactor`
- Compare results against current baseline and summarize deltas.

## Alternatives considered
1. `lexopt`
- Pros: low-level, tiny, expressive parser API.
- Cons: slightly more manual branching than needed for this simple CLI.
- Rejected in favor of `pico-args` for simpler positional/flag extraction.

2. `argh`
- Pros: simpler than clap derive in many cases.
- Cons: derive/proc-macro style may add overhead and reduce size gains.
- Rejected because goal is minimizing binary size impact.

3. Keep `clap` and tune build flags only
- Pros: zero behavior migration risk.
- Cons: leaves the biggest identified contributor in place.
- Rejected because user explicitly requested CLI library migration for size.

## Risks
- CLI compatibility drift (especially short/long flag handling and error text).
- Different exit code semantics from current parser behavior.
- Accidental acceptance/rejection of unexpected positional argument patterns.

## Test strategy
- Run existing build/tests affected by `rfactor` target.
- Execute explicit CLI behavior checks listed above.
- Re-run `cargo bloat` and compare crate-level contribution changes.
- Confirm output format stability by spot-checking text and JSON outputs.

## Assumptions
- Exact `clap` help text compatibility is not required.
- Preserving user-visible invocation semantics is required.
- `rust-version = 1.74` compatibility must be maintained.

## Open questions
- Should we preserve a `--help` output shape similar to current behavior, or is minimal custom help acceptable?
- Should parse failures keep panic-like behavior for invalid integers, or be normalized to clean error messages + exit code?

## Implementation Checklist
- [x] Remove `clap` and add `pico-args` in `Cargo.toml`.
- [x] Replace derive-based CLI parsing in `src/bin/rfactor.rs`.
- [x] Preserve existing stdin fallback and output formatting behavior.
- [x] Run functional CLI checks for positional, flags, stdin, and unknown args.
- [x] Run size checks (`cargo bloat` + binary file size) and record deltas.
- [x] Update checklist with completed items during implementation.
