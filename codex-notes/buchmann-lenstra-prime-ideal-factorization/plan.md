# Plan: Buchmann-Lenstra Prime Ideal Factorization Completion

## Overview
Implement the missing Buchmann-Lenstra decomposition path and integrate routing logic so `prime_decomp::decompose` works when `p | (Z_K : Z[theta])`.

Approach:
1. Implement `bl::decompose` with a mathematically sound decomposition pipeline compatible with existing `Vec<(Ideal, usize)>` output.
2. Route from `prime_decomp::decompose` to `simple` or `bl` based on the index divisibility condition.
3. Add targeted correctness tests, including an index-dividing-prime case that currently panics.
4. Add benchmark(s) to measure decomposition performance on representative inputs.

## Files To Change
- `src/prime_decomp/bl.rs`
- `src/prime_decomp/mod.rs`
- `src/prime_decomp/simple.rs` (tests only, if shared fixture helpers are useful)
- `src/bin/bench-prime-decomp.rs` (new, standard-library timing harness)
- `Cargo.toml` (update `rust-version` if current toolchain target is too old)

## Detailed Implementation Steps
1. Implement `src/prime_decomp/bl.rs::decompose`.
   - Use existing ring/order primitives (`Order`, `MultTable`, `Ideal`) and polynomial factorization helpers.
   - Ensure output format matches `simple::decompose`: list of `(prime_ideal, exponent)`.
   - Preserve compatibility with integral basis coordinate conversion and existing ideal representation.

2. Add algorithm dispatch in `src/prime_decomp/mod.rs`.
   - Compute `index = order::index(int_basis, order::trivial_order_monic(theta))`.
   - If `index % p == 0`, call `bl::decompose`; otherwise call `simple::decompose`.
   - Keep public API unchanged.

3. Add/expand tests for correctness and regression.
   - Non-dividing-index case: decomposition remains consistent with current behavior.
   - Dividing-index case: no panic; decomposition result satisfies structural expectations:
     - ideals are non-zero,
     - norms multiply to `p^n` with exponents,
     - each listed exponent is positive.
   - Add at least one deterministic example where `simple` would panic but BL path succeeds.

4. Add benchmark harness and decomposition benchmarks.
   - Add a small benchmark executable using `std::time::Instant` on stable Rust.
   - Add benchmarks for:
     - baseline simple-path input (`p ∤ index`),
     - BL-trigger input (`p | index`).
   - Benchmark loop measures only decomposition call path (setup done once before timing).
   - Collect per-iteration durations and report:
     - mean runtime,
     - sample standard deviation.

5. Verify with test + bench compilation.
   - `cargo test` (or targeted test modules if full suite is too heavy).
   - `cargo run --release --bin bench-prime-decomp` to collect benchmark output.

6. Align stable toolchain target.
   - Check current compiler baseline in `Cargo.toml` (`rust-version`).
   - If older than about one year, bump to an approximately one-year-old stable release and ensure code/tests/bench still pass on that baseline.

## Alternatives Considered
1. Keep using `simple::decompose` and suppress panic with fallback behavior.
   - Rejected because it is mathematically invalid for `p | (Z_K : Z[theta])` and preserves incorrect decomposition semantics.

2. Implement only dispatch without BL internals (temporary TODO fallback).
   - Rejected because user requested implementation of missing parts and this would still fail for critical cases.

3. Add benchmarks later in a separate PR.
   - Rejected because user explicitly requested benchmarks as part of this work.

## Risks
- Mathematical correctness risk if BL step translation is incomplete.
- Behavior mismatch risk (ordering/content differences) relative to existing callers/tests.
- Benchmark instability due to random or variable setup costs if fixtures are not deterministic.

## Test Strategy
- Unit tests under `src/prime_decomp/` for decomposition path behavior.
- Deterministic cases with explicit polynomial and primes.
- Validate invariants on decomposition output rather than exact ideal basis rows where unnecessary.
- Compile-check benchmark target to avoid CI breakage.

## Benchmark Strategy
- Use deterministic test polynomials and fixed primes.
- Separate setup and measurement where practical:
  - precompute `theta`, `int_basis`, `mult_table` outside per-iteration closure.
- Report decomposition throughput/time for both simple path and BL-trigger path.
- Use repeated runs with fixed iteration counts and print:
  - total elapsed,
  - mean time per iteration,
  - standard deviation per iteration.

## Assumptions
- Existing helper APIs (`ideal`, `order`, `poly_mod`) are sufficient for a first complete BL implementation.
- Project can move from `rust-version = 1.74` to a newer stable baseline (about one year old) if needed for your environment preference.
- Benchmark runtime is local-dev focused; CI may only compile benches.

## Open Questions
- Should benchmark coverage include multiple degree classes (e.g., quadratic + cubic), or is one per path sufficient for this change?

## Implementation Checklist
- [ ] Implement `src/prime_decomp/bl.rs::decompose`.
- [ ] Add BL/simple dispatch in `src/prime_decomp/mod.rs` based on `p | (Z_K : Z[theta])`.
- [ ] Add decomposition regression tests for both path types.
- [ ] Add deterministic test case that previously panicked in simple path.
- [ ] Add `src/bin/bench-prime-decomp.rs` with standard-library timing for simple and BL paths.
- [ ] Compute and print mean and standard deviation for each benchmark scenario.
- [ ] Update `Cargo.toml` `rust-version` to an approximately one-year-old stable release if the current value is too old.
- [ ] Run test verification (`cargo test` or targeted equivalent).
- [ ] Run benchmark executable in release mode and capture timings.
