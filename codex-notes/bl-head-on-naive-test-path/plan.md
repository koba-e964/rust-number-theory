# Plan: Direct Buchmann-Lenstra Implementation with Test-Only Naive Preservation

## Overview
Replace production `bl::decompose` with a direct Buchmann-Lenstra/Cohen-6.2-style workflow and remove brute-force `O(p^d)` fallback from runtime decomposition.

Preserve the current naive decomposition search logic for verification/testing only (behind test cfg), so we can cross-check direct BL outputs on small cases without exposing that path to production calls.

## Files To Change
- `src/prime_decomp/bl.rs`
- `src/prime_decomp/mod.rs` (tests if needed)
- `codex-notes/bl-head-on-naive-test-path/plan.md` (progress checklist updates)

## Detailed Implementation Steps
1. Refactor `src/prime_decomp/bl.rs` structure.
   - Separate production decomposition logic from naive helper logic.
   - Keep references/comments aligned with Cohen 6.2.

2. Implement direct BL production path in `bl::decompose`.
   - Build `(p)` and candidate ideals from mod-`p` factor data as BL requires.
   - Use ideal/ring operations and linear algebra primitives as needed to derive final prime ideal powers without exhaustive enumeration over all residues.
   - Ensure production path has no `p^d` brute-force enumeration.

3. Preserve naive algorithm under test-only scope.
   - Move existing fallback helpers (`is_power_of_p`, decomposition search, splitting helpers) into `#[cfg(test)]` helper section/module.
   - Expose only test helper entrypoints needed to compare direct BL output vs naive output on bounded examples.

4. Add/adjust tests.
   - Keep existing regression tests for:
     - index-dividing-prime dispatch case,
     - `x^3 + 9x + 1` at `p=3` expecting `(3)=P1*P2^2` and both norms 3.
   - Add cross-check test(s) on small inputs where naive test helper is feasible:
     - compare multiset of `(norm, exponent)` between direct BL and naive helper.

5. Verification.
   - Run `cargo test prime_decomp -- --nocapture`.
   - Run `cargo clippy --all-targets --all-features` and fix findings in touched files.

## Alternatives Considered
1. Keep current fallback in production with tighter cutoffs.
   - Rejected: still not direct BL and still pathological complexity.

2. Delete naive implementation entirely.
   - Rejected: user explicitly wants to keep naive logic for testing utility.

3. Hide naive logic in a separate non-test runtime flag.
   - Rejected for now: unnecessary production surface area and risk of accidental use.

## Risks
- Translating BL steps into existing abstractions may require careful handling of ideal equality/canonicalization.
- Existing tests may rely on representative ordering of factors.
- If direct BL implementation differs in factor order, tests must compare invariant data (norm/exponent multisets).

## Test Strategy
- Primary: `cargo test prime_decomp -- --nocapture`.
- Include targeted assertions on known difficult cases (index-dividing primes, cubic ramification pattern).
- Add direct-BL vs naive-helper comparison on bounded input(s) to validate correctness while avoiding production fallback use.

## Assumptions
- Existing ideal arithmetic and mod-`p` factorization primitives are sufficient for direct BL implementation.
- Naive helper can remain test-only without affecting public API.

## Open Questions
- Should we expose a dedicated internal comparison utility for future BL regression tests, or keep helper private to `bl.rs` tests only?

## Implementation Checklist
- [ ] Refactor `src/prime_decomp/bl.rs` to isolate production vs test-only helper logic.
- [ ] Implement direct BL production decomposition without brute-force residue enumeration.
- [ ] Move naive fallback logic behind `#[cfg(test)]` and keep it callable from tests.
- [ ] Update/extend tests for required decomposition expectations.
- [ ] Add at least one BL-vs-naive cross-check test on a small feasible case.
- [ ] Run `cargo test prime_decomp -- --nocapture`.
- [ ] Run `cargo clippy --all-targets --all-features` and resolve issues in touched code.
