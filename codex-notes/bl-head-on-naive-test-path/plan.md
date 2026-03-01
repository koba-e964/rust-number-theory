# Plan: Direct Buchmann-Lenstra Implementation with Ideal-Property Validation

## Overview
Replace production `bl::decompose` with a direct Buchmann-Lenstra/Cohen-6.2-style workflow and remove brute-force `O(p^d)` fallback from runtime decomposition.

Validation will use ideal properties directly:
- each returned factor should be prime in tested cases (norm check),
- multiplying returned factors with exponents should reconstruct `(p)`.

## Files To Change
- `src/prime_decomp/bl.rs`
- `src/prime_decomp/mod.rs` (tests + helper functions)
- `src/ideal.rs` (ideal division API/implementation if required by direct BL)
- `codex-notes/bl-head-on-naive-test-path/plan.md` (progress checklist updates)

## Detailed Implementation Steps
1. Refactor `src/prime_decomp/bl.rs` structure.
   - Separate production decomposition logic from naive helper logic.
   - Keep references/comments aligned with Cohen 6.2.

2. Implement direct BL production path in `bl::decompose`.
   - Build `(p)` and candidate ideals from mod-`p` factor data as BL requires.
   - Use ideal/ring operations and linear algebra primitives as needed to derive final prime ideal powers without exhaustive enumeration over all residues.
   - Ensure production path has no `p^d` brute-force enumeration.
   - If BL steps require ideal quotient operations, implement ideal division in `ideal.rs` (for integral ideals), not only inversion.

3. Remove naive fallback logic entirely.
   - Delete production fallback helpers (`is_power_of_p`, decomposition search, splitting helpers, `fallback_decompose`).
   - Keep `bl.rs` focused on direct BL decomposition only.

4. Add/adjust tests.
   - Keep existing regression tests for:
     - index-dividing-prime dispatch case,
     - `x^3 + 9x + 1` at `p=3` expecting `(3)=P1*P2^2` and both norms 3.
   - Add decomposition validation checks:
     - each returned ideal has prime norm (sufficient for tested cases),
     - product of returned ideals with exponents equals `(p)` as an ideal.
   - Add small test helper(s) to make intent explicit:
     - helper to construct `(p)` as an ideal,
     - helper to multiply decomposition factors with exponents.
   - If ideal division is added, include focused unit tests in `ideal.rs` validating quotient identities on known examples.

5. Verification.
   - Run `cargo test prime_decomp -- --nocapture`.
   - Run `cargo clippy --all-targets --all-features` and fix findings in touched files.

## Alternatives Considered
1. Keep current fallback in production with tighter cutoffs.
   - Rejected: still not direct BL and still pathological complexity.

2. Keep a hidden naive implementation under test cfg.
   - Rejected: no longer necessary and adds maintenance burden.

## Risks
- Translating BL steps into existing abstractions may require careful handling of ideal equality/canonicalization.
- Existing tests may rely on representative ordering of factors.
- If direct BL implementation differs in factor order, tests must compare invariant data (norm/exponent multisets).
- Norm-primality checks are sufficient for requested validation but are not a full primality proof in all settings; tests should scope this clearly.
- Implementing ideal division incorrectly would silently invalidate BL steps; quotient tests must cover nontrivial examples.

## Test Strategy
- Primary: `cargo test prime_decomp -- --nocapture`.
- Include targeted assertions on known difficult cases (index-dividing primes, cubic ramification pattern).
- For each checked case:
  - assert prime-ideal norm condition,
  - build `lhs = product(P_i^{e_i})` and `rhs = (p)` and assert `lhs == rhs`.

## Assumptions
- Existing ideal arithmetic and mod-`p` factorization primitives are sufficient for direct BL implementation.
- Ideal multiplication is implemented (`Mul` for `&Ideal`) and is sufficient for reconstruction checks; ideal division is not required for this validation plan.
- If implementation work reveals missing primitive operations for direct BL itself, add them in `ideal.rs` with focused scope.
- Prefer implementing integral ideal division API in `ideal.rs` when quotient operations are needed by BL (`I / J` over integral ideals).

## Open Questions
- None for this plan revision; test helper functions will be added for clarity.

## Implementation Checklist
- [ ] Refactor `src/prime_decomp/bl.rs` to direct BL-only production logic.
- [ ] Implement direct BL production decomposition without brute-force residue enumeration.
- [ ] Implement ideal division in `src/ideal.rs` if required by the direct BL flow.
- [ ] Remove naive fallback logic (`fallback_decompose` and associated helpers).
- [ ] Update/extend tests for required decomposition expectations.
- [ ] Add helper function(s) in prime decomposition tests for `(p)` construction and decomposition-product reconstruction.
- [ ] Add/extend `ideal.rs` tests for ideal division identities if division is introduced.
- [ ] Add assertions that each factor has prime norm in requested test cases.
- [ ] Add assertions that multiplying factors with exponents reconstructs `(p)`.
- [ ] Run `cargo test prime_decomp -- --nocapture`.
- [ ] Run `cargo clippy --all-targets --all-features` and resolve issues in touched code.
