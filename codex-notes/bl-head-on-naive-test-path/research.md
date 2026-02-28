# Research: Replace Brute-Force Fallback with Direct Buchmann-Lenstra, Keep Naive Path for Testing

## Task Scope
User request: implement Buchmann-Lenstra decomposition "head-on" (no production fallback enumeration), while preserving the current naive fallback as a testing aid.

## Relevant Files and Modules
- `src/prime_decomp/mod.rs`
  - Public dispatch entrypoint.
  - Routes to `simple::decompose` when `p ∤ (Z_K : Z[theta])`, else to `bl::decompose`.
- `src/prime_decomp/simple.rs`
  - Cohen 4.8.2-style decomposition for the easy index-coprime case.
- `src/prime_decomp/bl.rs`
  - Current BL-labeled path.
  - Contains direct `(p, f_i(theta))` construction plus `fallback_decompose` brute-force enumeration over `O/pO` residues.
- `src/ideal.rs`
  - Ideal representation and arithmetic (HNF-based lattice model).
- `number-theory-linear/src/subspace.rs`
  - `image_mod_p` and other linear algebra primitives over finite fields/rationals.
  - Historically relevant for Cohen 6.2.5-style multiplication in `O/pO`.

## Current Execution Flow / Call Graph
1. CLI / library users call `prime_decomp::decompose(theta, int_basis, mult_table, p)`.
2. `prime_decomp::mod::decompose` computes index condition:
   - if `index % p != 0` -> `simple::decompose`.
   - else -> `bl::decompose`.
3. In current `bl::decompose`:
   - Build `pz = (p)`.
   - Factor minimal polynomial mod `p`.
   - Build ideals `(p, f_i(theta))` from factors.
   - If norm product and trivial checks pass, return this result.
   - Else run `fallback_decompose`:
     - enumerate all nonzero vectors in `(Z/pZ)^d` (bounded by hard cutoff),
     - build candidates `(p, a)`,
     - solve combinatorially for product equal to `(p)`,
     - refine some `p^f` norm factors by splitting against norm-`p` candidates.

## Data Structures and Invariants
- Ideals are represented by `Ideal { hnf, mult_table }`; equality is structural via HNF canonicalization.
- `Ideal::norm()` is determinant of HNF.
- Decomposition output expects vector of `(ideal, exponent)` with positive exponents.
- For correct prime decomposition, multiplicative norm relation should satisfy:
  - `prod_i N(P_i)^{e_i} = p^d`, where `d = [K:Q]`.
- In production BL path, algorithmic correctness should not depend on exhaustive enumeration of all residue vectors.

## Observed Problems in Current BL Path
- `fallback_decompose` has worst-case complexity proportional to `p^d` (enumeration of residue vectors), then additional combinatorial search over candidates.
- For practical bounds such as `p` up to `10^4`, this is infeasible except tiny `d`.
- Current comments already acknowledge it is not a literal implementation of BL/Cohen 6.2 substeps.

## Architectural Patterns
- Module-level split by algorithm (`simple` vs `bl`) under a stable public API.
- Tests are colocated in `src/prime_decomp/mod.rs` and `src/prime_decomp/simple.rs`.
- Mathematical precondition failures are often handled by panic in old code paths, but current API returns vectors directly.

## Error Handling / Typing Conventions
- `BigInt` / `BigRational` across algebraic number theory logic.
- Minimal use of custom `Result` in prime decomposition entrypoints (mostly panic or deterministic return).

## Constraints
- Preserve public signature:
  - `pub fn decompose(...) -> Vec<(Ideal<'mul>, usize)>`.
- Keep existing simple-path behavior intact.
- Remove brute-force fallback from production path.
- Preserve naive method for test-only usage (or explicitly non-production usage), as requested.
- Keep PR scope to single logical change: BL correctness + test-only naive preservation.

## Potential Pitfalls
- Replacing fallback with direct BL may change output ordering and ideal representatives while remaining mathematically equivalent.
- Existing tests with exact shapes may need adaptation if canonicalization differs.
- Care is needed to avoid introducing heavy intermediate conversions that regress performance.

## Unknowns / Gaps to Resolve in Planning
- Exact in-repo representation strategy for `O/pO` ideals (currently no dedicated type).
- Whether to keep naive code under `#[cfg(test)]` in `bl.rs` or extract into test helper module.
- Which BL sub-steps map cleanly onto existing primitives vs requiring new helper functions.

## Summary of Current Reality
- Current BL path is partially direct and partially brute-force.
- User requirement is explicit: production path must be direct Buchmann-Lenstra without fallback enumeration.
- Naive path should remain available only for tests/validation.
