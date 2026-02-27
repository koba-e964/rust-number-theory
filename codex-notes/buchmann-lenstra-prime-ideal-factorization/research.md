# Research: Buchmann-Lenstra Prime Ideal Factorization Gaps

## Task Scope
Implement missing parts of prime ideal decomposition using Buchmann-Lenstra path in this repository.

## Relevant Files and Modules
- `src/prime_decomp/mod.rs`
  - Public entry point `prime_decomp::decompose`.
  - Currently always calls `simple::decompose`.
  - Contains TODO: support the case `p | (Z_K : Z[theta])`.
- `src/prime_decomp/simple.rs`
  - Working decomposition algorithm (Cohen 4.8.2) for the easy case `p ∤ (Z_K : Z[theta])`.
  - Explicitly panics when `p | index`.
- `src/prime_decomp/bl.rs`
  - Intended Buchmann-Lenstra implementation (comment references Cohen 6.2.2).
  - `decompose` is currently `panic!()`.
  - Contains helper `multiply` (Cohen 6.2.5) that multiplies ideals modulo `p` via image computation over `F_p`.
- `src/order.rs`
  - Provides `order::trivial_order_monic(theta)` for `Z[theta]`.
  - Provides `order::index(a, b)` used to compute `(a:b)`.
  - Provides conversion from algebraic expression to integral basis coordinates (`to_z_basis_int`).
- `src/ideal.rs`
  - Defines ideal representation with HNF (`Ideal`).
  - Supports creation by principal generator and ideal sum/product.
  - Norm is HNF determinant.
- `src/poly_mod/factorize_mod_p.rs`
  - Factors polynomial over `F_p`; returns `Vec<(factor, multiplicity)>`.
- `number-theory-linear/src/subspace.rs`
  - `image_mod_p` computes row-image basis mod `p` (used by `bl::multiply`).

## Current Execution Flow and Call Graph
1. CLI path (`src/main.rs`, command `prime-decomposition`) builds:
   - `theta = Algebraic::new(polynomial)`
   - `int_basis = integral_basis::find_integral_basis(&theta)`
   - `mult_table = int_basis.get_mult_table(&theta)`
2. Calls `prime_decomp::decompose(&theta, &int_basis, &mult_table, &p)`.
3. `prime_decomp::decompose` currently routes to `simple::decompose` unconditionally.
4. `simple::decompose` checks `index = (Z_K : Z[theta])`; if `p | index`, panics.

Implication: Prime decomposition fails for valid inputs when the prime divides the index.

## Data Structures and Invariants
- `Order` basis vectors are reduced to HNF-like normal form over rationals.
- `MultTable` assumes basis element `w_0 = 1` in several methods (explicit note in `inv`).
- `Ideal` stores full-rank lattice over the integral basis via integer HNF.
- `Ideal::new(HNF, mult_table)` does not validate ideal closure; callers must preserve invariants.
- `simple::decompose` constructs each prime ideal as `(p, f(theta))` mapped to integral basis coordinates.
- `poly_mod::factorize_mod_p` normalizes factors mod `p` and reports multiplicities.

## Existing Architectural and Coding Patterns
- Prime decomposition code is split by algorithm module (`simple`, `bl`) and routed from `mod.rs`.
- Algebraic elements are represented as polynomials in `theta` and converted to `Z_K` basis via `Order` methods.
- Matrix/lattice canonicalization is done via `HNF::new` and linear-algebra helper crate.
- Panics are used in precondition violations (not Result-based APIs).

## Error Handling and Typing Conventions
- The crate uses panic on invalid preconditions for many math routines.
- No dedicated error enum for prime decomposition currently.
- Integer arithmetic primarily uses `num::BigInt`; rational basis uses `BigRational`.

## Potential Pitfalls
- Switching dispatch from `simple` to `bl` for some inputs can change behavior and output ordering.
- `bl::multiply` operates in `O/pO` semantics but returns `Ideal`; invariants for this representation are implicit.
- Mapping polynomial factors to ideals must respect basis conversion and degree constraints.
- If BL algorithm path depends on assumptions not encoded in current helpers, decomposition may compile but be mathematically incorrect.

## Constraints
- Need to preserve public API: `pub fn decompose(...) -> Vec<(Ideal, usize)>`.
- Must integrate with existing `main.rs` prime-decomposition workflow.
- Keep compatibility with current BigInt/HNF abstractions.

## Unknowns
- No existing tests cover `prime_decomp::bl::decompose` behavior.
- No in-repo reference implementation for Cohen 6.2.2 steps.
- It is unclear whether intended first implementation should be fully general BL or a pragmatic bridge that removes panic path while preserving output contract.
- `bl::multiply` currently has no callsites; intended role in final decomposition pipeline is not encoded.

## Summary of Gap (Current Reality)
- Missing implementation: `src/prime_decomp/bl.rs::decompose`.
- Missing integration: `src/prime_decomp/mod.rs` does not route to BL when `p | (Z_K : Z[theta])`.
- Missing verification assets: tests for index-dividing primes and BL-specific correctness.
