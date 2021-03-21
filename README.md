# rust-number-theory [![Rust](https://github.com/koba-e964/rust-number-theory/actions/workflows/rust.yml/badge.svg)](https://github.com/koba-e964/rust-number-theory/actions/workflows/rust.yml)
Rust implementation of algorithms in number theory.
Implemented algorithms include:

- HNF
- Discriminant / resultant of polyomials


## Examples
Take a look at the file `data/input-discriminant.yml`.

```
# Find the discriminant of 2x^3 + x^2 - 2x + 3.
# To feed this file to the executable, run `cargo run data/input-discriminant.yml'.
# The output should be {"discriminant":"-1132"}.
input:
  polynomials:
    -
      - '3'
      - '-2'
      - '1'
      - '2'

to_find:
  - 'discriminant'
```

If you feed this file to the executable, the discriminant of 2x^3 + x^2 - 2x + 3 will be found.
```
$ cargo run data/input-discriminant.yml
    Finished dev [unoptimized + debuginfo] target(s) in 0.37s
     Running `target/debug/rust-number-theory data/input-discriminant.yml`
{
  "discriminant": "-1132"
}
```