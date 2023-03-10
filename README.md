# Comparison of Barnes-Hut vs brute force

The Barnes-Hut simulation is a divide-and-conquer approximation algorithm that performs 
n-body simulations in O(n log n) time complexity, compared to O(n2) for a brute-force 
approach.

For UBC's PHYS 410 term paper.

## Bench test

Run bench test with (requires [rust](https://www.rust-lang.org/tools/install))

```
cargo bench
```

A report will be generated in `./target/criterion/report/index.html`

## Results

![results.svg](https://raw.githubusercontent.com/lucasvanmol/barnes-hut-bench/master/results.svg)

## References

Barnes, J. and Hut, P., 1986. A hierarchical O (N log N) force-calculation algorithm. nature, 324(6096), pp.446-449.
