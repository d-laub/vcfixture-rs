# Bulk generation measurements (issue #22)

Harness: `cargo run --release --example bulk_bench --features bulk`
Machine: 8 cores, 4.18.0-553.36.1.el8_10.x86_64

## Baseline — v0.4.0 (`9e64a0c`), before any change

```
workers=8
 samples   records        cells       secs       s/cell peakRSS_MB
     500      5000      5000000      1.443     2.886e-7       27.2
     500     20000     20000000      5.606     2.803e-7       45.4
    2000      5000     20000000      5.511     2.755e-7       47.9
    2000     20000     80000000     21.746     2.718e-7       78.5
    8000      5000     80000000     21.316     2.664e-7       83.9
    8000     20000    320000000     87.470     2.733e-7      200.9
```

## After the change

(filled in by Task 8)
