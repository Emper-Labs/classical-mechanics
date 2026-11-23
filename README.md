# classical-mechanics

Classical-mechanics modules for Emper (`classical::`).

This repository hosts independent, CPU-based scientific modules that are
intended to be integrated with the Emper engine later. See
[`modules/nbody/README.md`](modules/nbody/README.md) for the first module: an
N-body gravitational simulator.

## Layout

```text
CMakeLists.txt
modules/
  nbody/        # N-body gravitational simulation (module + unit tests)
samples/
  nbody/        # two-body circular orbit sample
```

## Build & test

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```
