# N-body

An independent, CPU-based Newtonian gravitational N-body simulation module.

> **Status:** This is a standalone scientific module. It has no dependency on the
> Emper engine (Emper's Storage and other infrastructure are still under
> development). The public interface is designed to be integrated with Emper
> Storage later.

## Physical model

Each particle is a point mass with a 3D position `r`, velocity `v`, and mass `m`.
The mutual gravitational acceleration obeys Newton's law of universal
gravitation with a softening parameter:

```
a_i = G * sum_{j != i} m_j * (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
```

All quantities use SI units:

| Quantity | Unit            |
| -------- | --------------- |
| position | metres (m)      |
| velocity | metres/second   |
| mass     | kilograms (kg)  |
| time     | seconds (s)     |
| `G`      | `6.67430e-11` m^3 kg^-1 s^-2 |

The gravitational constant is exposed as
`classical::nbody::GRAVITATIONAL_CONSTANT`.

## Gravitational softening

`eps` (metres) softens the 1/r² singularity that would otherwise arise when two
bodies approach coincidence. With `eps = 0` the classical Newtonian field is
recovered exactly. Setting `eps > 0` caps the pairwise acceleration, which is
useful when close encounters / collisions are not modelled explicitly. `eps`
must be finite and non-negative.

## Supported integration algorithms

Both algorithms keep the force calculation separate from the stepping scheme, so
additional integrators can be added later without touching the physics.

- **Semi-implicit Euler (Symplectic Euler):** first-order, symplectic.
  ```
  v += a * dt
  x += v * dt
  ```
- **Velocity Verlet:** second-order, symplectic. Accelerations are evaluated at
  both the beginning and the end of each timestep.
  ```
  a0 = F(x)
  x += v * dt + 0.5 * a0 * dt^2
  a1 = F(x_new)
  v += 0.5 * (a0 + a1) * dt
  ```

A timestep `dt` must be finite and non-negative; `dt == 0` is a no-op.

## Computational complexity

Direct pairwise summation: **O(N²)** per force evaluation, `O(N²)` per step.
There is no Barnes-Hut tree, no GPU code, and no parallel execution yet — this
build is focused on correctness and a clean, testable interface, not peak
performance.

## Building and testing

Requires CMake ≥ 3.16 and a C++17 compiler. GoogleTest is fetched automatically
via CMake `FetchContent` (pinned to a known release). An internet connection is
needed on the first configure.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Running the sample

A two-body circular orbit sample (CPU only, no rendering dependencies):

```bash
cmake --build build --target two_body -j
./build/samples/nbody/two_body
```

It prints per-body positions at regular intervals and the final total energy and
momentum.

## Public API overview

Namespace `classical::nbody`:

- `Vector3` — value-type 3D vector (header-only).
- `Particle` — position, velocity, mass; mass must be finite and `> 0`,
  positions/velocities must be finite.
- `NBodySystem` — owns a mutable `std::vector<Particle>` (no ECS).
- `ForceCalculator` — configurable `G` and `eps`; `computeAccelerations(system)`.
- `integration::semiImplicitEuler(...)` and `integration::velocityVerlet(...)`.

## Limitations

- No Barnes-Hut or hierarchical force evaluation; `O(N²)` only.
- No GPU acceleration or parallelism.
- Direct summation can accumulate round-off with large particle counts.
- Softening must be configured explicitly for collision-dominated problems.
- Integrators are single-step; no adaptive substepping between the chosen `dt`.