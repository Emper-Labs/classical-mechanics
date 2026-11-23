#pragma once

#include "classical/nbody/ForceCalculator.hpp"
#include "classical/nbody/NBodySystem.hpp"

namespace classical {
namespace nbody {
namespace integration {

// Symplectic (semi-implicit) Euler step:
//   v += a * dt
//   x += v * dt
// First-order accurate and symplectic, so energy error does not grow
// monotonically over long oscillatory runs.
void semiImplicitEuler(NBodySystem& system, double dt,
                       const ForceCalculator& force);

// Velocity Verlet step:
//   a0 = F(x)
//   x += v * dt + 0.5 * a0 * dt^2
//   a1 = F(x_new)
//   v += 0.5 * (a0 + a1) * dt
// Second-order accurate and symplectic. Accelerations are evaluated at both
// the beginning and the end of each timestep.
void velocityVerlet(NBodySystem& system, double dt,
                    const ForceCalculator& force);

// dt must be finite and non-negative. dt == 0 leaves the system unchanged.
// A negative, NaN, or infinite dt throws std::invalid_argument.
//
// Force evaluation is kept out of the integrators (supplied as a
// ForceCalculator) so additional integrators can be added later without
// changing the physics.

}  // namespace integration
}  // namespace nbody
}  // namespace classical