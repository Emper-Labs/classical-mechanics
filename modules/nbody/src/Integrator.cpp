#include "classical/nbody/Integrator.hpp"

#include <cmath>
#include <stdexcept>

namespace classical {
namespace nbody {
namespace integration {

namespace {
void validateDt(double dt) {
  if (!std::isfinite(dt)) {
    throw std::invalid_argument("Timestep dt must be finite (non-NaN, non-infinite).");
  }
  if (dt < 0.0) {
    throw std::invalid_argument("Timestep dt must be non-negative.");
  }
}
}  // namespace

void semiImplicitEuler(NBodySystem& system, double dt,
                       const ForceCalculator& force) {
  validateDt(dt);
  const std::size_t n = system.particleCount();
  if (n == 0 || dt == 0.0) {
    return;
  }

  const auto acceleration = force.computeAccelerations(system);
  for (std::size_t i = 0; i < n; ++i) {
    Particle& p = system.particle(i);
    p.velocity += acceleration[i] * dt;
    p.position += p.velocity * dt;
  }
}

void velocityVerlet(NBodySystem& system, double dt,
                    const ForceCalculator& force) {
  validateDt(dt);
  const std::size_t n = system.particleCount();
  if (n == 0 || dt == 0.0) {
    return;
  }

  // Beginning-of-step acceleration.
  const auto a_start = force.computeAccelerations(system);

  const double half_dt2 = 0.5 * dt * dt;
  for (std::size_t i = 0; i < n; ++i) {
    Particle& p = system.particle(i);
    p.position += p.velocity * dt + a_start[i] * half_dt2;
  }

  // End-of-step acceleration evaluated at the new positions.
  const auto a_end = force.computeAccelerations(system);

  const double half_dt = 0.5 * dt;
  for (std::size_t i = 0; i < n; ++i) {
    system.particle(i).velocity += (a_start[i] + a_end[i]) * half_dt;
  }
}

}  // namespace integration
}  // namespace nbody
}  // namespace classical