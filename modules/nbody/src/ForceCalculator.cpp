#include "classical/nbody/ForceCalculator.hpp"

#include <cmath>
#include <stdexcept>

namespace classical {
namespace nbody {

ForceCalculator::ForceCalculator(double gravitational_constant, double softening)
    : g_(gravitational_constant), eps_(softening) {
  if (!std::isfinite(g_) || g_ < 0.0) {
    throw std::invalid_argument(
        "gravitational_constant must be finite and non-negative.");
  }
  if (!std::isfinite(eps_) || eps_ < 0.0) {
    throw std::invalid_argument("softening must be finite and non-negative.");
  }
}

std::vector<Vector3> ForceCalculator::computeAccelerations(
    const NBodySystem& system) const {
  const std::size_t n = system.particleCount();
  std::vector<Vector3> acceleration(n);

  // Empty systems and a zero gravitational constant produce zero accelerations.
  if (n < 2 || g_ == 0.0) {
    return acceleration;
  }

  const double eps2 = eps_ * eps_;
  const auto& particles = system.particles();

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const Vector3 displacement = particles[j].position - particles[i].position;
      const double r2 = squaredNorm(displacement);
      const double denominator = r2 + eps2;

      // Coincident particles with eps == 0 produce an undefined 1/0 -> skip
      // the pair (softening is recommended to avoid this).
      if (denominator == 0.0) {
        continue;
      }

      // inv3 = (|r_j - r_i|^2 + eps^2)^(-3/2)
      const double inv3 = 1.0 / (denominator * std::sqrt(denominator));
      const double factor = g_ * inv3;

      // Newton's third law: the two per-particle accelerations are evaluated
      // symmetrically so the sum_i m_i * a_i is exactly zero.
      acceleration[i] += (particles[j].mass * displacement) * factor;
      acceleration[j] += (particles[i].mass * (particles[i].position -
                                               particles[j].position)) *
                         factor;
    }
  }

  return acceleration;
}

}  // namespace nbody
}  // namespace classical