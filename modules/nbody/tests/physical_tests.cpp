#include "classical/nbody/Constants.hpp"
#include "classical/nbody/ForceCalculator.hpp"
#include "classical/nbody/Integrator.hpp"
#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Particle.hpp"
#include "classical/nbody/Vector3.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

using classical::nbody::ForceCalculator;
using classical::nbody::GRAVITATIONAL_CONSTANT;
using classical::nbody::NBodySystem;
using classical::nbody::Particle;
using classical::nbody::Vector3;
using classical::nbody::integration::velocityVerlet;
using classical::nbody::norm;
using classical::nbody::squaredNorm;

namespace classical {
namespace nbody {
inline Vector3 cross(const Vector3& a, const Vector3& b) {
  return Vector3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x};
}
}  // namespace nbody
}  // namespace classical

namespace {

// Total mechanical energy of the system (kinetic + potential over all pairs).
double totalEnergy(const NBodySystem& sys, double g) {
  double kinetic = 0.0;
  for (const auto& p : sys.particles()) {
    kinetic += 0.5 * p.mass * squaredNorm(p.velocity);
  }
  double potential = 0.0;
  const auto& ps = sys.particles();
  for (std::size_t i = 0; i < ps.size(); ++i) {
    for (std::size_t j = i + 1; j < ps.size(); ++j) {
      potential -= g * ps[i].mass * ps[j].mass /
                   norm(ps[i].position - ps[j].position);
    }
  }
  return kinetic + potential;
}

Vector3 totalMomentum(const NBodySystem& sys) {
  Vector3 p{0, 0, 0};
  for (const auto& body : sys.particles()) {
    p += body.velocity * body.mass;
  }
  return p;
}

Vector3 totalAngularMomentum(const NBodySystem& sys) {
  Vector3 l{0, 0, 0};
  for (const auto& body : sys.particles()) {
    l += classical::nbody::cross(body.position, body.velocity) * body.mass;
  }
  return l;
}

// Builds a two-body circular orbit in the centre-of-mass frame with equal
// masses m, separated by distance r.
NBodySystem circularTwoBody(double m, double r) {
  const double omega =
      std::sqrt(GRAVITATIONAL_CONSTANT * (2.0 * m) / (r * r * r));
  NBodySystem sys;
  sys.addParticle(Particle{{-r / 2.0, 0, 0}, {0, -omega * r / 2.0, 0}, m});
  sys.addParticle(Particle{{r / 2.0, 0, 0}, {0, omega * r / 2.0, 0}, m});
  return sys;
}

}  // namespace

TEST(PhysicalValidation, CircularOrbitRunsStably) {
  constexpr double m = 1e12;   // kg
  constexpr double r = 10.0;   // m
  constexpr double dt = 0.01;  // s
  constexpr int steps = 2000;

  NBodySystem sys = circularTwoBody(m, r);
  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.0);

  const double e0 = totalEnergy(sys, GRAVITATIONAL_CONSTANT);
  const Vector3 p0 = totalMomentum(sys);
  const Vector3 l0 = totalAngularMomentum(sys);

  double e_min = e0, e_max = e0;
  double p_deviations[3] = {0.0, 0.0, 0.0};
  double l_deviations[3] = {0.0, 0.0, 0.0};

  for (int k = 0; k < steps; ++k) {
    velocityVerlet(sys, dt, fc);

    const double e = totalEnergy(sys, GRAVITATIONAL_CONSTANT);
    e_min = std::min(e_min, e);
    e_max = std::max(e_max, e);

    const Vector3 p = totalMomentum(sys);
    const Vector3 l = totalAngularMomentum(sys);
    p_deviations[0] = std::max(p_deviations[0], std::abs(p.x - p0.x));
    p_deviations[1] = std::max(p_deviations[1], std::abs(p.y - p0.y));
    p_deviations[2] = std::max(p_deviations[2], std::abs(p.z - p0.z));
    l_deviations[0] = std::max(l_deviations[0], std::abs(l.x - l0.x));
    l_deviations[1] = std::max(l_deviations[1], std::abs(l.y - l0.y));
    l_deviations[2] = std::max(l_deviations[2], std::abs(l.z - l0.z));
  }

  // Bounded energy error over the interval.
  const double energy_span = (e_max - e_min) / std::abs(e0);
  EXPECT_LT(energy_span, 1e-3);

  // Approximate momentum conservation (should be near machine precision).
  const double p_scale = m;  // mass ~ kg m/s magnitude baseline
  for (int c = 0; c < 3; ++c) {
    EXPECT_LT(p_deviations[c], 1e-6 * p_scale);
  }

  // Approximate angular-momentum conservation.
  const double l_magnitude = norm(l0);
  for (int c = 0; c < 3; ++c) {
    EXPECT_LT(l_deviations[c], 1e-6 * l_magnitude);
  }

  // Momentum was set up to be ~zero in the COM frame.
  EXPECT_LT(norm(p0), 1e-6 * m);
}