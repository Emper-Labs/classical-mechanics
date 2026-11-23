#include "classical/nbody/Constants.hpp"
#include "classical/nbody/ForceCalculator.hpp"
#include "classical/nbody/Integrator.hpp"
#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Particle.hpp"
#include "classical/nbody/Vector3.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using classical::nbody::ForceCalculator;
using classical::nbody::GRAVITATIONAL_CONSTANT;
using classical::nbody::NBodySystem;
using classical::nbody::Particle;
using classical::nbody::norm;
using classical::nbody::squaredNorm;
using classical::nbody::integration::velocityVerlet;

namespace {

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

}  // namespace

int main() {
  // Two equal masses on a circular orbit about their common centre of mass.
  const double mass = 1e12;   // kg
  const double separation = 10.0;  // m
  const double omega = std::sqrt(GRAVITATIONAL_CONSTANT * (2.0 * mass) /
                                 (separation * separation * separation));

  NBodySystem system;
  system.addParticle(Particle{{-separation / 2.0, 0, 0},
                              {0, -omega * separation / 2.0, 0}, mass});
  system.addParticle(Particle{{separation / 2.0, 0, 0},
                              {0, omega * separation / 2.0, 0}, mass});

  ForceCalculator force(GRAVITATIONAL_CONSTANT, /*softening=*/0.0);

  const double dt = 0.05;      // seconds
  const int steps = 400;
  const int printEvery = 50;

  std::printf("Two-body circular orbit (Velocity Verlet, G = %.4e)\n",
              GRAVITATIONAL_CONSTANT);
  std::printf("%6s %11s %11s %11s %11s\n", "t(s)", "x1", "y1", "x2", "y2");

  for (int k = 0; k <= steps; ++k) {
    if (k % printEvery == 0) {
      const auto& b1 = system.particle(0);
      const auto& b2 = system.particle(1);
      std::printf("%6.2f %11.6f %11.6f %11.6f %11.6f\n", k * dt,
                  b1.position.x, b1.position.y, b2.position.x, b2.position.y);
    }
    if (k < steps) {
      velocityVerlet(system, dt, force);
    }
  }

  const double e = totalEnergy(system, GRAVITATIONAL_CONSTANT);
  const auto& b1 = system.particle(0);
  const auto& b2 = system.particle(1);
  const double p_x = mass * (b1.velocity.x + b2.velocity.x);
  const double p_y = mass * (b1.velocity.y + b2.velocity.y);

  std::printf("\nFinal energy     : %.6e J\n", e);
  std::printf("Total momentum   : (%.3e, %.3e, 0) kg m/s\n", p_x, p_y);
  return EXIT_SUCCESS;
}