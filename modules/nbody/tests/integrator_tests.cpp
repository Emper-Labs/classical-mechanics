#include "classical/nbody/Constants.hpp"
#include "classical/nbody/ForceCalculator.hpp"
#include "classical/nbody/Integrator.hpp"
#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Particle.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

using classical::nbody::ForceCalculator;
using classical::nbody::GRAVITATIONAL_CONSTANT;
using classical::nbody::NBodySystem;
using classical::nbody::Particle;
using classical::nbody::Vector3;
using classical::nbody::integration::semiImplicitEuler;
using classical::nbody::integration::velocityVerlet;

namespace {

constexpr double kTol = 1e-9;

NBodySystem twoBodySystem() {
  NBodySystem sys;
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, 1.0});
  sys.addParticle(Particle{{10.0, 0, 0}, {0, 1.0, 0}, 2.0});
  return sys;
}

}  // namespace

TEST(SemiImplicitEulerTest, UpdatesVelocityThenPosition) {
  NBodySystem sys = twoBodySystem();
  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.0);
  const double dt = 0.1;

  const auto a0 = fc.computeAccelerations(sys);
  const Vector3 v0_0 = sys.particle(0).velocity;
  const Vector3 p0_0 = sys.particle(0).position;

  semiImplicitEuler(sys, dt, fc);

  // Velocity updated first by a*dt.
  const Vector3 expect_v = v0_0 + a0[0] * dt;
  EXPECT_NEAR(sys.particle(0).velocity.x, expect_v.x, kTol);
  EXPECT_NEAR(sys.particle(0).velocity.y, expect_v.y, kTol);

  // Position updated using the NEW velocity.
  const Vector3 expect_p = p0_0 + expect_v * dt;
  EXPECT_NEAR(sys.particle(0).position.x, expect_p.x, kTol);
  EXPECT_NEAR(sys.particle(0).position.y, expect_p.y, kTol);
}

TEST(SemiImplicitEulerTest, FreeParticleMovesAtConstantVelocity) {
  NBodySystem sys;
  sys.addParticle(Particle{{1.0, 2.0, 3.0}, {0.5, -1.0, 2.0}, 1.0});
  ForceCalculator fc(0.0, 0.0);  // no gravity
  const double dt = 2.0;

  semiImplicitEuler(sys, dt, fc);
  EXPECT_NEAR(sys.particle(0).position.x, 1.0 + 0.5 * dt, kTol);
  EXPECT_NEAR(sys.particle(0).position.y, 2.0 - 1.0 * dt, kTol);
  EXPECT_NEAR(sys.particle(0).position.z, 3.0 + 2.0 * dt, kTol);
  // Velocity unchanged.
  EXPECT_NEAR(sys.particle(0).velocity.z, 2.0, kTol);
}

TEST(VelocityVerletTest, FreeParticleMovesAtConstantVelocity) {
  NBodySystem sys;
  sys.addParticle(Particle{{0.0, 0.0, 0.0}, {2.0, 3.0, 4.0}, 1.0});
  ForceCalculator fc(0.0, 0.0);
  const double dt = 0.5;

  velocityVerlet(sys, dt, fc);
  EXPECT_NEAR(sys.particle(0).position.x, 2.0 * dt, kTol);
  EXPECT_NEAR(sys.particle(0).position.y, 3.0 * dt, kTol);
  EXPECT_NEAR(sys.particle(0).position.z, 4.0 * dt, kTol);
  EXPECT_NEAR(sys.particle(0).velocity.x, 2.0, kTol);
}

TEST(VelocityVerletTest, SingleStepMatchesFormulaExactly) {
  NBodySystem sys = twoBodySystem();
  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.0);
  const double dt = 0.1;

  const auto a0 = fc.computeAccelerations(sys);
  const Vector3 p0 = sys.particle(0).position;
  const Vector3 v0 = sys.particle(0).velocity;

  velocityVerlet(sys, dt, fc);

  // Position: x += v*dt + 0.5*a0*dt^2  (exact for one step).
  const Vector3 expect_p = p0 + v0 * dt + a0[0] * (0.5 * dt * dt);
  EXPECT_NEAR(sys.particle(0).position.x, expect_p.x, kTol);
  EXPECT_NEAR(sys.particle(0).position.y, expect_p.y, kTol);

  // Velocity: v += 0.5*(a0 + a1)*dt where a1 is evaluated at new positions.
  const auto a1 = fc.computeAccelerations(sys);
  const Vector3 expect_v = v0 + (a0[0] + a1[0]) * (0.5 * dt);
  EXPECT_NEAR(sys.particle(0).velocity.x, expect_v.x, kTol);
  EXPECT_NEAR(sys.particle(0).velocity.y, expect_v.y, kTol);
}

TEST(ZeroTimestepTest, LeavesSystemUnchanged) {
  NBodySystem sys = twoBodySystem();
  const NBodySystem copy = sys;
  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.0);

  semiImplicitEuler(sys, 0.0, fc);
  velocityVerlet(sys, 0.0, fc);

  EXPECT_EQ(sys.particleCount(), copy.particleCount());
  EXPECT_NEAR(sys.particle(0).position.x, copy.particle(0).position.x, kTol);
  EXPECT_NEAR(sys.particle(1).velocity.y, copy.particle(1).velocity.y, kTol);
}

TEST(InvalidTimestepTest, RejectsNegativeNanInfinite) {
  NBodySystem sys = twoBodySystem();
  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.0);

  EXPECT_THROW(semiImplicitEuler(sys, -1.0, fc), std::invalid_argument);
  EXPECT_THROW(semiImplicitEuler(sys, std::nan(""), fc), std::invalid_argument);
  EXPECT_THROW(semiImplicitEuler(sys, std::numeric_limits<double>::infinity(), fc),
               std::invalid_argument);

  EXPECT_THROW(velocityVerlet(sys, -1.0, fc), std::invalid_argument);
  EXPECT_THROW(velocityVerlet(sys, std::nan(""), fc), std::invalid_argument);
  EXPECT_THROW(velocityVerlet(sys, -std::numeric_limits<double>::infinity(), fc),
               std::invalid_argument);
}