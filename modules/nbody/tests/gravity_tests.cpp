#include "classical/nbody/Constants.hpp"
#include "classical/nbody/ForceCalculator.hpp"
#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Particle.hpp"
#include "classical/nbody/Vector3.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using classical::nbody::ForceCalculator;
using classical::nbody::GRAVITATIONAL_CONSTANT;
using classical::nbody::NBodySystem;
using classical::nbody::Particle;
using classical::nbody::Vector3;

namespace {

constexpr double kRelTol = 1e-9;
using classical::nbody::norm;
using classical::nbody::normalized;

NBodySystem twoBodySystem(double m1, double m2, double separation) {
  NBodySystem sys;
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, m1});
  sys.addParticle(Particle{{separation, 0, 0}, {0, 0, 0}, m2});
  return sys;
}

}  // namespace

TEST(GravityTest, TwoBodyAccelerationMagnitudeAndDirection) {
  const double source_mass = 1.0;
  const double test_mass = 1.0;
  const double d = 5.0;
  NBodySystem sys = twoBodySystem(source_mass, test_mass, d);
  ForceCalculator fc;  // default G, no softening

  const auto acc = fc.computeAccelerations(sys);
  ASSERT_EQ(acc.size(), std::size_t{2});

  // Test particle at x = d is pulled toward the source at the origin (-x).
  EXPECT_NEAR(acc[1].x, -GRAVITATIONAL_CONSTANT * source_mass / (d * d), kRelTol * GRAVITATIONAL_CONSTANT);
  EXPECT_NEAR(acc[1].y, 0.0, 1e-12);
  EXPECT_NEAR(acc[1].z, 0.0, 1e-12);

  // The source particle is pulled in the opposite (+x) direction.
  EXPECT_NEAR(acc[0].x, GRAVITATIONAL_CONSTANT * test_mass / (d * d), kRelTol * GRAVITATIONAL_CONSTANT);

  // The test particle's acceleration points back toward the origin.
  const Vector3 dir = normalized(acc[1] * -1.0);
  const Vector3 expected_dir{1.0, 0.0, 0.0};
  EXPECT_NEAR(norm(dir - expected_dir), 0.0, 1e-9);
}

TEST(GravityTest, AccelerationDependsOnSourceNotSelfMass) {
  // Acceleration of particle 0 must not depend on its own mass.
  const double d = 4.0;
  NBodySystem sys;
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, 1000.0});
  sys.addParticle(Particle{{d, 0, 0}, {0, 0, 0}, 2.0});
  ForceCalculator fc;

  const auto acc = fc.computeAccelerations(sys);
  // a_0 = G * m_1 / d^2 -- only the OTHER particle's mass appears.
  EXPECT_NEAR(acc[0].x, GRAVITATIONAL_CONSTANT * 2.0 / (d * d), 1e-9);
}

TEST(GravityTest, SelfInteractionExcluded) {
  NBodySystem sys;
  sys.addParticle(Particle{{3, 4, 5}, {0, 0, 0}, 1.0});
  ForceCalculator fc;

  const auto acc = fc.computeAccelerations(sys);
  ASSERT_EQ(acc.size(), std::size_t{1});
  EXPECT_NEAR(acc[0].x, 0.0, 1e-12);
  EXPECT_NEAR(acc[0].y, 0.0, 1e-12);
  EXPECT_NEAR(acc[0].z, 0.0, 1e-12);
}

TEST(GravityTest, EmptySystemProducesNoAccelerations) {
  NBodySystem sys;
  ForceCalculator fc;
  EXPECT_TRUE(fc.computeAccelerations(sys).empty());
}

TEST(GravityTest, ZeroGravitationalConstantProducesNoAcceleration) {
  NBodySystem sys = twoBodySystem(1.0, 1.0, 5.0);
  ForceCalculator fc(0.0, 0.0);
  const auto acc = fc.computeAccelerations(sys);
  ASSERT_EQ(acc.size(), std::size_t{2});
  EXPECT_NEAR(acc[0].x, 0.0, 1e-12);
  EXPECT_NEAR(acc[1].x, 0.0, 1e-12);
}

TEST(GravityTest, SofteningReducesMagnitude) {
  const double m1 = 1.0, m2 = 2.0, d = 5.0, eps = 1.0;
  NBodySystem sys = twoBodySystem(m1, m2, d);

  ForceCalculator newtonian;  // eps = 0
  ForceCalculator softened(GRAVITATIONAL_CONSTANT, eps);

  const auto a_newton = newtonian.computeAccelerations(sys);
  const auto a_soft = softened.computeAccelerations(sys);

  // Softened magnitude < Newtonian magnitude.
  EXPECT_GT(norm(a_newton[1]), norm(a_soft[1]));

  // Matches the analytical softened formula: G*m1*d/(d^2+eps^2)^(3/2).
  const double expected = GRAVITATIONAL_CONSTANT * m1 * d /
                          std::pow(d * d + eps * eps, 1.5);
  EXPECT_NEAR(a_soft[1].x, -expected, 1e-9);
}

TEST(GravityTest, SofteningRemovesSingularityAtCoincidence) {
  NBodySystem sys;
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, 1.0});
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, 1.0});

  ForceCalculator softened(GRAVITATIONAL_CONSTANT, 0.5);
  const auto acc = softened.computeAccelerations(sys);
  ASSERT_EQ(acc.size(), std::size_t{2});
  // Zero separation -> zero net displacement -> zero acceleration, and finite.
  EXPECT_TRUE(std::isfinite(acc[0].x));
  EXPECT_NEAR(acc[0].x, 0.0, 1e-12);
  EXPECT_NEAR(acc[1].y, 0.0, 1e-12);
  EXPECT_NEAR(acc[1].z, 0.0, 1e-12);
}

TEST(GravityTest, TotalInternalForceConserved) {
  // Asymmetric masses; the weighted sum m_i * a_i must be ~ 0.
  NBodySystem sys;
  sys.addParticle(Particle{{0, 0, 0}, {0, 0, 0}, 5.0});
  sys.addParticle(Particle{{2, 1, 0}, {0, 0, 0}, 3.0});
  sys.addParticle(Particle{{-1, -1, 2}, {0, 0, 0}, 7.0});

  ForceCalculator fc(GRAVITATIONAL_CONSTANT, 0.01);
  const auto acc = fc.computeAccelerations(sys);
  ASSERT_EQ(acc.size(), std::size_t{3});

  Vector3 total{0, 0, 0};
  for (std::size_t i = 0; i < acc.size(); ++i) {
    total += acc[i] * sys.particle(i).mass;
  }
  EXPECT_NEAR(total.x, 0.0, 1e-9);
  EXPECT_NEAR(total.y, 0.0, 1e-9);
  EXPECT_NEAR(total.z, 0.0, 1e-9);
}

TEST(GravityTest, RejectsInvalidParameters) {
  EXPECT_THROW(ForceCalculator(-1.0, 0.0), std::invalid_argument);
  EXPECT_THROW(ForceCalculator(std::nan(""), 0.0), std::invalid_argument);
  EXPECT_THROW(ForceCalculator(GRAVITATIONAL_CONSTANT, -0.1), std::invalid_argument);
  EXPECT_THROW(ForceCalculator(GRAVITATIONAL_CONSTANT, std::numeric_limits<double>::infinity()),
               std::invalid_argument);
}