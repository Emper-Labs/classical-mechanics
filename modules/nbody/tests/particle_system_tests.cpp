#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Particle.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

using classical::nbody::NBodySystem;
using classical::nbody::Particle;
using classical::nbody::Vector3;

namespace {

constexpr double kTol = 1e-12;

}  // namespace

TEST(ParticleTest, ValidCreation) {
  const Particle p{Vector3{1.0, 2.0, 3.0}, Vector3{4.0, 5.0, 6.0}, 7.5};
  EXPECT_NEAR(p.position.x, 1.0, kTol);
  EXPECT_NEAR(p.position.y, 2.0, kTol);
  EXPECT_NEAR(p.position.z, 3.0, kTol);
  EXPECT_NEAR(p.velocity.x, 4.0, kTol);
  EXPECT_NEAR(p.velocity.y, 5.0, kTol);
  EXPECT_NEAR(p.velocity.z, 6.0, kTol);
  EXPECT_NEAR(p.mass, 7.5, kTol);
}

TEST(ParticleTest, RejectsNonPositiveMass) {
  EXPECT_THROW(Particle({0, 0, 0}, {0, 0, 0}, 0.0), std::invalid_argument);
  EXPECT_THROW(Particle({0, 0, 0}, {0, 0, 0}, -3.0), std::invalid_argument);
}

TEST(ParticleTest, RejectsInvalidMassValues) {
  EXPECT_THROW(Particle({0, 0, 0}, {0, 0, 0}, std::nan("")), std::invalid_argument);
  EXPECT_THROW(Particle({0, 0, 0}, {0, 0, 0}, std::numeric_limits<double>::infinity()),
               std::invalid_argument);
}

TEST(ParticleTest, RejectsInvalidPositionComponents) {
  EXPECT_THROW(Particle({std::nan(""), 0, 0}, {0, 0, 0}, 1.0),
               std::invalid_argument);
  EXPECT_THROW(Particle({0, std::numeric_limits<double>::infinity(), 0}, {0, 0, 0}, 1.0),
               std::invalid_argument);
}

TEST(ParticleTest, RejectsInvalidVelocityComponents) {
  EXPECT_THROW(Particle({0, 0, 0}, {std::nan(""), 0, 0}, 1.0),
               std::invalid_argument);
  EXPECT_THROW(Particle({0, 0, 0}, {0, -std::numeric_limits<double>::infinity(), 0}, 1.0),
               std::invalid_argument);
}

TEST(ParticleTest, SetMassValidation) {
  Particle p;
  EXPECT_THROW(p.setMass(0.0), std::invalid_argument);
  EXPECT_THROW(p.setMass(-1.0), std::invalid_argument);
  EXPECT_THROW(p.setMass(std::nan("")), std::invalid_argument);
  p.setMass(2.0);
  EXPECT_NEAR(p.mass, 2.0, kTol);
}

TEST(SystemTest, EmptySystemBehavior) {
  NBodySystem sys;
  EXPECT_TRUE(sys.empty());
  EXPECT_EQ(sys.particleCount(), std::size_t{0});
}

TEST(SystemTest, AddAndAccessParticles) {
  NBodySystem sys;
  const std::size_t i0 = sys.addParticle(Particle{{1, 2, 3}, {0, 0, 0}, 2.0});
  const std::size_t i1 = sys.addParticle(Particle{{-1, 0, 0}, {1, 0, 0}, 5.0});
  EXPECT_EQ(i0, std::size_t{0});
  EXPECT_EQ(i1, std::size_t{1});
  EXPECT_EQ(sys.particleCount(), std::size_t{2});

  EXPECT_NEAR(sys.particle(0).mass, 2.0, kTol);
  EXPECT_NEAR(sys.particle(1).position.x, -1.0, kTol);

  // Mutable access.
  sys.particle(1).velocity = Vector3{7, 7, 7};
  EXPECT_NEAR(sys.particle(1).velocity.y, 7.0, kTol);
}

TEST(SystemTest, OutOfRangeAccessThrows) {
  NBodySystem sys;
  sys.addParticle(Particle{});
  EXPECT_EQ(sys.particleCount(), std::size_t{1});
  EXPECT_THROW(sys.particle(5), std::out_of_range);
}