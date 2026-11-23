#include "classical/nbody/Vector3.hpp"

#include <gtest/gtest.h>

#include <cmath>

using classical::nbody::Vector3;
using classical::nbody::dot;
using classical::nbody::norm;
using classical::nbody::normalized;
using classical::nbody::squaredNorm;

namespace {

constexpr double kTol = 1e-12;

}  // namespace

TEST(Vector3Test, Addition) {
  const Vector3 a{1.0, 2.0, 3.0};
  const Vector3 b{4.0, 5.0, 6.0};
  const Vector3 c = a + b;
  EXPECT_NEAR(c.x, 5.0, kTol);
  EXPECT_NEAR(c.y, 7.0, kTol);
  EXPECT_NEAR(c.z, 9.0, kTol);
}

TEST(Vector3Test, Subtraction) {
  const Vector3 a{7.0, 5.0, 3.0};
  const Vector3 b{1.0, 2.0, 3.0};
  const Vector3 c = a - b;
  EXPECT_NEAR(c.x, 6.0, kTol);
  EXPECT_NEAR(c.y, 3.0, kTol);
  EXPECT_NEAR(c.z, 0.0, kTol);
}

TEST(Vector3Test, UnaryMinus) {
  const Vector3 v{1.0, -2.0, 3.0};
  const Vector3 n = -v;
  EXPECT_NEAR(n.x, -1.0, kTol);
  EXPECT_NEAR(n.y, 2.0, kTol);
  EXPECT_NEAR(n.z, -3.0, kTol);
}

TEST(Vector3Test, ScalarMultiplication) {
  const Vector3 v{2.0, -1.0, 3.0};
  const Vector3 a = v * 2.5;
  const Vector3 b = 2.5 * v;
  EXPECT_NEAR(a.x, 5.0, kTol);
  EXPECT_NEAR(a.y, -2.5, kTol);
  EXPECT_NEAR(a.z, 7.5, kTol);
  EXPECT_NEAR(b.x, 5.0, kTol);
  EXPECT_NEAR(b.y, -2.5, kTol);
  EXPECT_NEAR(b.z, 7.5, kTol);
}

TEST(Vector3Test, ScalarDivision) {
  const Vector3 v{4.0, 6.0, 8.0};
  const Vector3 d = v / 2.0;
  EXPECT_NEAR(d.x, 2.0, kTol);
  EXPECT_NEAR(d.y, 3.0, kTol);
  EXPECT_NEAR(d.z, 4.0, kTol);
}

TEST(Vector3Test, DotProduct) {
  EXPECT_NEAR(dot(Vector3{1.0, 0.0, 0.0}, Vector3{0.0, 1.0, 0.0}), 0.0, kTol);
  EXPECT_NEAR(dot(Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0}), 1.0, kTol);
  EXPECT_NEAR(dot(Vector3{2.0, 3.0, 4.0}, Vector3{3.0, 4.0, 5.0}), 38.0, kTol);
  EXPECT_NEAR(dot(Vector3{1.0, 1.0, 1.0}, Vector3{1.0, -1.0, 1.0}), 1.0, kTol);
}

TEST(Vector3Test, SquaredMagnitudeAndMagnitude) {
  const Vector3 v{3.0, 4.0, 0.0};
  EXPECT_NEAR(squaredNorm(v), 25.0, kTol);
  EXPECT_NEAR(norm(v), 5.0, kTol);
  EXPECT_NEAR(norm(Vector3{0.0, 0.0, 0.0}), 0.0, kTol);
}

TEST(Vector3Test, Normalization) {
  const Vector3 v{3.0, 4.0, 0.0};
  const Vector3 u = normalized(v);
  EXPECT_NEAR(u.x, 0.6, kTol);
  EXPECT_NEAR(u.y, 0.8, kTol);
  EXPECT_NEAR(u.z, 0.0, kTol);
  EXPECT_NEAR(norm(u), 1.0, kTol);
}

TEST(Vector3Test, NormalizationOfZeroVectorIsZero) {
  const Vector3 u = normalized(Vector3{0.0, 0.0, 0.0});
  EXPECT_NEAR(u.x, 0.0, kTol);
  EXPECT_NEAR(u.y, 0.0, kTol);
  EXPECT_NEAR(u.z, 0.0, kTol);
  EXPECT_TRUE(std::isfinite(u.x));
  EXPECT_TRUE(std::isfinite(u.y));
  EXPECT_TRUE(std::isfinite(u.z));
}