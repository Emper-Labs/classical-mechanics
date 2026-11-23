#pragma once

#include <cmath>

namespace classical {
namespace nbody {

// A simple 3D vector with double-precision components.
// Positions are in metres, velocities in metres/second. This is a plain value
// type: elements are copied by value so ordinary C++ value semantics apply.
struct Vector3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  constexpr Vector3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
      : x(x_), y(y_), z(z_) {}

  Vector3& operator+=(const Vector3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  Vector3& operator-=(const Vector3& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }
  Vector3& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  Vector3& operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }
};

inline Vector3 operator+(Vector3 a, const Vector3& b) { return a += b; }
inline Vector3 operator-(Vector3 a, const Vector3& b) { return a -= b; }
inline Vector3 operator-(const Vector3& v) { return Vector3{-v.x, -v.y, -v.z}; }
inline Vector3 operator*(Vector3 v, double s) { return v *= s; }
inline Vector3 operator*(double s, Vector3 v) { return v *= s; }
inline Vector3 operator/(Vector3 v, double s) { return v /= s; }

inline double dot(const Vector3& a, const Vector3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline double squaredNorm(const Vector3& v) { return dot(v, v); }
inline double norm(const Vector3& v) { return std::sqrt(dot(v, v)); }

// Returns the vector normalized to unit length. The zero vector is returned
// unchanged, so this never produces NaN/Inf for a well-formed input. Callers
// that require an explicit error on a degenerate input should check
// squaredNorm(v) == 0 first.
inline Vector3 normalized(const Vector3& v) {
  const double n = norm(v);
  if (n == 0.0) {
    return Vector3{0.0, 0.0, 0.0};
  }
  return v / n;
}

}  // namespace nbody
}  // namespace classical