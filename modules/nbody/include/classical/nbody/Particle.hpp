#pragma once

#include "classical/nbody/Vector3.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace classical {
namespace nbody {

// A single point mass. All quantities use SI units:
//   position -> metres
//   velocity -> metres/second
//   mass     -> kilograms
//
// Mass must be finite and strictly positive. Position and velocity must
// contain only finite components. These invariants are validated on
// construction and when mass is changed via setMass().
struct Particle {
  Vector3 position;    // metres
  Vector3 velocity;    // metres/second
  double mass = 1.0;   // kilograms

  Particle() = default;

  Particle(const Vector3& p, const Vector3& v, double m)
      : position(p), velocity(v) {
    validateVector(position, "position");
    validateVector(velocity, "velocity");
    setMass(m);
  }

  void setMass(double m) {
    if (!std::isfinite(m)) {
      throw std::invalid_argument(
          "Particle mass must be finite (non-NaN, non-infinite).");
    }
    if (m <= 0.0) {
      throw std::invalid_argument("Particle mass must be strictly positive.");
    }
    mass = m;
  }

 private:
  static void validateVector(const Vector3& v, const char* name) {
    if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) {
      throw std::invalid_argument(
          std::string("Particle ") + name +
          " must contain only finite components.");
    }
  }
};

}  // namespace nbody
}  // namespace classical