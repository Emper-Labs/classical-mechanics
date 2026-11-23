#pragma once

#include "classical/nbody/Constants.hpp"
#include "classical/nbody/NBodySystem.hpp"
#include "classical/nbody/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace classical {
namespace nbody {

// Computes pairwise gravitational accelerations by direct O(N^2) summation:
//
//   a_i = G * sum_{j != i} m_j * (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
//
// eps is the gravitational softening length, in metres. It removes the 1/r^2
// singularity as two bodies approach coincidence; with eps == 0 the classical
// Newtonian acceleration is recovered exactly. eps must be finite and >= 0.
//
// The pairwise force is evaluated symmetrically (Newton's third law), so the
// total internal force over all particles is zero by construction.
class ForceCalculator {
 public:
  // gravitational_constant defaults to the CODATA value of G; setting it to 0
  // disables all gravitational interaction (useful for tests / free particles).
  // Both arguments must be finite and non-negative.
  explicit ForceCalculator(double gravitational_constant = GRAVITATIONAL_CONSTANT,
                           double softening = 0.0);

  double gravitationalConstant() const noexcept { return g_; }
  double softeningLength() const noexcept { return eps_; }

  // Returns one acceleration vector (m/s^2) per particle, evaluated from the
  // system's current state. Self-interaction is excluded. Returns an empty
  // vector for an empty system.
  std::vector<Vector3> computeAccelerations(const NBodySystem& system) const;

 private:
  double g_;
  double eps_;
};

}  // namespace nbody
}  // namespace classical