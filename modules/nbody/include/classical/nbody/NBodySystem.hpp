#pragma once

#include "classical/nbody/Particle.hpp"

#include <cstddef>
#include <vector>

namespace classical {
namespace nbody {

// Owns a mutable collection of particles, indexed by insertion order.
// This deliberately stays a simple container (std::vector backing); no ECS,
// generic storage framework, or registry is used.
class NBodySystem {
 public:
  NBodySystem() = default;

  // Appends a particle and returns its index (0-based).
  std::size_t addParticle(const Particle& p) {
    particles_.push_back(p);
    return particles_.size() - 1;
  }

  std::size_t particleCount() const noexcept { return particles_.size(); }
  bool empty() const noexcept { return particles_.empty(); }

  // Indexed access; throws std::out_of_range for invalid indices.
  const Particle& particle(std::size_t index) const { return particles_.at(index); }
  Particle& particle(std::size_t index) { return particles_.at(index); }

  const std::vector<Particle>& particles() const noexcept { return particles_; }

 private:
  std::vector<Particle> particles_;
};

}  // namespace nbody
}  // namespace classical