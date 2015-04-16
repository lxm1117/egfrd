#include "Particle.hpp"
#include "binding_common.hpp"

namespace binding {

void register_particle_class()
// Parameterize the particle wrapper with the correct Particle class and registers it
{
    ParticleWrapper<Particle>::__register_class("Particle");
}

} // namespace binding
