#include "MultiParticleContainer.hpp"
#include "binding_common.hpp"

namespace binding {

void register_multi_particle_container_class()
{
    register_multi_particle_container_class<MultiParticleContainer, World>("MultiParticleContainer");
}

} // namespace binding
