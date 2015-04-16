#include "binding_common.hpp"
#include "ParticleSimulator.hpp"

namespace binding {

void register_particle_simulator_classes()
{
    register_particle_simulator_class<ParticleSimulator>("_ParticleSimulator");
}

} // namespace binding
