#include "ParticleModel.hpp"
#include "binding_common.hpp"

namespace binding {

void register_particle_model_class()
{
    register_particle_model_class<ParticleModel>("ParticleModel");
}

} // namesapce binding
