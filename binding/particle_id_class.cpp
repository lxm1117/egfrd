#include "Identifier.hpp"
#include "SerialIDGenerator.hpp"
#include "binding_common.hpp"
#include "peer/utils.hpp"

namespace binding {

void register_particle_id_class()
{
    IdentifierWrapper<ParticleID>::__register_class("ParticleID");
    register_serial_id_generator_class<ParticleID>("ParticleIDGenerator");
}

} // namespace binding
