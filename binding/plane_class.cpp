#include "binding_common.hpp"
#include "Plane.hpp"

namespace binding {

void register_plane_class()
{
    register_plane_class<Plane>("Plane");
}

} // namespace binding
