#include "Cylinder.hpp"
#include "binding_common.hpp"

namespace binding {

void register_cylinder_class()
{
    register_cylinder_class<Cylinder>("Cylinder");
}

} // namespace binding
