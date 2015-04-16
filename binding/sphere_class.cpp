#include "binding_common.hpp"
#include "Sphere.hpp"

namespace binding {

void register_sphere_class()
{
    register_sphere_class<Sphere>("Sphere");
}

} // namespace binding

