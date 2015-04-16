#include "binding_common.hpp"
#include "Box.hpp"

namespace binding {

void register_box_class()
{
    register_box_class<Box>("Box");
}

} // namespace binding
