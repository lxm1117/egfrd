#include "binding_common.hpp"
#include "BDSimulator.hpp"

namespace binding {

void register_bd_simulator_classes()
{
    register_bd_simulator_class<BDSimulator>("_BDSimulator");
}

} // namespace binding
