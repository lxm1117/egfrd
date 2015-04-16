#include "binding_common.hpp"
#include "EGFRDSimulator.hpp"

namespace binding {

void register_egfrd_simulator_classes()
{
    register_egfrd_simulator_class<EGFRDSimulator>("_EGFRDSimulator");
}

} // namespace binding
