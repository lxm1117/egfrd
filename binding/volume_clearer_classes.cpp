#include <boost/python.hpp>
#include "binding_common.hpp"
#include "volume_clearer_converter.hpp"

namespace binding {

void register_volume_clearer_classes()
{
    register_volume_clearer_converter<VolumeClearer>();
}

} // namespace binding
