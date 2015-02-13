#ifndef BIND_LOGGER_HPP
#define BIND_LOGGER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_logger_class(char const* name);

} // namespace binding

#endif /* BIND_LOGGER_HPP */
