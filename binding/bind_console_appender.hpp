#ifndef BIND_CONSOLE_APPENDER_HPP
#define BIND_CONSOLE_APPENDER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_console_appender_class(char const* name);

} // namespace binding

#endif /* BIND_CONSOLE_APPENDER_HPP */
