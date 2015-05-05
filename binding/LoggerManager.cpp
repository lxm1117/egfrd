#include <boost/python.hpp>
#include "peer/converters/sequence.hpp"
#include "../Logger.hpp"

namespace binding {

    // defining the loglevel
    boost::python::objects::enum_base register_logger_level_enum(char const* name)
    {
        using namespace boost::python;
        return enum_<Logger::loglevel>(name)
            .value("OFF", Logger::loglevel::L_OFF)
            .value("DEBUG", Logger::loglevel::L_DEBUG)
            .value("INFO", Logger::loglevel::L_INFO)
            .value("WARNING", Logger::loglevel::L_WARNING)
            .value("ERROR", Logger::loglevel::L_ERROR)
            .value("FATAL", Logger::loglevel::L_FATAL)
            ;
    }


    ////// Registering master function
    boost::python::objects::class_base register_logger_manager_class(char const* name)
    {
        using namespace boost::python;
        typedef LoggerManager impl_type;

        peer::converters::register_range_to_tuple_converter<std::vector<boost::shared_ptr<LogAppender> > >();

        // defining the python class
        return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
            .add_property("level", static_cast<Logger::loglevel(impl_type::*)() const>(&impl_type::level), static_cast<void(impl_type::*)(Logger::loglevel)>(&impl_type::level))
            .add_property("name", &impl_type::name)
            .add_property("appenders", make_function(&impl_type::appenders, return_value_policy<return_by_value>()))
            .def("add_appender", &impl_type::add_appender)
            .def("register_logger_manager", &impl_type::register_logger_manager)
            .staticmethod("register_logger_manager")
            .def("get_logger_manager", &impl_type::get_logger_manager)
            .staticmethod("get_logger_manager")
            ;
    }

} // namespace binding
