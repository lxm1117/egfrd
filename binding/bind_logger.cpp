#include <boost/python.hpp>
#include "../Logger.hpp"


namespace binding {


    ////// Registering master function
    boost::python::objects::class_base register_logger_class(char const* name)
    {
        using namespace boost::python;
        typedef Logger impl_type;

        // defining the python class
        return class_<impl_type, boost::noncopyable>(name, no_init)
            .add_property("level", static_cast<impl_type::loglevel(impl_type::*)() const>(&impl_type::level), static_cast<void (impl_type::*)(impl_type::loglevel)>(&impl_type::level))
            .add_property("name", &impl_type::name)
            .add_property("manager", make_function(&impl_type::manager, return_value_policy<return_by_value>()))
            .def("flush", &impl_type::flush)
            .def("get_logger", &impl_type::get_logger, return_value_policy<reference_existing_object>())
            .staticmethod("get_logger")
            .def("stringize_error_level", &impl_type::stringize_error_level)
            .staticmethod("stringize_error_level")
            ;
    }
} // namespace binding
