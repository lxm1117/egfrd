#include <boost/python.hpp>
#include "../Logger.hpp"

namespace binding {

template<typename Timpl>
static void log_appender_call(Timpl& self, Logger::loglevel lv, char const* name, PyObject* _chunks)
{
    using namespace boost::python;

    handle<> iter(allow_null(PyObject_GetIter(_chunks)));
    if (!iter)
    {
        PyErr_SetString(PyExc_TypeError, "chunks must be iterable");
        throw_error_already_set();
    }

    std::vector<char const*> chunks;
    for (;;)
    {
        handle<> item(allow_null(PyIter_Next(iter.get())));
        if (!item) break;
        chunks.push_back(extract<char const*>(item.get()));
    }

    chunks.push_back(0);
    self(lv, name, &chunks.front());
}


////// Registering master function
boost::python::objects::class_base register_log_appender_class(char const* name)
{
    using namespace boost::python;
    typedef LogAppender impl_type;

    // defining the python class
    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .def("flush", &impl_type::flush)
        .def("__call__", &log_appender_call<impl_type>)
        ;
}

} // namespace binding
