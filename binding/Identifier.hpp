#ifndef BINDING_IDENTIFIER_HPP
#define BINDING_IDENTIFIER_HPP

#include <boost/python.hpp>

namespace binding
{
    using namespace boost::python;

    template<typename Timpl_>
    class IdentifierWrapper
    {

    public:
        static void __register_class(char const* name)
        {
            class_<Timpl_>(name)
                .def(init<typename Timpl_::value_type>())
                .add_property("serial", &Timpl_::serial)
                .def(self_ns::repr(self))
                .def(self_ns::str(self))
                .def("__hash__", &Timpl_::serial)
                .def(self == self)
                .def(self != self)
                .def(self > self)
                .def(self < self)
                .def(self >= self)
                .def(self <= self)
                ;
        }
    };

} //namespace binding

#endif /* BINDING_IDENTIFIER_HPP */
