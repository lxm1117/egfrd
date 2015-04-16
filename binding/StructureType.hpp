#ifndef BINDING_STRUCTURE_TYPE_HPP
#define BINDING_STRUCTURE_TYPE_HPP

#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/lexical_cast.hpp>
#include "../Model.hpp"

namespace binding {


////// Defining structure with some helper functions.
template<typename Timpl_>
struct StructureTypeExtras
{
    typedef Timpl_ impl_type;

    static std::string const& __getitem__(impl_type* impl, std::string const& key)
    {
        return (*impl)[key];
    }

    static void __setitem__(impl_type* impl, std::string const& key, std::string const& val)
    {
        (*impl)[key] = val;
    }

    static std::string __str__(impl_type* impl)
    {
        return boost::lexical_cast<std::string>(*impl);
    }
};


////// Registering master function
template<typename Timpl_>
static boost::python::objects::class_base register_structure_type_class(
        char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;
    typedef StructureTypeExtras<impl_type> extras_type;

    // defining the python class
    return class_<impl_type, boost::shared_ptr<impl_type> >(name, init<>())
        .add_property("id",
                make_function((typename impl_type::identifier_type const&(impl_type::*)() const)&impl_type::id,
                    return_value_policy<copy_const_reference>()))
        .add_property("model",
                make_function((Model*(impl_type::*)() const)&impl_type::model,
                    return_value_policy<reference_existing_object>()))
        .def("__str__", &extras_type::__str__)
        .def("__getitem__", &extras_type::__getitem__,
                return_value_policy<copy_const_reference>())
        .def("__setitem__", &extras_type::__setitem__)
        ;
}

} // namespace binding

#endif /* BINDING_STRUCTURE_TYPE_HPP */
