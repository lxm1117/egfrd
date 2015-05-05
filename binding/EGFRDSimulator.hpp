#ifndef EGFRD_SIMULATOR_HPP
#define EGFRD_SIMULATOR_HPP

#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include "peer/util/to_native_converter.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/util/shared_const_ptr.hpp"

namespace binding {

//// Converters. These convert C++ types to something that Python can handle.

//
template<typename Timpl_>
struct shell_variant_converter
{
    typedef Timpl_ impl_type;
    typedef typename impl_type::shell_variant_type native_type;
    typedef typename impl_type::spherical_shell_type spherical_shell_type;
    typedef typename impl_type::cylindrical_shell_type cylindrical_shell_type;

    struct visitor: boost::static_visitor<boost::python::object>
    {
        boost::python::object operator()(boost::none_t const&) const
        {
            return boost::python::object();
        }

        boost::python::object operator()(spherical_shell_type const& shell) const
        {
            return boost::python::object(shell);
        }

        boost::python::object operator()(cylindrical_shell_type const& shell) const
        {
            return boost::python::object(shell);
        }
    };

    static PyObject* convert(native_type const& val)
    {
        return boost::python::incref(boost::apply_visitor(visitor(), val).ptr());
    }

    // Also need to register the converter.
    static void __register()
    {
        boost::python::to_python_converter<native_type, shell_variant_converter>();
    }
};


////// Registering master function
template<typename Timpl>
void register_egfrd_simulator_class(char const* name)
{
    using namespace boost::python;
    // typedefs
    typedef Timpl impl_type;
    typedef std::pair<typename impl_type::shell_id_type, typename impl_type::shell_variant_type> get_shell_result_type;
    enum_<typename impl_type::domain_kind>("DomainKind")
        .value("NONE", impl_type::domain_kind::NONE)
        .value("SPHERICAL_SINGLE", impl_type::domain_kind::SPHERICAL_SINGLE)
        .value("CYLINDRICAL_SINGLE", impl_type::domain_kind::CYLINDRICAL_SINGLE)
        .value("SPHERICAL_PAIR", impl_type::domain_kind::SPHERICAL_PAIR)
        .value("CYLINDRICAL_PAIR", impl_type::domain_kind::CYLINDRICAL_PAIR)
        .value("MULTI", impl_type::domain_kind::MULTI)
        ;
    enum_<typename impl_type::single_event_kind>("SingleEventKind")
        .value("REACTION", impl_type::single_event_kind::SINGLE_EVENT_REACTION)
        .value("ESCAPE", impl_type::single_event_kind::SINGLE_EVENT_ESCAPE)
        ;
    enum_<typename impl_type::pair_event_kind>("PairEventKind")
        .value("SINGLE_REACTION_0", impl_type::pair_event_kind::PAIR_EVENT_SINGLE_REACTION_0)
        .value("SINGLE_REACTION_1", impl_type::pair_event_kind::PAIR_EVENT_SINGLE_REACTION_1)
        .value("COM_ESCAPE", impl_type::pair_event_kind::PAIR_EVENT_COM_ESCAPE)
        .value("IV_UNDETERMINED", impl_type::pair_event_kind::PAIR_EVENT_IV_UNDETERMINED)
        .value("IV_ESCAPE", impl_type::pair_event_kind::PAIR_EVENT_IV_ESCAPE)
        .value("IV_REACTION", impl_type::pair_event_kind::PAIR_EVENT_IV_REACTION)
        ;
    enum_<typename impl_type::multi_type::event_kind>("MultiEventKind")
        .value("NONE", impl_type::multi_type::event_kind::NONE)
        .value("ESCAPE", impl_type::multi_type::event_kind::ESCAPE)
        .value("REACTION", impl_type::multi_type::event_kind::REACTION)
        ;

    // defining the python class
    class_<impl_type, bases<typename impl_type::base_type>, boost::noncopyable>(
            name,
            init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
                 typename impl_type::rng_type&, int>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
             typename impl_type::rng_type&, int, double>())
        .def(init<boost::shared_ptr<typename impl_type::world_type>,
                 boost::shared_ptr<typename impl_type::network_rules_type const>,
             typename impl_type::rng_type&, int, double,
                 typename impl_type::length_type>())
        .def("get_shell", static_cast<get_shell_result_type(impl_type::*)(typename impl_type::shell_id_type const&)>(&impl_type::get_shell))
        .def("num_domains_per_type", &impl_type::num_domains_per_type)
        .def("num_single_steps_per_type", &impl_type::num_single_steps_per_type)
        .def("num_pair_steps_per_type", &impl_type::num_pair_steps_per_type)
        .def("num_multi_steps_per_type", &impl_type::num_multi_steps_per_type)
        .def("check", &impl_type::check)
        .def("__len__", &impl_type::num_domains)
        .def("__getitem__", &impl_type::get_domain)
        .def("__iter__", &impl_type::get_domains,
                return_value_policy<return_by_value>())
        ;

    // register wrappers and converters
    peer::wrappers::generator_wrapper<
        ptr_generator<typename impl_type::domain_id_pair_generator,
                      std::auto_ptr<
                        typename impl_type::domain_id_pair_generator> > >::__register_class("DomainIDPairGenerator");

    peer::util::register_shared_const_ptr_from_python<typename impl_type::network_rules_type>();
    peer::converters::register_tuple_converter<typename impl_type::domain_id_pair>();
    peer::converters::register_tuple_converter<get_shell_result_type>();
    shell_variant_converter<Timpl>::__register();
}

} // namespace binding
#endif /* EGFRD_SIMULATOR_HPP */
