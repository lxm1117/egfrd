#ifndef SHELL_ID_HPP
#define SHELL_ID_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif
#include "Identifier.hpp"


struct ShellID : public Identifier < ShellID, unsigned int >
{
    typedef Identifier<ShellID, unsigned int> base_type;
    ShellID(value_type const& value = value_type(0)) : base_type(value) {}
};

template<> struct std::hash < ShellID >
{
    std::size_t operator()(ShellID const& id) const
    {
        return static_cast<std::size_t>(id());
    }
};

inline std::ostream& operator<<(std::ostream& stream, const ShellID& id)
{
    stream << "ShellID(" << id() << ")";
    return stream;
};

#endif /* SHELL_ID_HPP */
