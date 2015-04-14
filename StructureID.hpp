#ifndef STRUCTURE_ID_HPP
#define STRUCTURE_ID_HPP

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


struct StructureID : public Identifier < StructureID, unsigned int >
{
    typedef Identifier<StructureID, unsigned int> base_type;
    StructureID(value_type const& value = value_type(0)) : base_type(value) {}
};

template<> struct std::hash < StructureID >
{
    std::size_t operator()(StructureID const& id) const
    {
        return static_cast<std::size_t>(id());
    }
};

inline std::ostream& operator<<(std::ostream& stream, const StructureID& id)
{
    stream << "StructID(" << id() << ")";
    return stream;
};

#endif /* STRUCTURE_ID_HPP */
