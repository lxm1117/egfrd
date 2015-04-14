#ifndef SPECIES_TYPE_ID_HPP
#define SPECIES_TYPE_ID_HPP

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


struct SpeciesTypeID : public Identifier < SpeciesTypeID, unsigned int >
{
    typedef Identifier<SpeciesTypeID, unsigned int> base_type;
    SpeciesTypeID(value_type const& value = value_type(0)) : base_type(value) {}
};

template<> struct std::hash < SpeciesTypeID >
{
    std::size_t operator()(SpeciesTypeID const& id) const
    {
        return static_cast<std::size_t>(id());
    }
};

inline std::ostream& operator<<(std::ostream& stream, const SpeciesTypeID& id)
{
    stream << "SpeciesID(" << id() << ")";
    return stream;
};

#endif /* SPECIES_TYPE_ID_HPP */
