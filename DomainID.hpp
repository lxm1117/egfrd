#ifndef DOMAIN_ID_HPP
#define DOMAIN_ID_HPP

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


struct DomainID : public Identifier < DomainID, unsigned int >
{
    typedef Identifier<DomainID, unsigned int> base_type;
    DomainID(value_type const& value = value_type(0)) : base_type(value) {}
};

template<> struct std::hash < DomainID >
{
    std::size_t operator()(DomainID const& id) const
    {
        return static_cast<std::size_t>(id());
    }
};

inline std::ostream& operator<<(std::ostream& stream, const DomainID& id)
{
    stream << "DID(" << id() << ")";
    return stream;
};

#endif /* DOMAIN_ID_HPP */
