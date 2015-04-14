#ifndef PARTICLE_ID_HPP
#define PARTICLE_ID_HPP

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


struct ParticleID : public Identifier < ParticleID, unsigned int >
{
    typedef Identifier<ParticleID, unsigned int> base_type;
    ParticleID(value_type const& value = value_type(0)) : base_type(value) {}
};

template<> struct std::hash < ParticleID >
{
    std::size_t operator()(ParticleID const& id) const
    {
        return static_cast<std::size_t>(id());
    }
};

inline std::ostream& operator<<(std::ostream& stream, const ParticleID& id)
{
    stream << "PID(" << id() << ")";
    return stream;
};

#endif /* PARTICLE_ID_HPP */
