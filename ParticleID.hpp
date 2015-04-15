#ifndef PARTICLE_ID_HPP
#define PARTICLE_ID_HPP

#include <ostream>
#include <functional>
#include "Identifier.hpp"

struct ParticleID : public Identifier < ParticleID, unsigned int >
{
    typedef Identifier<ParticleID, unsigned int> base_type;
    ParticleID(value_type const& value = value_type(0)) : base_type(value) {}
};

namespace std
{
    template<>
    struct hash < ParticleID >
    {
        inline size_t  operator()(const ParticleID & id) const
        {
            return hash<ParticleID::value_type>()(id());
        }
    };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const ParticleID& id)
{
    stream << "PID(" << id() << ")";
    return stream;
};

#endif /* PARTICLE_ID_HPP */
