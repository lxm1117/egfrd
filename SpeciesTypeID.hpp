#ifndef SPECIES_TYPE_ID_HPP
#define SPECIES_TYPE_ID_HPP

#include <ostream>
#include <functional>
#include "Identifier.hpp"

struct SpeciesTypeID : public Identifier < SpeciesTypeID, unsigned int >
{
    typedef Identifier<SpeciesTypeID, unsigned int> base_type;
    SpeciesTypeID(value_type const& value = value_type(0)) : base_type(value) {}
};

namespace std
{
    template<>
    struct hash < SpeciesTypeID >
    {
        inline size_t  operator()(const SpeciesTypeID & id) const
        {
            return hash<SpeciesTypeID::value_type>()(id());
        }
    };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const SpeciesTypeID& id)
{
    stream << "SpeciesTypeID(" << id() << ")";
    return stream;
};

#endif /* SPECIES_TYPE_ID_HPP */
