#ifndef STRUCTURE_ID_HPP
#define STRUCTURE_ID_HPP

#include <ostream>
#include <functional>
#include "Identifier.hpp"

struct StructureID : public Identifier < StructureID, unsigned int >
{
    typedef Identifier<StructureID, unsigned int> base_type;
    StructureID(value_type const& value = value_type(0)) : base_type(value) {}
};

namespace std
{
    template<>
    struct hash < StructureID >
    {
        inline size_t  operator()(const StructureID & id) const
        {
            return hash<StructureID::value_type>()(id());
        }
    };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const StructureID& id)
{
    stream << "StructID(" << id() << ")";
    return stream;
};

#endif /* STRUCTURE_ID_HPP */
