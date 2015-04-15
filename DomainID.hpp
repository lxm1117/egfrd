#ifndef DOMAIN_ID_HPP
#define DOMAIN_ID_HPP

#include <ostream>
#include <functional>
#include "Identifier.hpp"

struct DomainID : public Identifier < DomainID, unsigned int >
{
    typedef Identifier<DomainID, unsigned int> base_type;
    DomainID(value_type const& value = value_type(0)) : base_type(value) {}
};

namespace std
{
    template<>
    struct hash < DomainID >
    {
        inline size_t  operator()(const DomainID & id) const
        {
            return hash<DomainID::value_type>()(id());
        }
    };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const DomainID& id)
{
    stream << "DID(" << id() << ")";
    return stream;
};

#endif /* DOMAIN_ID_HPP */
