#ifndef SHELL_ID_HPP
#define SHELL_ID_HPP

#include <ostream>
#include <functional>
#include "Identifier.hpp"

struct ShellID : public Identifier < ShellID, unsigned int >
{
    typedef Identifier<ShellID, unsigned int> base_type;
    ShellID(value_type const& value = value_type(0)) : base_type(value) {}
};

namespace std
{
    template<>
    struct hash < ShellID >
    {
        inline size_t  operator()(const ShellID & id) const
        {
            return hash<ShellID::value_type>()(id());
        }
    };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const ShellID& id)
{
    stream << "ShellID(" << id() << ")";
    return stream;
};

#endif /* SHELL_ID_HPP */
