#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include "Defs.hpp"

class GF_CLASS GreensFunction
{
public:
    enum EventKind
    {
        IV_ESCAPE,
        IV_REACTION
    };

public:
    GreensFunction(const Real D) : D(D) {}

    Real getD() const
    {
        return this->D;
    }

protected:
    const Real D;
};

#endif // GREENSFUNCTION_HPP
