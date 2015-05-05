#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include "Defs.hpp"

class GF_CLASS GreensFunction
{
public:
    enum EventKind { IV_ESCAPE, IV_REACTION };

    GreensFunction(const Real D) : D(D) {}
    virtual ~GreensFunction() { }

    // prevent default construction, assignment and move operations, allow copy construction
    GreensFunction() = delete;
    GreensFunction(const GreensFunction&) = default;
    GreensFunction& operator=(const GreensFunction&) = delete;
    GreensFunction(const GreensFunction&&) = delete;
    GreensFunction& operator=(const GreensFunction&&) = delete;

    Real getD() const { return D; }

    virtual std::string dump() const = 0;;
    virtual const char* getName() const = 0;

protected:
    const Real D;
};

#endif // GREENSFUNCTION_HPP
