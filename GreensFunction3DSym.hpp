#ifndef GREENSFUNCTION3DSYM_HPP
#define GREENSFUNCTION3DSYM_HPP

#include "compat.h"
#include <gsl/gsl_integration.h>
#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

/**
  Green's Function for a free diffusion particle.
  */

class GF_CLASS GreensFunction3DSym : public GreensFunction
{
private:

    static const Real TOLERANCE;

public:

    GreensFunction3DSym(const Real D) : GreensFunction(D) { }

    Real drawTime(const Real) const { return INFINITY; }

    Real drawR(const Real rnd, const Real t) const;

    Real p_r(const Real r, const Real t) const;

    Real ip_r(const Real r, const Real t) const;

    std::string dump() const;

    const char* getName() const { return "GreensFunction3DSym"; }

private:

    static Logger& log_;
};

#endif // GREENSFUNCTION3DSYM_HPP
