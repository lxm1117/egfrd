#ifndef GREENSFUNCTION3DSYM_HPP
#define GREENSFUNCTION3DSYM_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

/**
  Green's Function for a free diffusion particle.
  */

class GF_CLASS GreensFunction3DSym : public GreensFunction
{
    static const Real TOLERANCE;

public:

    GreensFunction3DSym(const Real D) : GreensFunction(D) { }

    virtual ~GreensFunction3DSym(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override{ return "GreensFunction3DSym"; }

    Real drawTime(const Real) const { return INFINITY; }

    Real drawR(const Real rnd, const Real t) const;

    Real p_r(const Real r, const Real t) const;

    Real ip_r(const Real r, const Real t) const;

private:

    static Logger& log_;
};

#endif // GREENSFUNCTION3DSYM_HPP
