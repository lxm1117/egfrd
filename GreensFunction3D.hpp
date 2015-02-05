#ifndef GREENSFUNCTION3D_HPP
#define GREENSFUNCTION3D_HPP

#include "compat.h"
#include <gsl/gsl_integration.h>
#include "Logger.hpp"
#include "PairGreensFunction.hpp"

/**
   Pair Green's function for the case where the pair never interact.

   Therefore, drawTime() always returns +INFINITY.
   kf == sigma == 0.
   */

class GreensFunction3D : public PairGreensFunction
{
private:
    static const Real TOLERANCE;
    static const Real H;

public:

    GreensFunction3D(Real D, Real r0) : PairGreensFunction(D, 0.0, r0, 0.0)
    {
    }

    virtual Real drawTime(const Real rnd) const;

    Real drawR(const Real rnd, const Real t) const;

    Real drawTheta(const Real rnd, const Real r, const Real t) const;

    Real p_r(Real r, Real t) const;

    Real ip_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t) const;

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3D";
    }

private:
    static Logger& log_;
};

#endif // GREENSFUNCTION3D_HPP
