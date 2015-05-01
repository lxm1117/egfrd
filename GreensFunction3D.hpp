#ifndef GREENSFUNCTION3D_HPP
#define GREENSFUNCTION3D_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "PairGreensFunction.hpp"

/**
   Pair Green's function for the case where the pair never interact.

   Therefore, drawTime() always returns +INFINITY.
   kf == sigma == 0.
   */

class GF_CLASS GreensFunction3D : public PairGreensFunction
{
    static const Real TOLERANCE;
    static const Real H;

public:

    GreensFunction3D(Real D, Real r0) : PairGreensFunction(D, 0.0, r0, 0.0) { }

    virtual ~GreensFunction3D(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction3D"; }

    virtual Real drawR(const Real rnd, const Real t) const override;

    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const override;

    Real drawTime(const Real rnd) const;

    Real p_r(Real r, Real t) const;

    Real ip_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t) const;

private:
    static Logger& log_;
};

#endif // GREENSFUNCTION3D_HPP
