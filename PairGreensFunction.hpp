#ifndef PAIRGREENSFUNCTION_HPP
#define PAIRGREENSFUNCTION_HPP

#include "Defs.hpp"
#include "GreensFunction.hpp"

class GF_CLASS PairGreensFunction : public GreensFunction
{
public:
    PairGreensFunction(Real D, Real kf, Real r0, Real sigma) : GreensFunction(D), kf(kf), r0(r0), sigma(sigma) {}

    Real getkf() const { return kf; }
    Real getSigma() const { return sigma; }
    Real getr0() const { return r0; }

    virtual Real drawR(const Real rnd, const Real t) const = 0;
    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const = 0;

protected:
    const Real kf;
    const Real r0;
    const Real sigma;
};

#endif /* PAIRGREENSFUNCTION_HPP */
