#ifndef FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP
#define FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP

#include "PairGreensFunction.hpp"

class GreensFunction3DRadAbsBase : public PairGreensFunction
{
public:
    GreensFunction3DRadAbsBase(Real D, Real kf, Real r0, Real Sigma) : PairGreensFunction(D, kf, r0, Sigma) {}

    virtual Real drawTime(const Real rnd) const = 0;

    virtual EventKind drawEventType(const Real rnd, const Real t) const = 0;

    virtual Real drawR(const Real rnd, const Real t) const = 0;

    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const = 0;
};

#endif /* FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP */
