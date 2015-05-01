#ifndef GREENSFUNCTION3DABSSYM_HPP
#define GREENSFUNCTION3DABSSYM_HPP

#include <ostream>
#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

class GF_CLASS GreensFunction3DAbsSym : public GreensFunction
{
    static const Real CUTOFF;       // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    static const Real CUTOFF_H;     // 5.6: ~1e-8, 6.0: ~1e-9

public:
    GreensFunction3DAbsSym(Real D, Real a) : GreensFunction(D), a(a) {}

    virtual ~GreensFunction3DAbsSym(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction3DAbsSym"; }

    Real geta() const { return a; }

    Real drawTime(const Real rnd) const;

    Real drawR(Real rnd, Real t) const;

    Real p_survival(Real t) const;
    Real p_int_r(Real r, Real t) const;
    Real p_int_r_free(Real r, Real t) const;
    Real p_r_fourier(Real r, Real t) const;

private:

    static Real ellipticTheta4Zero(Real q);

    const Real a;
    static Logger& log_;
};

#endif // GREENSFUNCTION3DABSSYM_HPP
