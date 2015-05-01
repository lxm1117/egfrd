#ifndef GREENSFUNCTION2DABSSYM_HPP
#define GREENSFUNCTION2DABSSYM_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

class GF_CLASS GreensFunction2DAbsSym : public GreensFunction
{
    static const Real CUTOFF;           // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7, 
    static const Real CUTOFF_H;         // 5.6: ~1e-8, 6.0: ~1e-9

public:

    GreensFunction2DAbsSym(const Real D, const Real a) : GreensFunction(D), a(a) {}

    virtual ~GreensFunction2DAbsSym(){};

    virtual std::string dump() const override;

    virtual const char* getName() const override     { return "GreensFunction2DAbsSym"; }

    Real geta() const     { return a; }

    Real drawTime(const Real rnd) const;

    Real drawR(const Real rnd, const Real t) const;

    Real p_survival(const Real t) const;
    Real p_int_r(const Real r, const Real t) const;
    Real p_int_r_free(const Real r, const Real t) const;

private:

    struct p_survival_params
    {
        const GreensFunction2DAbsSym* const gf;
        const Real rnd;
    };

    struct p_r_params
    {
        const GreensFunction2DAbsSym* const gf;
        const Real t;
        const Real target;
    };

    static Real p_survival_F(const Real t, const p_survival_params* params);

    static Real p_r_free_F(const Real r, const p_r_params* params);

    static Real p_r_F(const Real r, const p_r_params* params);
    
private:
    const Real a;
    static Logger& log_;
};

#endif // GREENSFUNCTION2DABSSYM_HPP
