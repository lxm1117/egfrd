#ifndef GREENSFUNCTION3DABSSYM_HPP
#define GREENSFUNCTION3DABSSYM_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"
#include <ostream>

class GF_CLASS GreensFunction3DAbsSym : public GreensFunction
{
    static const Real CUTOFF;

    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    // 5.6: ~1e-8, 6.0: ~1e-9
    static const Real CUTOFF_H;

public:
    GreensFunction3DAbsSym(Real D, Real a) : GreensFunction(D), a(a) {}

    Real getD() const
    {
        return this->D;
    }

    Real geta() const
    {
        return this->a;
    }

    Real p_survival(Real t) const;

    Real drawTime(const Real rnd) const;

    Real drawR(Real rnd, Real t) const;

    Real p_int_r(Real r, Real t) const;
    Real p_int_r_free(Real r, Real t) const;

    Real p_r_fourier(Real r, Real t) const;

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3DAbsSym";
    }

private:

    static Real ellipticTheta4Zero(Real q);

private:

    const Real a;
    static Logger& log_;
};

template<typename Tstrm, typename Ttraits>
inline std::basic_ostream<Tstrm, Ttraits>& operator <<(std::basic_ostream<Tstrm, Ttraits>& strm, GreensFunction3DAbsSym const& gf)
{
    strm << "GreensFunction3DAbsSym("
        << "D=" << gf.getD() << ", "
        << "a=" << gf.geta() << ")";
    return strm;
}

#endif // GREENSFUNCTION3DABSSYM_HPP
