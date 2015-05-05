#ifndef GREENSFUNCTION3DABS_HPP
#define GREENSFUNCTION3DABS_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "PairGreensFunction.hpp"

class GF_CLASS GreensFunction3DAbs : public PairGreensFunction
{
    static const Real TOLERANCE;            // Error tolerance used by default.
    static const Real THETA_TOLERANCE;      // SphericalBesselGenerator's accuracy, used by some theta-related calculations.
    static const Real MIN_T;
    static const uint MAX_ORDER;
    static const uint MAX_ALPHA_SEQ;

public:

    GreensFunction3DAbs(Real D, Real r0, Real a);

    virtual ~GreensFunction3DAbs(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction3DAbs"; }

    virtual Real drawR(const Real rnd, const Real t) const override;

    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const override;

    Real geta() const { return a; }

    Real drawTime(const Real rnd) const;

    EventKind drawEventType(const Real rnd, const Real t) const;

    Real p_survival(Real t) const;

    Real dp_survival(Real t) const;

    Real p_int_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t) const;

    Real dp_theta(Real theta, Real r, Real t) const;

    Real idp_theta(Real theta, Real r, Real t) const;

    Real p_n(int n, Real r, Real t) const;

    Real dp_n(int n, Real t) const;

    Real p_n_alpha(uint i, uint n, Real r, Real t) const;

    Real dp_n_alpha(uint i, uint n, Real t) const;

private:

    Real p_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const;

    Real ip_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const;

    void makep_nTable(RealVector& p_nTable, Real r, Real t) const;

    void makedp_nTable(RealVector& p_nTable, Real t) const;

    struct ip_theta_params;
    static Real ip_theta_F(Real theta, ip_theta_params const* params);

    const Real a;
    static Logger& log_;
};

#endif // GREENSFUNCTION3DABS_HPP
