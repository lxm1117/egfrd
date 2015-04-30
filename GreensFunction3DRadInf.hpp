#if !defined(GREENSFUNCTION3DRADINF_HPP)
#define GREENSFUNCTION3DRADINF_HPP 

#include "Defs.hpp"
#include "Logger.hpp"
#include "PairGreensFunction.hpp"
#include <gsl/gsl_integration.h>

class GF_CLASS GreensFunction3DRadInf : public PairGreensFunction
{
private:
    struct p_corr_R_params;
    struct p_theta_params;

    // Error tolerance used by default.
    static const Real TOLERANCE;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.
    static const Real THETA_TOLERANCE;

    static const Real MIN_T;

    static const uint MAX_ORDER;

    static const Real H;

public:

    GreensFunction3DRadInf(Real D, Real kf, Real r0, Real Sigma);

    Real drawTime(const Real rnd) const;

    Real drawR(const Real rnd, const Real t) const;

    Real drawTheta(const Real rnd, const Real r, const Real t) const;

    Real getkD() const
    {
        return this->kD;
    }

    Real getalpha() const
    {
        return this->alpha;
    }

    Real p_reaction(Real t) const;
    Real p_survival(Real t) const;
    Real p_int_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real time) const;

    Real ip_theta(Real theta, Real r, Real time) const;

    Real p_free(Real theta, Real r, Real t) const;

    Real ip_free(Real theta, Real r, Real t) const;

    Real p_corr(Real theta, Real r, Real t) const;

    Real ip_corr(Real theta, Real r, Real t) const;

    std::string dump() const;

    const char* getName() const { return "GreensFunction3DRadInf"; }

private:
    Real p_corr_R(Real alpha, uint n, Real r, Real t) const;

    Real p_corr_n(uint n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real ip_corr_n(uint n, RealVector const& RnTable, RealVector const& lgndTable) const;

    Real p_corr_table(Real theta, Real r, Real t, RealVector const& RnTable) const;

    Real ip_corr_table(Real theta, Real r, Real t, RealVector const& RnTable) const;

    Real p_theta_table(Real r, Real theta, Real time, RealVector const& RnTable) const;

    Real ip_theta_table(Real r, Real theta, Real time, RealVector const& RnTable) const;

    void makeRnTable(RealVector& RnTable, Real r, Real t) const;

    Real Rn(uint order, Real r, Real t, gsl_integration_workspace* workspace, Real tol) const;

private:
    static Real p_corr_R_F(Real, p_corr_R_params*);
    static Real ip_theta_F(Real theta, p_theta_params* params);

    const Real kD;
    const Real alpha;

    static Logger& log_;
};

#endif // GREENSFUNCTION3DRADINF_HPP 
