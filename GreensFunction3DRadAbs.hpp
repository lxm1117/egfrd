#ifndef GREENSFUNCTION3DRADABS_HPP
#define GREENSFUNCTION3DRADABS_HPP

#include <vector>
#include <boost/array.hpp>
#include <gsl/gsl_roots.h>
#include "Logger.hpp"
#include "PairGreensFunction.hpp"

class GF_CLASS GreensFunction3DRadAbs : public PairGreensFunction
{
    static const Real TOLERANCE;            // Error tolerance used by default.
    static const Real THETA_TOLERANCE;      // SphericalBesselGenerator's accuracy, used by some theta-related calculations.
    static const Real MIN_T_FACTOR;
    static const uint MAX_ORDER;
    static const uint MAX_ALPHA_SEQ;

public:

    GreensFunction3DRadAbs(Real D, Real kf, Real r0, Real Sigma, Real a);

    virtual ~GreensFunction3DRadAbs(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction3DRadAbs"; }

    Real geth() const { return h; }

    Real geta() const { return a; }

    virtual Real drawR(const Real rnd, const Real t) const override;

    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const override;

    Real drawTime(const Real rnd) const;

    EventKind drawEventType(const Real rnd, const Real t) const;

    Real f_alpha0(Real alpha) const;
    Real f_alpha0_aux(Real alpha) const;

    Real f_alpha(Real alpha, Integer n) const;
    Real f_alpha_aux(Real alpha, Integer n) const;

    Real p_0(Real t, Real r) const;

    Real p_survival(Real t) const;

    Real p_survival_table(Real t, RealVector& table) const;

    Real p_leave_table(Real t, RealVector const& table) const;

    Real dp_survival(Real t) const;

    Real leaves(Real t) const;

    Real leavea(Real t) const;

    Real p_leaves(Real t) const;

    Real p_leavea(Real t) const;

    Real p_int_r(Real r, Real t) const;

    Real p_theta(Real theta, Real r, Real t) const;

    Real ip_theta(Real theta, Real r, Real t) const;

    Real dp_theta(Real theta, Real r, Real t) const;

    Real idp_theta(Real theta, Real r, Real t) const;

    Real p_n(Integer n, Real r, Real t, Real max_alpha) const;

    Real dp_n_at_a(Integer n, Real t, Real max_alpha) const;

    Real p_n_alpha(uint i, uint n, Real r, Real t) const;

    Real dp_n_alpha_at_a(uint i, uint n, Real t) const;

    // methods below are kept public for debugging purpose.

    uint alphaOffset(uint n) const;

    Real alpha0_i(Integer i) const;

    Real alpha_i(Integer i, Integer n, gsl_root_fsolver* solver) const;

    Real p_survival_i(Real alpha) const;

    Real p_0_i(Real alpha, Real r) const;

    Real dp_survival_i(Real alpha) const;

    Real leavea_i(Real alpha) const;

    Real leaves_i(Real alpha) const;

    Real p_leavea_i(Real alpha, Real pleave_factor) const;

    Real p_leaves_i(Real alpha, Real pleave_factor) const;

    Real p_survival_den(Real alpha) const;

    Real p_int_r_i(Real r, Real alpha, Real num_r0) const;

    Real p_0_i_exp(uint i, Real t, Real r) const;

    Real p_survival_i_exp(uint i, Real t) const;

    Real p_survival_i_alpha(Real alpha, Real t) const;

    Real p_survival_2i_exp(uint i, Real t) const;

protected:

    void clearAlphaTable() const;

    RealVector& getAlphaTable(size_t n) const
    {
        assert(alphaTable.size() > n);
        return alphaTable[n];
    }

    Real getAlpha(size_t n, RealVector::size_type i) const
    {
        assert(alphaTable.size() > n);

        RealVector& aTable(alphaTable[n]);
        RealVector::size_type oldSize(aTable.size());

        if (oldSize <= i)
        {
            aTable.resize(i + 1);
            uint offset(alphaOffset(n));
            gsl_root_fsolver* solver(gsl_root_fsolver_alloc(gsl_root_fsolver_brent));
            for (RealVector::size_type m(oldSize); m <= i; ++m)
                aTable[m] = alpha_i(m + offset, n, solver);
            gsl_root_fsolver_free(solver);
        }

        return aTable[i];

    }

    Real getAlpha0(RealVector::size_type i) const
    {
        RealVector& aTable(alphaTable[0]);
        RealVector::size_type oldSize(aTable.size());

        if (oldSize <= i)
        {
            aTable.resize(i + 1);
            for (RealVector::size_type m(oldSize); m <= i; ++m)
                aTable[m] = alpha0_i(m);
        }
        return aTable[i];
    }

    Real p_int_r_table(Real r, Real t, RealVector const& num_r0Table) const;

    Real ip_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const;

    Real p_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const;

    void make_p_thetaTable(RealVector& pTable, Real r, Real t, uint n, RealVector const& p_nTable) const;

    Real p_survival_i_exp_table(uint i, Real t, RealVector const& table) const;

    Real p_leave_i_exp_table(uint i, Real t, RealVector const& table) const;

    Real dp_survival_i_exp(uint i, Real alpha) const;

    Real leavea_i_exp(uint i, Real alpha) const;

    Real leaves_i_exp(uint i, Real alpha) const;

    Real p_leavea_i_exp(uint i, Real alpha) const;

    Real p_leaves_i_exp(uint i, Real alpha) const;

    Real p_int_r_i_exp(uint i, Real t, Real r) const;

    Real p_int_r_i_exp_table(uint i, Real t, Real r, RealVector& num_r0Table) const;

    void updateAlphaTable0(Real t) const;
    void updateAlphaTable(uint n, Real t) const;

    void createPsurvTable(RealVector& table) const;
    void createNum_r0Table(RealVector& table) const;

    void createPleaveFactorTable(RealVector& table) const;
    void createPleavesTable(RealVector& table, RealVector const& pleaveFactorTable) const;
    void createPleaveaTable(RealVector& table, RealVector const& pleaveFactorTable) const;

    void makep_nTable(RealVector& p_nTable, Real r, Real t) const;

    void makedp_n_at_aTable(RealVector& p_nTable, Real t) const;

    uint guess_maxi(Real t) const;

    Real drawPleaves(gsl_function const& F, gsl_root_fsolver* solver, Real t_guess, RealVector& pleaveFactorTable, RealVector& pleavesTable) const;

    Real drawPleavea(gsl_function const& F, gsl_root_fsolver* solver, Real t_guess, RealVector& pleaveFactorTable, RealVector& pleavesTable) const;

    Real num_r0(Real alpha) const;

    Real pleaveFactor(Real alpha) const;

    struct ip_theta_params;
    static Real ip_theta_F(Real, ip_theta_params const*);

private:
    const Real a;
    const Real h;
    const Real hsigma_p_1;

    mutable boost::array<Integer, GF_MAX_ORDER> alphaOffsetTable;
    mutable boost::array<RealVector, GF_MAX_ORDER> alphaTable;

    static Logger& log_;
};

#endif // GREENSFUNCTION3DRADABS_HPP
