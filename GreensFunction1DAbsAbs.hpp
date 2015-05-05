#ifndef GREENSFUNCTION1DABSABS_HPP
#define GREENSFUNCTION1DABSABS_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

class GF_CLASS GreensFunction1DAbsAbs : public GreensFunction
{
    static const Real L_TYPICAL;                // This is a typical length scale of the system, may not be true!
    static const Real T_TYPICAL;                // The typical timescale of the system, may also not be true!!
    static const Real EPSILON;                  // measure of 'sameness' when comparing floating points numbers
    static const Real PDENS_TYPICAL;            //E3; Is 1E3 a good measure for the probability density?!
    static const uint MAX_TERMS;                // The maximum number of terms in the sum
    static const uint MIN_TERMS;                // The minimum
    static const Real CUTOFF_H;                 // Cutoff distance: When H * sqrt(2Dt) < 1/2*L, use free greens function instead of absorbing.

public:
    GreensFunction1DAbsAbs(Real D, Real r0, Real sigma, Real a)
        : GreensFunction(D), v(0.0), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
    }

    // The constructor is overloaded and can be called with or without drift v
    GreensFunction1DAbsAbs(Real D, Real v, Real r0, Real sigma, Real a) // copy constructor including drift variable v
        : GreensFunction(D), v(v), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
    }

    virtual ~GreensFunction1DAbsAbs(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction1DAbsAbs"; }

    Real getsigma() const { return sigma; }

    Real geta() const     { return a; }

    Real getv() const     { return v; }

    Real getr0() const    { return r0; }

    // Draws the first passage time from the propensity function
    Real drawTime(Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the domain
    Real drawR(Real rnd, Real t) const;

    // Calculates the amount of flux leaving the left boundary at time t
    Real leaves(Real t) const;

    // Calculates the amount of flux leaving the right boundary at time t
    Real leavea(Real t) const;

    // Determines based on the flux ratios if the particle left the left or right boundary
    EventKind drawEventType(Real rnd, Real t) const;

    // Calculates the probability of finding the particle inside the 
    // domain at time t so, the survival probability
    Real p_survival(Real t) const;

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the domain.
    Real calcpcum(Real r, Real t) const;

    // Calculates the probability density of finding the particle at location r at time t.
    Real prob_r(Real r, Real t) const;

private:

    uint guess_maxi(Real const& t) const;

    Real p_survival_table(Real  t, RealVector& psurvTable) const;
    Real p_survival_i(uint i, Real const& t, RealVector const& table) const;

    Real p_survival_table_i_v(uint const& i) const;

    Real p_survival_table_i_nov(uint const& i) const;

    void createPsurvTable(uint const& maxi, RealVector& table) const;

    Real p_int_r_table(Real const& r, Real const& t, RealVector& table) const;

    Real p_int_r_i(uint i, Real const& r, Real const& t, RealVector& table) const;

    void create_p_int_r_Table(Real const& t, uint const& maxi, RealVector& table) const;

    Real get_p_int_r_Table_i(uint& i, Real const& t, RealVector& table) const
    {
        if (i >= table.size())
            create_p_int_r_Table(t, i + 1, table);
        return table[i];
    }

private:
    const Real v;         // The diffusion constant and drift velocity
    const Real sigma;     // These are the dimensions of our domain; L is calculated as a-sigma
    const Real a;
    const Real r0;
    const Real l_scale;       // This is the 'length scale' of your system (1e-14 or 1e6), Although rescaling is discontinued, we use it to check whether a is well-chosen
    const Real t_scale;       // This is the time scale of the system, used by drawTime_f

    static Logger& log_;
};

#endif // GREENSFUNCTION1DABSABS_HPP
