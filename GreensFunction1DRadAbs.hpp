#ifndef GREENSFUNCTION1DRADABS_HPP
#define GREENSFUNCTION1DRADABS_HPP

#include "Defs.hpp"
#include "Logger.hpp"
#include "GreensFunction.hpp"

class GF_CLASS GreensFunction1DRadAbs : public GreensFunction
{
    static const Real L_TYPICAL;        // This is a typical length scale of the system, may not be true!
    static const Real T_TYPICAL;        // The typical timescale of the system, may also not be true!!
    static const Real EPSILON;          // measure of 'sameness' when comparing floating points numbers
    static const Real PDENS_TYPICAL;    // Is 1E3 a good measure for the probability density?!
    static const uint MAX_TERMS;        // The maximum number of terms used in calculating the sum
    static const uint MIN_TERMS;        // The minimum number of terms
    static const Real CUTOFF_H;         // Cutoff distance: When H * sqrt(2Dt) < a - r0 OR ro - sigma  use free greens function instead of absorbing. 

public:
    GreensFunction1DRadAbs(Real D, Real k, Real r0, Real sigma, Real a)
        : GreensFunction(D), v(0.0), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
        //set first root.
        calculate_n_roots(1);
    }

    // The constructor is overloaded and can be called with or without drift v
    // copy constructor including drift variable v
    GreensFunction1DRadAbs(Real D, Real v, Real k, Real r0, Real sigma, Real a)
        : GreensFunction(D), v(v), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
        //set first root;
        calculate_n_roots(1);
    }

    virtual ~GreensFunction1DRadAbs(){}

    virtual std::string dump() const override;

    virtual const char* getName() const override { return "GreensFunction1DRadAbs"; }

    Real geta() const     { return a; }

    Real getsigma() const { return sigma; }

    Real getr0() const    { return r0; }

    Real getk() const    { return k; }

    Real getv() const    { return v; }

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    Real calcpcum(Real r, Real t) const;

    // Determine which event has occurred, an escape or a reaction. Based 
    // on the fluxes through the boundaries at the given time. Beware: if 
    // t is not a first passage time you still get an answer!
    EventKind drawEventType(Real rnd, Real t) const;

    // Draws the first passage time from the propensity function
    Real drawTime(Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the domain
    Real drawR(Real rnd, Real t) const;

    // These methods are both public and private, they are used by public methods 
    // but can also be called from the 'outside'. This is mainly because of 
    // debugging purposes.

    // Calculates the probability of finding the particle inside the 
    // domain at time t -> the survival probability
    Real p_survival(Real t) const;

    // Calculates the total probability flux leaving the domain at time t
    Real flux_tot(Real t) const;

    // Calculates the probability flux leaving the domain through the 
    // radiative boundary at time t
    Real flux_rad(Real t) const;

    // Calculates the flux leaving the domain through the radiative 
    // boundary as a fraction of the total flux. This is the probability 
    // that the particle left the domain through the radiative
    // boundary instead of the absorbing boundary.
    Real fluxRatioRadTot(Real t) const;

    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r(Real r, Real t) const;

private:

    Real An(Real a_n) const;
    Real Bn(Real a_n) const;
    Real Cn(Real a_n, Real t) const;

    struct tan_f_params
    {
        Real a;
        Real h;
    };

    struct drawT_params
    {
        GreensFunction1DRadAbs const* gf;
        RealVector& psurvTable;
        Real rnd;
    };

    struct drawR_params
    {
        GreensFunction1DRadAbs const* gf;
        const Real t;
        RealVector table;
        Real rnd;
    };

    /* Functions managing the rootList */

    /* return the rootList size */
    uint rootList_size() const { return rootList.size(); }

    /* return the n + 1'th root */
    Real get_root(uint n) const;

    /* Check the rootList for the first n roots. */
    void calculate_n_roots(uint n) const;

    /* Fills the rootList from i to n. */
    void fill_table_to_n(uint i, uint n);

    /* Guess the number of terms needed for convergence, given t. */
    uint guess_maxi(Real const& t) const;

    /* this is the appropriate definition of the function in gsl. */
    static Real tan_f(Real x, void *p);

    /* functions for drawTime / p_survival */

    static Real drawT_f(Real t, void *p);

    Real p_survival_table(Real  t, RealVector& psurvTable) const;

    Real p_survival_i(uint i, Real const& t, RealVector const& table) const;

    Real p_survival_table_i_v(uint const& i) const;

    Real p_survival_table_i_nov(uint const& i) const;

    void createPsurvTable(RealVector& table) const;

    /* functions for drawR */

    static Real drawR_f(Real z, void* p);

    Real p_int_r_table(Real const& r, Real const& t, RealVector& table) const;

    Real p_int_r_i(uint i, Real const& r, Real const& t, RealVector& table) const;

    void create_p_int_r_Table(Real const& t, RealVector& table) const;

    Real get_p_int_r_Table_i(uint& i, Real const& t, RealVector& table) const
    {
        if (i >= table.size())
        {
            calculate_n_roots(i + 1);
            create_p_int_r_Table(t, table);
        }

        return table[i];
    }

private:
    const Real v;         // The diffusion constant and drift velocity
    const Real k;         // The reaction constant
    const Real r0;
    const Real sigma;     // The left and right boundary of the domain (sets the l_scale, see below)
    const Real a;
    const Real l_scale;   // This is the length scale of the system
    const Real t_scale;   // This is the time scale of the system.

    RealVector rootList;            /* vector containing the roots 0f tan_f. */
    static Logger& log_;
};
#endif // GREENSFUNCTION1DRADABS_HPP
