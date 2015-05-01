#include <sstream>
#include <exception>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include "findRoot.hpp"
#include "GreensFunction2DAbsSym.hpp"

const Real GreensFunction2DAbsSym::CUTOFF = 1e-10;
const Real GreensFunction2DAbsSym::CUTOFF_H = 6.0;

Logger& GreensFunction2DAbsSym::log_(Logger::get_logger("GreensFunction2DAbsSym"));

// an alternative form, which is not very convergent.
Real GreensFunction2DAbsSym::p_survival(const Real t) const
{
    const Real Dt(-D * t);
    const int N(100);	// number of terms to use
    Real sum(0.);
    const Real threshold(CUTOFF);
    for (int n(1); n <= N; ++n)
    {
        Real aAn = gsl_sf_bessel_zero_J0(n);		// gsl roots of J0(aAn)
        Real An = aAn / a;
        Real J1_aAn = gsl_sf_bessel_J1(aAn);
        Real term = (exp(An*An*Dt)) / (An*J1_aAn);
        sum += term;

        if (fabs(term / sum) < threshold) break;
    }
    return (2.0 / a) * sum;
}

Real GreensFunction2DAbsSym::p_int_r_free(const Real r, const Real t) const
{
    const Real Dt(D * t);
    const Real sqrtDt(sqrt(Dt));
    const Real sqrtPI(sqrt(M_PI));
    return erf(r / (sqrtDt + sqrtDt)) - r * exp(-r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

Real GreensFunction2DAbsSym::p_int_r(const Real r, const Real t) const
{
    const Real Dt(-D * t);
    Real term;
    Real sum(0.0);
    int n(1);

    //    const Real maxn( ( a / M_PI ) * sqrt( log( exp( DtPIsq_asq ) / CUTOFF ) / 
    //                                          ( D * t ) ) );

    const int N_MAX(10000);
    const Real threshold(CUTOFF);

    do
    {
        Real aAn = gsl_sf_bessel_zero_J0(n);         // gsl roots of J0(aAn)
        Real An = aAn / a;
        Real rAn = r*An;
        Real J1_aAn = gsl_sf_bessel_J1(aAn);
        Real J1_rAn = gsl_sf_bessel_J1(rAn);
        term = (exp(An*An*Dt) * r * J1_rAn) / (An*J1_aAn*J1_aAn);
        sum += term;
        n++;
    } while (fabs(term / sum) > threshold && n <= N_MAX);

    return (2.0 / (a*a)) * sum;
}

Real GreensFunction2DAbsSym::p_survival_F(const Real t, const p_survival_params* params)
{
    const GreensFunction2DAbsSym* const gf(params->gf);
    return 1 - gf->p_survival(t) - params->rnd;
}

Real GreensFunction2DAbsSym::drawTime(const Real rnd) const
{
    THROW_UNLESS(std::invalid_argument, rnd < 1.0 && rnd >= 0.0);
    if (D == 0.0 || a == INFINITY) return INFINITY;
    if (a == 0.0) return 0.0;

    p_survival_params params = { this, rnd };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_survival_F), &params };


    // Find a good interval to determine the first passage time in
    const Real t_guess(a * a / (4. * D));   // construct a guess: msd = sqrt (2*d*D*t)
    Real value(GSL_FN_EVAL(&F, t_guess));
    Real low(t_guess);
    Real high(t_guess);

    // scale the interval around the guess such that the function straddles
    if (value < 0.0)               // if the guess was too low
    {
        do
        {
            high *= 10;     // keep increasing the upper boundary until the function straddles
            value = GSL_FN_EVAL(&F, high);

            if (fabs(high) >= t_guess * 1e6)
            {
                log_.warn("Couldn't adjust high. F(%.16g) = %.16g", high, value);
                throw std::exception();
            }
        } while (value <= 0.0);
    }
    else                            // if the guess was too high
    {
        Real value_prev(value);
        do
        {
            low *= .1;      // keep decreasing the lower boundary until the function straddles
            value = GSL_FN_EVAL(&F, low);     // get the accompanying value

            if (fabs(low) <= t_guess * 1e-6 || fabs(value - value_prev) < CUTOFF)
            {
                log_.warn("Couldn't adjust high. F(%.16g) = %.16g", low, value);
                return low;
            }
            value_prev = value;
        } while (value >= 0.0);
    }

    // find the root
    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);  // a new solver type brent
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));   // make a new solver instance
    const Real t(findRoot(F, solver, low, high, 1e-18, 1e-12, "GreensFunction2DAbsSym::drawTime"));
    gsl_root_fsolver_free(solver);
    return t;
}

Real GreensFunction2DAbsSym::p_r_free_F(const Real r, const p_r_params* params)
{
    const GreensFunction2DAbsSym* const gf(params->gf);
    return gf->p_int_r_free(r, params->t) - params->target;
}

Real GreensFunction2DAbsSym::p_r_F(const Real r, const p_r_params* params)
{
    const GreensFunction2DAbsSym* const gf(params->gf);
    return gf->p_int_r(r, params->t) - params->target;
}

Real GreensFunction2DAbsSym::drawR(const Real rnd, const Real t) const
{
    THROW_UNLESS(std::invalid_argument, rnd <= 1.0 && rnd >= 0.0);
    THROW_UNLESS(std::invalid_argument, t >= 0.0);

    if (a == 0.0 || t == 0.0 || D == 0.0) return 0.0;

    //const Real thresholdDistance( CUTOFF_H * sqrt( 4.0 * D * t ) );

    //  if( a <= thresholdDistance )	// if the domain is not so big, the boundaries are felt
    //  {
    Real psurv = p_survival(t);
    //psurv = p_int_r( a, t );
    //printf("dr %g %g\n",psurv, p_survival( t ));
    //assert( fabs(psurv - p_int_r( a, t )) < psurv * 1e-8 );

    assert(psurv > 0.0);

    gsl_function F;
    F.function = reinterpret_cast<double(*)(double, void*)>(&p_r_F);
    /*  }
        else				// if the domain is very big, just use the free solution
        {
        // p_int_r < p_int_r_free
        if( p_int_r_free( a, t ) < rnd )	// if the particle is outside the domain?
        {
        std::cerr << "p_int_r_free( a, t ) < rnd, returning a."
        << std::endl;
        return a;
        }

        psurv = 1.0;
        F.function = reinterpret_cast<double(*)(double, void*)>( &p_r_free_F );
        }

        */
    const Real target(psurv * rnd);
    p_r_params params = { this, t, target };
    F.params = &params;

    const Real low(0.0);
    const Real high(a);
    //const Real high( std::min( thresholdDistance, a ) );

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    const Real r(findRoot(F, solver, low, high, 1e-18, 1e-12, "GreensFunction2DAbsSym::drawR"));
    gsl_root_fsolver_free(solver);
    return r;
}

std::string GreensFunction2DAbsSym::dump() const
{
    std::ostringstream ss;
    ss << "D = " << D << ", a = " << a << std::endl;
    return ss.str();
}
