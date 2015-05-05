#include <stdexcept>
#include <sstream>
#include <boost/format.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "GreensFunction3DSym.hpp"

const Real GreensFunction3DSym::TOLERANCE = 1e-8;

Logger& GreensFunction3DSym::log_(Logger::get_logger("GreensFunction3DSym"));

Real GreensFunction3DSym::p_r(Real r, Real t) const
{
    const Real Dt(D * t);
    const Real Dt4(4.0 * Dt);
    const Real Dt4Pi(Dt4 * M_PI);
    const Real term1(1.0 / sqrt(gsl_pow_3(Dt4Pi)));
    const Real term2(exp(-r * r / Dt4));
    const Real jacobian(4.0 * r * r * M_PI);
    return jacobian * term1 * term2;
}

Real GreensFunction3DSym::ip_r(Real r, Real t) const
{
    const Real Dt(D * t);
    const Real sqrtDt_r(1.0 / sqrt(D * t));
    const Real sqrtPi_r(1.0 / sqrt(M_PI));
    const Real term1(exp(-r * r / (4.0 * Dt)) * r * sqrtDt_r * sqrtPi_r);
    const Real term2(erf(r * 0.5 * sqrtDt_r));
    return term2 - term1;
}

struct ip_r_params
{
    GreensFunction3DSym const* const gf;
    const Real t;
    const Real value;
};

static Real ip_r_F(Real r, ip_r_params const* params)
{
    const GreensFunction3DSym* const gf(params->gf);
    return gf->ip_r(r, params->t) - params->value;
}

Real GreensFunction3DSym::drawR(Real rnd, Real t) const
{
    // input parameter range checks.
    if (!(rnd <= 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DSym: rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DSym: t >= 0.0 : t=%.16g") % t).str());
    }

    // t == 0 or D == 0 means no move.
    if (t == 0.0 || D == 0.0) return 0.0;

    ip_r_params params = { this, t, rnd };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&ip_r_F), &params };
    Real max_r(4.0 * sqrt(6.0 * D * t));

    while (GSL_FN_EVAL(&F, max_r) < 0.0)
    {
        max_r *= 10;
    }

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, 0.0, max_r);

    const uint maxIter(100);
    for (uint i(0);;++i)
    {
        gsl_root_fsolver_iterate(solver);
        const Real low(gsl_root_fsolver_x_lower(solver));
        const Real high(gsl_root_fsolver_x_upper(solver));
        const int status(gsl_root_test_interval(low, high, 1e-15, TOLERANCE));

        if (status != GSL_CONTINUE) break;
        if (i >= maxIter) throw std::runtime_error("GreensFunction3DSym: drawR: failed to converge");
    }

    const Real r(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);
    return r;
}

std::string GreensFunction3DSym::dump() const
{
    std::ostringstream ss;
    ss << "D = " << D << std::endl;
    return ss.str();
}
