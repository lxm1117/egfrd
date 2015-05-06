#include <vector>
#include <cmath>
#include <gsl/gsl_sum.h>
#include "Logger.hpp"
#include "funcSum.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log(Logger::get_logger("funcSum"));

// --------------------------------------------------------------------------------------------------------------------------------

Real funcSum_all(std::function<Real(uint i)> f, size_t max_i)
{
    const Real p_0(f(0));
    if (p_0 == 0.0) return 0.0;

    Real sum = p_0;
    for (uint i(1); i < max_i; ++i)
        sum += f(i);
    return sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Will simply calculate the sum over a certain function f, until it converges
// (i.e. the sum > tolerance*current_term for a CONVERGENCE_CHECK number of 
// terms), or a maximum number of terms is summed (usually 2000).
//
// Input:
// - f: A       function object
// - max_i:     maximum number of terms it will evaluate
// - tolerance: convergence condition, default value 1-e8 (see .hpp)
Real funcSum(std::function<Real(uint i)> f, size_t max_i, Real tolerance)
{
    const uint CONVERGENCE_CHECK(4);
    const Real p_0(f(0));
    if (p_0 == 0.0) return 0.0;

    RealVector table;
    table.push_back(p_0);

    Real sum = p_0;
    bool extrapolationNeeded(true);
    uint convergenceCounter(0);

    uint i(1);
    for (; i < max_i;)
    {
        const Real p_i(f(i));
        table.push_back(p_i);
        sum += p_i;
        ++i;

        if (fabs(sum) * tolerance >= fabs(p_i)) // '=' is important
            ++convergenceCounter;
        else
            convergenceCounter = 0;

        if (convergenceCounter >= CONVERGENCE_CHECK)
        {
            extrapolationNeeded = false;
            break;
        }
    }

    if (extrapolationNeeded)
    {
        Real error;
        gsl_sum_levin_utrunc_workspace* workspace(gsl_sum_levin_utrunc_alloc(i));
        gsl_sum_levin_utrunc_accel(&table[0], table.size(), workspace, &sum, &error);
        if (fabs(error) >= fabs(sum * tolerance * 10))
            _log.error("series acceleration error: %.16g (rel error: %.16g), terms_used = %d (%d given)", fabs(error), fabs(error / sum), workspace->terms_used, table.size());
        gsl_sum_levin_utrunc_free(workspace);
    }
    return sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

