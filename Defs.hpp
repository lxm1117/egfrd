#ifndef DEFS_HPP 
#define DEFS_HPP

#include <vector>
#include <gsl/gsl_math.h>

typedef double Real;
typedef unsigned int uint;
typedef std::vector<Real> RealVector;
typedef std::pair<Real, Real> RealPair;

const uint GF_MAX_ORDER = 50;

#define THROW_UNLESS( CLASS, EXPRESSION )\
    if (!(EXPRESSION)) throw CLASS ("Check ["+std::string(#EXPRESSION) + "] failed.");

const Real SEPARATION_TOLERANCE(1e-07);
const Real MINIMAL_SEPARATION_FACTOR(1.0 + SEPARATION_TOLERANCE);


#ifdef _MSC_VER

#ifdef GF_EXPORT
#define GF_CLASS __declspec(dllexport)
#else
#define GF_CLASS __declspec(dllimport)
#endif

#ifdef GFRD_EXPORT
#define GFRD_CLASS __declspec(dllexport)
#else
#define GFRD_CLASS __declspec(dllimport)
#endif


#else
#define GF_CLASS
#define GFRD_CLASS
#endif



// Template class to map Lambda functions to GSL function pointers
// Uses params as this pointer, so no params avail to functions
template< typename _Lambda>
class gsl_lambda : public gsl_function
{
public:
    gsl_lambda(const _Lambda& func) : _func(func)
    {
        function = &gsl_lambda::invoke;
        params = this;
    }
private:
    const _Lambda& _func;
    static double invoke(double x, void *params)
    {
        return static_cast<gsl_lambda*>(params)->_func(x);
    }
};





#endif // DEFS_HPP
