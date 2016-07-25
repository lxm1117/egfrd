#if !defined( __FUNCSUM_HPP )
#define __FUNCSUM_HPP

#include <boost/function.hpp>
#include <cstddef>

#include "Defs.hpp"

static const Real TOLERANCE( 1e-8 );

Real funcSum_all(boost::function<Real(unsigned int i)> f, std::size_t max_i);

Real funcSum_all_accel(boost::function<Real(unsigned int i)> f,
                       std::size_t max_i, Real tolerance = TOLERANCE);

Real funcSum(boost::function<Real(unsigned int i)> f,
             std::size_t max_i, Real tolerance = TOLERANCE);


//#define ENABLE_GF_TESTFUNCTIONS
#if defined(ENABLE_GF_TESTFUNCTIONS)
extern "C"
{
   void GF_CLASS TestGreensFunction(unsigned int gfn, unsigned int gfc, unsigned int size, double *buffer, double D, double r0, double a, double kf, double sigma, double t, double r);
   double GF_CLASS DrawGreensFunction(unsigned int gfn, unsigned int gfc, double rnd, double D, double r0, double a, double kf, double sigma, double t, double r);
}
#endif


#endif /* __FUNCSUM_HPP */
