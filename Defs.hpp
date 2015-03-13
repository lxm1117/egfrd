#ifndef DEFS_HPP 
#define DEFS_HPP

#include <vector>

typedef double Real;
typedef unsigned int uint;
typedef long int Integer;
typedef unsigned long int UnsignedInteger;
typedef std::vector<Real> RealVector;
typedef std::pair<Real, Real> real_pair;
typedef std::vector<Real> RealVector;

const uint GF_MAX_ORDER = 50;

#define THROW_UNLESS( CLASS, EXPRESSION )\
    if (!(EXPRESSION)) throw CLASS ("Check ["+std::string(#EXPRESSION) + "] failed.");

const Real SEPARATION_TOLERANCE(1e-07);
const Real MINIMAL_SEPARATION_FACTOR(1.0 + SEPARATION_TOLERANCE);


#ifndef M_PI
const Real M_PI = 3.1415926535897932384626433832795;
#endif


#define USE_SPHERICALBESSELGENERATOR


#ifdef _MSC_VER
#ifdef GF_EXPORT
#define GF_CLASS __declspec(dllexport)
#else
#define GF_CLASS __declspec(dllimport)
#endif
#else
#define GF_CLASS
#endif


#endif // DEFS_HPP
