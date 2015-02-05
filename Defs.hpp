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

#define cMAX_ORDER			50

#define THROW_UNLESS( CLASS, EXPRESSION )\
    if (!(EXPRESSION)) throw CLASS ("Check ["+std::string(#EXPRESSION) + "] failed.");

#endif // DEFS_HPP
