#ifndef __SPHERICALBESSELGENERATOR_HPP
#define __SPHERICALBESSELGENERATOR_HPP

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include "Defs.hpp"

class SphericalBesselGenerator
{
public:
    SphericalBesselGenerator() { }

    Real j(uint n, Real z) const;
    Real y(uint n, Real z) const;

    static uint getMinNJ();
    static uint getMinNY();
    static uint getMaxNJ();
    static uint getMaxNY();

    static SphericalBesselGenerator const& instance();
};

#endif /* __SPHERICALBESSELGENERATOR_HPP */
