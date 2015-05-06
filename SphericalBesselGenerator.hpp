#ifndef SPHERICALBESSELGENERATOR_HPP
#define SPHERICALBESSELGENERATOR_HPP

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

#endif /* SPHERICALBESSELGENERATOR_HPP */
