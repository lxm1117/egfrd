#ifndef CYLINDRICALBESSELGENERATOR_HPP
#define CYLINDRICALBESSELGENERATOR_HPP

#include "Defs.hpp"

class CylindricalBesselGenerator
{
public:
    CylindricalBesselGenerator() {}

    Real J(uint n, Real z) const;
    Real Y(uint n, Real z) const;

    static uint getMinNJ();
    static uint getMinNY();
    static uint getMaxNJ();
    static uint getMaxNY();

    static CylindricalBesselGenerator const& instance();
};

#endif /* CYLINDRICALBESSELGENERATOR_HPP */
