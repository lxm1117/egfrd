#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "CylindricalBesselGenerator"
#include <boost/test/included/unit_test.hpp>
#endif

#include "CylindricalBesselGenerator.cpp"
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>


const unsigned int maxnCylindrical(51);
const CylindricalBesselGenerator& generatorCylindrical(CylindricalBesselGenerator::instance());
const Real rel_tol_cylindrical(1e-5);
const Real abs_tol_cylindrical(1e-11);



#define CHECK_ERROR_CYLINDRICAL(n, z, a, b, abs_tol, rel_tol)      \
{\
    const Real abs_error(fabs(a - b));              \
    const Real rel_error(abs_error / fabs(a));\
    BOOST_CHECK_MESSAGE(abs_error < abs_tol_cylindrical || \
                        rel_error < rel_tol_cylindrical, \
                        "n " << n        \
                         << " z " << z                            \
                        << " abs error " << abs_error           \
                         << " rel error " << rel_error);\
}


BOOST_AUTO_TEST_CASE(testCylinderJ)
{
    const UnsignedInteger resolution(300);
    const Real maxz(std::max(1000., static_cast<Real>(maxnCylindrical * maxnCylindrical)) * 2);

    for(UnsignedInteger i(0); i <= resolution; ++i)
    {
        const Real z(maxz * i / resolution);
        
        for(UnsignedInteger n(0); n <= maxnCylindrical; ++n)
        {
            const Real tJ(generatorCylindrical.J(n, z));
            const Real J(gsl_sf_bessel_Jn(n, z));
            
            //BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
            CHECK_ERROR_CYLINDRICAL(n, z, J, tJ, abs_tol_cylindrical, rel_tol_cylindrical);

            //printf("%d %g\n",n,z);
        }
    }

}


BOOST_AUTO_TEST_CASE(testCylinderY)
{
    const UnsignedInteger resolution(300);
    const Real maxz(std::max(1000., static_cast<Real>(maxnCylindrical * maxnCylindrical)) * 2);

    // it is unstable around z==0, so we test for i in [1...resolution]
    for(UnsignedInteger i(1); i <= resolution; ++i)
    {
        const Real z(maxz * i / resolution);
        
        for(UnsignedInteger n(0); n <= maxnCylindrical; ++n)
        {
            const Real tY(generatorCylindrical.Y(n, z));
            const Real Y(gsl_sf_bessel_Yn(n, z));
            
            CHECK_ERROR_CYLINDRICAL(n, z, Y, tY, abs_tol, rel_tol);

            //printf("y %d %g\n",n,z);
        }
    }

}
