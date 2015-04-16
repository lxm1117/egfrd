#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "SphericalBesselGenerator"
#include <boost/test/included/unit_test.hpp>
#endif

#include "SphericalBesselGenerator.hpp"
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>

const unsigned int maxnSpherical( 51 );
const SphericalBesselGenerator& generatorSpherical(SphericalBesselGenerator::instance());

const Real rel_tol_spherical( 1e-5 );
const Real abs_tol_spherical( 1e-11 );


#define CHECK_ERROR_SPHERICAL( n, z, a, b, abs_tol, rel_tol )\
{\
    const Real abs_error( fabs( a - b ) );              \
    const Real rel_error( abs_error / fabs( a ) );\
    BOOST_CHECK_MESSAGE( abs_error < abs_tol_spherical || \
                         rel_error < rel_tol_spherical, \
                        "n " << n        \
                         << " z " << z                            \
                        << " abs error " << abs_error           \
                         << " rel error " << rel_error );\
}


BOOST_AUTO_TEST_CASE( testSphericalJ )
{
    const UnsignedInteger resolution( 300 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxnSpherical * maxnSpherical ) ) * 2 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxnSpherical; ++n )
        {
	  const Real tj(generatorSpherical.j(n, z));
          const Real j( gsl_sf_bessel_jl( n, z));
            
          //BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
          CHECK_ERROR_SPHERICAL( n, z, j, tj, abs_tol, rel_tol );

          //printf("%d %g\n",n,z);
        }
    }

}


BOOST_AUTO_TEST_CASE( testSphericalY )
{
    const UnsignedInteger resolution( 300 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxnSpherical * maxnSpherical ) ) * 2 );

    // it is unstable around z==0, so we test for i in [1...resolution]
    for( UnsignedInteger i( 1 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxnSpherical; ++n )
        {
            const Real ty( generatorSpherical.y( n, z ) );
            const Real y( gsl_sf_bessel_yl( n, z ) );
            
            CHECK_ERROR_SPHERICAL( n, z, y, ty, abs_tol, rel_tol );
            //printf("y %d %g\n",n,z);
        }
    }

    }
