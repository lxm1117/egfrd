#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE SphericalBesselTable
#include <boost/test/included/unit_test.hpp>
#endif

#include "SphericalBesselGenerator.cpp"
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>


const unsigned int maxnSphericalTable( 50 );
const unsigned int tableResolution( 200 );
const Real TOLERANCE( 1e-4 );

const SphericalBesselGenerator& generator_spherical_table(SphericalBesselGenerator::instance());


BOOST_AUTO_TEST_CASE( testSphericalTableJ )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxnSphericalTable * maxnSphericalTable ) ) * 2 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( i * maxz / resolution );
        

        for( UnsignedInteger n( 0 ); n <= maxnSpherical; ++n )
        {
            const Real tj( generator_spherical_table.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
        }
    }

}
