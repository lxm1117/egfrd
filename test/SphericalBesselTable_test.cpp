#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE SphericalBesselTable

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>


#include "SphericalBesselGenerator.cpp"


const unsigned int maxn( 50 );
const unsigned int tableResolution( 200 );
const Real TOLERANCE( 1e-4 );

const SphericalBesselGenerator& generator(SphericalBesselGenerator::instance());


BOOST_AUTO_TEST_CASE( testTableJ )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxn * maxn ) ) * 2 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( i * maxz / resolution );
        

        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real tj( generator.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
        }
    }

}
