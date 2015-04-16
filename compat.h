#if !defined( COMPAT_HPP )
#define COMPAT_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include <limits>

#if !defined( HAVE_DECL_INFINITY )
#if !defined( INFINITY )
#    define INFINITY ( std::numeric_limits< double >::infinity() )
#endif
#endif /* HAVE_DECL_INFINITY */

#if !defined( HAVE_SINCOS )
inline void sincos( double x, double* s, double* c )
{
    *s = sin( x );
    *c = cos( x );
}
#endif /* !HAVE_SINCOS */

#endif // __COMPAT_HPP
