import unittest
import numpy
import _greens_functions as mod


class GreensFunction1DRadAbsTestCase( unittest.TestCase ):    
    def test_Instantiation( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
	r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )
        self.failIf( gf == None )


    def test_DrawTime( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = 5e-8
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 1 - 1e-8 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_a_equal_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )


    def test_DrawTime_a_near_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-14
        r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_r0_equal_a( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )


    def test_DrawTime_r0_equal_sigma_kf_zero( self ):
        D = 1e-12
	v = -3e-8
        kf = 0.0
        L = 1e-7
        r0 = 0
        gf = mod.GreensFunction1DRadAbs( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t < 0.0 or t >= numpy.inf )


    def no_test_DrawTime_r0_equal_sigma_kf_large( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 20e-7
        r0 = 1e-12
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t < 0.0 or t >= numpy.inf )


    def test_DrawEventType( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = gf.drawTime( 0.5 )
        eventType = gf.drawEventType( 0.999999, t )
        self.assertEqual(eventType, mod.PairEventKind.IV_ESCAPE)

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual(eventType, mod.PairEventKind.IV_REACTION)


    def test_DrawEventType_smallt( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-6
        r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )


        t = gf.drawTime( 0.001 )
        eventType = gf.drawEventType( 0.5, t )
        self.failIf( eventType is not mod.PairEventKind.IV_ESCAPE  
                     and eventType is not mod.PairEventKind.IV_REACTION)

        eventType = gf.drawEventType( 0.9999, t )
        self.assertEqual(eventType, mod.PairEventKind.IV_ESCAPE)

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual(eventType, mod.PairEventKind.IV_REACTION)


    def test_DrawR( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-6
        r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = 1e-3
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        r1 = gf.drawR( 0.00001, t )
        r2 = gf.drawR( 0.99999, t )
        self.assertAlmostEqual( r1, 0, 5)
        self.assertAlmostEqual( r2, L, 5)


    def test_DrawR_zerot( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 1e-7
        r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = 0.0
        r = gf.drawR( 0.5, t )
        self.assertEqual( r0, r )


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = 0
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )
        
        t = 1e-3
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )


    def test_DrawR_squeezed( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-5
	r0 = L/2
        gf = mod.GreensFunction1DRadAbs( D, kf, r0, v, L )

        t = 1e-3
        r0 = 0
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        print "hello r: {}".format(r)
        self.failIf( r < 0 or r > L )

        r0 = 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        r0 = L - 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )
        self.failIf( r < 0 or r > L )


if __name__ == "__main__":
    unittest.main()
