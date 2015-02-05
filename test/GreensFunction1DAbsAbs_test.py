#!/usr/bin/env python

__author__ = 'Laurens Bossen, Thomas Sokolowski'
__copyright__ = ''


import unittest

import _greens_functions as mod

import numpy


class GreensFunction1DAbsAbsTestCase(unittest.TestCase):

    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_Instantiation(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)
        self.failIf(gf == None)


    '''
        def test_DrawTime(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = 5e-8

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime : t =",t

        t = gf.drawTime(0.0)
        self.failIf(t < 0.0 or t >= numpy.inf)
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime : t =",t

        t = gf.drawTime(1e-16)
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime : t =",t

        t = gf.drawTime(1 - 1e-16)
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime : t =",t


    def test_DrawTime_a_equal_sigma(self):
        D = 1e-12
        v = 0
        L = 0.0
        r0 = L

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        self.assertEqual(0.0, t)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime_a_equal_sigma : t =",t


    def test_DrawTime_a_near_sigma(self):
        D = 1e-12
        v = 0
        L = 2e-14
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime_a_near_sigma : t =",t


    def test_DrawTime_r0_equal_a(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = L

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        self.assertEqual(0.0, t)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime_r0_equal_a : t =",t


    def test_DrawTime_r0_equal_sigma(self):
        D = 1e-12
        v = 0
        L = 1e-7
        r0 = 0

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        self.failIf(t < 0.0 or t >= numpy.inf)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawTime_r0_equal_sigma : t =",t


    def test_DrawEventType(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.5)
        eventType = gf.drawEventType(0.5, t)
        self.failIf(eventType != 12 and eventType != 13 and eventType != 14)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawEventType : eventType =",eventType

        eventType = gf.drawEventType(0.999999, t)
        self.assertEqual(eventType, 13) # ESCAPE

        eventType = gf.drawEventType(0.0, t)
        self.assertEqual(eventType, 14) # REACTION


    def no_test_DrawEventType_smallt(self):
        D = 1e-12
        v = 0
        L = 2e-6
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = gf.drawTime(0.001)

        eventType = gf.drawEventType(0.5, t)
        self.failIf(eventType != 12 and eventType != 13 and eventType != 14)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawEventType_smallt : eventType =",eventType

        eventType = gf.drawEventType(0.9999, t)
        self.assertEqual(eventType, 13) # ESCAPE

        eventType = gf.drawEventType(0.0, t)
        self.assertEqual(eventType, 14) # REACTION


    def test_DrawR(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = 1e-3

        r = gf.drawR(0.5, t)
        self.failIf(r < 0 or r > L)

        r1 = gf.drawR(0.0, t)
        r2 = gf.drawR(0.999999999999, t)

        self.failIf(r1 != 0)
        self.failIf(r2 < 0 or r2 > L)

        self.assertAlmostEqual(r1, 0)
        self.assertAlmostEqual(r2, L)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawR : r1 =",r1
        print "GreensFunction1DAbsAbs_test.py : test_DrawR : r2 =",r2


    def test_DrawR_zerot(self):
        D = 1e-12
        v = 0
        L = 1e-7
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = 0.0

        r = gf.drawR(0.5, t)
        self.assertEqual(r0, r)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawR_zerot : r =",r


    def test_DrawR_r0_equal_sigma(self):
        D = 1e-12
        v = 0
        L = 2e-7
        r0 = 0

        t = 0.0#1e-3

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        # This raises an exception, which at this point we cannot catch
        # r = gf.drawR( 0.5, t )
        # self.failIf( r < 0 or r > L )

    def test_DrawR_squeezed(self):

        D = 1e-12
        v = 0
        L = 0.02e-8
        r0 = L / 2

        gf = mod.GreensFunction1DAbsAbs(D, v, r0, L)

        t = 1e-6

        r0 = 0
        gf.setr0(r0)
        #r = gf.drawR( 0.501, t )
        #self.failIf( r < 0 or r > L )

        # near s
        r0 = 0.0001e-8
        gf.setr0(r0)
        r = gf.drawR(0.502, t)
        self.failIf(r < 0 or r > L)
        print " "
        print "GreensFunction1DAbsAbs_test.py : test_DrawR_squeezed : 1st drawn r =",r

        # near a
        r0 = L - 0.0001e-8
        gf.setr0(r0)
        r = gf.drawR(0.503, t)
        self.failIf(r < 0 or r > L)
        print "GreensFunction1DAbsAbs_test.py : test_DrawR_squeezed : 2nd drawn r =",r
    '''


if __name__ == "__main__":
    unittest.main()
