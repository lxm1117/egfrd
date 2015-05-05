import unittest
import _greens_functions as mod
import math
import numpy


class GreensFunction3DRadInfTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_instantiation(self):
        D = 1e-12
        kf = 1e8
        sigma = 1e-8
        a = 1e-7
        r0 = 5e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
        self.failIf(gf == None)


    def test_draw_time(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 5e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = gf.drawTime(0.5)
        self.failIf(t <= 0.0)

        t = gf.drawTime(0.0)
        self.failIf(t < 0.0)

        t = gf.drawTime(0.9999999)
        self.failIf(t <= 0.0)


    def test_draw_time_r0_equal_sigma_kf_zero(self):
        D = 1e-12
        kf = 0.0 # note this
        sigma = 1e-8
        r0 = sigma

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = gf.drawTime(0.5)
        self.failIf(t < 0.0)


    def test_drawR(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 2e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = 1e-3

        r = gf.drawR(0.5, t)
        self.failIf(r < sigma)

        r1 = gf.drawR(0.0, t)
        r2 = gf.drawR(0.9999999, t)

        self.failIf(r1 < sigma)
        self.failIf(r2 < sigma)

        self.failIf(abs(r1 - sigma) > 1e-15)


    def test_drawR_zerot(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 2e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = 0.0

        r = gf.drawR(0.5, t)
        self.assertEqual(r0, r)


    def test_drawR_smallt(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 1.000001e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = 1e-12

        r = gf.drawR(0.5, t)

        self.failIf(r < sigma)


    def test_drawR_r0_equal_sigma(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = sigma

        t = 1e-5

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        r = gf.drawR(0.5, t)
        self.failIf(r < sigma)


    def test_draw_theta(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 2e-8
        t = 1e-3
        r = 2.1e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        theta = gf.drawTheta(0.0, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        theta = gf.drawTheta(0.9999999, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

    '''
    def test_draw_theta2(self):
        D = 2e-12
        kf = 0
        sigma = 5e-9
        r0 = 5.064e-9
        r = 5.05e-9
        t = 1e-9

        gf = mod.GreensFunction3DRadInf(D, kf, sigma)

        theta = gf.drawTheta(0.5, r, r0, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        #theta = gf.drawTheta(0.0, r, r0, t)
        #self.failIf(theta < 0.0 or theta > numpy.pi)
        #theta = gf.drawTheta(0.9999999, r, r0, t)
        #self.failIf(theta < 0.0 or theta > numpy.pi)
'''


    def test_draw_theta_zerot(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r = 5e-8
        r0 = 5e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        t = 0.0
        theta = gf.drawTheta(0.5, r0, t)
        self.assertEqual(0.0, theta)


    def test_draw_theta_smallt(self):

        D = 1e-12
        kf = 1e-8
        sigma = 1e-8

        t = 1e-11

        disp = 3 * math.sqrt(6 * D * t)

        r = sigma + disp + disp
        r0 = sigma + disp

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
        theta = gf.drawTheta(0.5, r, t)

        self.failIf(theta < 0.0 or theta > numpy.pi)


    def test_draw_theta_r0_equal_sigma(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = sigma

        t = 1e-3
        r = r0

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)


    def test_p_int_r_at_s_is_zero(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-3
        r0 = 2e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        pintr = gf.p_int_r(sigma, t)
        self.assertEqual(0.0, pintr)


    def test_p_int_r0_at_s_zerot_is_zero(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 0
        r0 = 2e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        pintr = gf.p_int_r(sigma, t)

        self.assertEqual(0.0, pintr)


    def test_p_int_r_large_is_p_survival(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-3
        r0 = 2e-8

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        pintr = gf.p_int_r(sigma * 1e8, t)
        psurv = gf.p_survival(t)

        self.assertAlmostEqual(psurv, pintr)


    def test_ip_theta_is_int_p_theta(self):

        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-9

        t = 1e-4
        r0 = sigma*1.1

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
        r = r0

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)

        resolution = 10
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution
            ip = gf.ip_theta(theta, r, t)
            result = scipy.integrate.quad(gf.p_theta, 0.0, theta, args=(r, t))
            np = result[0]
            #print 'theta, np, ip', theta, np, ip
            self.assertAlmostEqual(0.0, (np-ip)/ip)


    def test_ip_theta_r0_is_sigma(self):

        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-3
        r0 = sigma

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
        r = 1.1 * r0

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)

        ip = gf.ip_theta(numpy.pi, r, t)
        #print 'ip', ip
        #self.assertEqual(0.0, ip)


    def test_ip_theta_pi_is_p_irr(self):

        D = 1e-12
        sigma = 1e-8

        kf = 0

        t = 1e-3
        r0 = 1.1e-8
        r = r0

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        ip = gf.ip_theta(numpy.pi, r, t) * (2 * numpy.pi * r * r)
        pirr = mod.p_irr(r, t, r0, kf, D, sigma)
        pcorr = gf.ip_corr(numpy.pi, r, t) * (2 * numpy.pi * r * r)
        pfree = gf.ip_free(numpy.pi, r, t) * (2 * numpy.pi * r * r)

        self.assertNotAlmostEqual(pirr, pfree, 6, 'pcorr estimated to be too small. test may not be valid.')

        #print 'PP', pirr, ip, pcorr, pfree

        self.assertNotEqual(0.0, ip)
        self.assertAlmostEqual(ip/pirr, 1)


    def test_ip_theta_pi_at_sigma_is_p_irr(self):

        import math

        D = 1e-12
        sigma = 1e-8

        kf = 0

        t = 1e-5
        r0 = sigma
        r = r0 + math.sqrt(6 * D * t)

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        ip = gf.ip_theta(numpy.pi, r, t) * (2 * numpy.pi * r * r)
        pirr = mod.p_irr(r, t, r0, kf, D, sigma)
        pcorr = gf.ip_corr(numpy.pi, r, t) * (2 * numpy.pi * r * r)
        pfree = gf.ip_free(numpy.pi, r, t) * (2 * numpy.pi * r * r)

        self.assertNotAlmostEqual(pirr, pfree, 7,'pcorr estimated to be too small. test may not be valid.')

        #print 'PP', pirr, ip, pcorr, pfree

        self.assertNotEqual(0.0, ip)
        self.assertAlmostEqual(ip/pirr, 1)


    def test_p_theta_never_negative(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        pint = gf.ip_theta(numpy.pi, r, t)

        pmin = 0.0
        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            p = gf.p_theta(theta, r, t) / pint / resolution
            pmin = min(pmin, p)
            #print 'theta: ', theta, '\tp: ', p

        self.failIf(pmin < 0.0, 'Negative p_theta; t= %g, %s' % (t, gf.dump()))


    def test_ip_theta_never_decrease(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)

        pint_prev = 0.0

        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            pint = gf.ip_theta(theta, r, t)
            #print pint
            self.failIf(pint < pint_prev)
            pint_prev = pint



if __name__ == "__main__":
    unittest.main()
