import unittest
import FKIAMToolbox as fkt
import numpy as np
import concurrent.futures
import time
from scipy import misc


def default_parameters():

    fkt.equationsmodule.stm_required = 0

    fkt.equationsmodule.atm = 0
    fkt.equationsmodule.j2 = 0
    fkt.equationsmodule.srp = 0
    fkt.equationsmodule.sun = 0
    fkt.equationsmodule.mercury = 0
    fkt.equationsmodule.venus = 0
    fkt.equationsmodule.earth = 0
    fkt.equationsmodule.mars = 0
    fkt.equationsmodule.jupiter = 0
    fkt.equationsmodule.saturn = 0
    fkt.equationsmodule.uranus = 0
    fkt.equationsmodule.neptune = 0
    fkt.equationsmodule.moon = 0
    fkt.equationsmodule.cmplxmoon = 0

    fkt.equationsmodule.jd_zero = 2459599.5
    fkt.equationsmodule.order = 50
    fkt.equationsmodule.area = 1
    fkt.equationsmodule.mass = 100


def hard_parameters():

    fkt.equationsmodule.stm_required = 1

    fkt.equationsmodule.atm = 1
    fkt.equationsmodule.j2 = 1
    fkt.equationsmodule.srp = 1
    fkt.equationsmodule.sun = 1
    fkt.equationsmodule.mercury = 1
    fkt.equationsmodule.venus = 1
    fkt.equationsmodule.earth = 1
    fkt.equationsmodule.mars = 1
    fkt.equationsmodule.jupiter = 1
    fkt.equationsmodule.saturn = 1
    fkt.equationsmodule.uranus = 1
    fkt.equationsmodule.neptune = 1
    fkt.equationsmodule.moon = 1
    fkt.equationsmodule.cmplxmoon = 1

    fkt.equationsmodule.jd_zero = 2459599.5
    fkt.equationsmodule.order = 50
    fkt.equationsmodule.area = 1
    fkt.equationsmodule.mass = 100


def propagate(central_body, t0, tf, x0, sources, dat, stm, variables):
    neq = 42 if stm else 6
    t, y = fkt.propagationmodule.propagatetrajectory(central_body, t0, tf, x0, sources, dat,
                                                     stm, variables, neq)
    return t, y


def rv2ee(rv, mu, grad_req):
    ee, _ = fkt.translations.krv2ee(rv, mu, grad_req)
    return ee


def ee2rv(ee, mu, grad_req):
    rv, _ = fkt.translations.kee2rv(ee, mu, grad_req)
    return rv


def rv2oe(rv, mu, grad_req):
    oe, _ = fkt.translations.krv2oe(rv, mu, grad_req)
    return oe


def oe2rv(oe, mu, grad_req):
    rv, _ = fkt.translations.koe2rv(oe, mu, grad_req)
    return rv


def derivative(func, x0, order):
    f0 = func(x0)
    eps = 1.0e-5
    dfdx = np.zeros((len(f0), len(x0)))
    for i in range(len(x0)):
        xp = x0.copy()
        xm = x0.copy()
        h = eps
        if order == 1:
            xp[i] = xp[i] + h
            fp = func(xp)
            dfdx[:, i] = (fp - f0) / h
        elif order == 2:
            xp[i] = xp[i] + h
            xm[i] = xm[i] - h
            fp = func(xp)
            fm = func(xm)
            dfdx[:, i] = (fp - fm) / h / 2.0
        elif order == 4:
            xm2 = x0.copy()
            xp2 = x0.copy()
            xm2[i] = xm2[i] - 2.0*h
            xm[i] = xm[i] - h
            xp[i] = xp[i] + h
            xp2[i] = xp2[i] + 2.0*h
            fm2 = func(xm2)
            fm = func(xm)
            fp = func(xp)
            fp2 = func(xp2)
            dfdx[:, i] = ((1/12)*fm2 - (2/3)*fm + (2/3)*fp - (1/12)*fp2)/h
    return dfdx


class Test(unittest.TestCase):

    def test_abouts(self):
        print('Running test_abouts...')
        default_parameters()
        result1 = fkt.basemeanstoolbox.basemeanstoolbox_about()
        result2 = fkt.ephemeris.ephemeris_about()
        result3 = fkt.equationsmodule.equationsmodule_about()
        result4 = fkt.linearalgebrainterfaces.linearalgebrainterfaces_about()
        result5 = fkt.linearalgebralowlevel.linearalgebralowlevel_about()
        result6 = fkt.odetoolbox.odetoolbox_about()
        result7 = fkt.translations.translations_about()
        self.assertEqual(result1, 0)
        self.assertEqual(result2, 0)
        self.assertEqual(result3, 0)
        self.assertEqual(result4, 0)
        self.assertEqual(result5, 0)
        self.assertEqual(result6, 0)
        self.assertEqual(result7, 0)
        print('Done test_abouts.')

    def test_constants_and_units(self):
        earth_equatorradius = fkt.constantsandunits.earth_equatorradius
        self.assertEqual(earth_equatorradius, 6378.1366)
        moon_units = fkt.constantsandunits.kunits_onebody('Moon')
        self.assertEqual(moon_units[0], 4902.800145956161)
        self.assertEqual(moon_units[1], 1737.4)
        self.assertEqual(moon_units[2], 1.679856509416075)
        self.assertEqual(moon_units[3], 0.01197054402181422)
        self.assertEqual(moon_units[4], 1.6242188858222395)
        earth_moon_units = fkt.constantsandunits.kunits_twobody('Earth', 'Moon')
        self.assertEqual(earth_moon_units[0], 403503.2357459562)
        self.assertEqual(earth_moon_units[1], 0.012150584460351)
        self.assertEqual(earth_moon_units[2], 384402.0)
        self.assertEqual(earth_moon_units[3], 1.0245441822513066)
        self.assertEqual(earth_moon_units[4], 4.342513772754916)
        self.assertEqual(earth_moon_units[5], 0.002730711030080485)

    def test_juliandate(self):
        print('Running test_juliandate...')
        default_parameters()
        result = fkt.ephemeris.juliandate(2022, 1, 20, 0, 0, 0)
        self.assertEqual(result, 2459599.5)
        print('Done test_juliandate.')

    def test_planetstate(self):

        print('Running test_planetstate...')

        default_parameters()

        result_earth_moon = fkt.ephemeris.planetstate(2459599.5, 'Moon', 'Earth')
        result_earth_moon_true = np.array([-3.119060309811642e+05,  2.071701572294455e+05,  1.2755218784496e+05,
                                           -5.8072501044365e-01, -7.474093014705249e-01, -3.1550783290824264e-01])
        np.testing.assert_array_equal(result_earth_moon, result_earth_moon_true)

        result_earth_sun = fkt.ephemeris.planetstate(2459599.5, 'Sun', 'Earth')
        result_earth_sun_true = np.array([7.267933677076249e+07, -1.1744076654990216e+08, -5.090989948165558e+07,
                                          2.638118621411191e+01, 1.359034418533267e+01,  5.891320964254108e+00])
        np.testing.assert_array_equal(result_earth_sun, result_earth_sun_true)

        result_earth_jupiter = fkt.ephemeris.planetstate(2459599.5, 'Jupiter', 'Earth')
        result_earth_jupiter_true = np.array([7.766805566345683e+08, -3.3848215372106415e+08, -1.6279092184591967e+08,
                                              3.0560182587866525e+01, 2.5534940593655627e+01,  1.0909383877203007e+01])
        np.testing.assert_array_equal(result_earth_jupiter, result_earth_jupiter_true)

        print('Done test_planetstate.')

    def test_moonlibration(self):
        print('Running test_moonlibration...')
        default_parameters()
        result = fkt.ephemeris.moonlibration(2459599.5)
        result_true = np.array([-5.9797312187837015e-02,  3.9576132531222225e-01,  4.416564077692134e+03])
        np.testing.assert_array_equal(result, result_true)
        print('Done test_moonlibration.')

    def test_r2bp_equations(self):

        print('Running test_r2bp_equations...')

        default_parameters()

        result = fkt.equationsmodule.kr2bp(0, np.array([1, 0, 0, 0, 1, 0]))
        result_true = np.array([0.,  1.,  0., -1., -0., -0.])
        np.testing.assert_array_equal(result_true, result)

        for _ in range(10):
            x = np.random.randn(6)
            result = fkt.equationsmodule.kr2bp(0, x)
            result_true = np.concatenate([x[3:6], -x[0:3]/np.linalg.norm(x[0:3])**3])
            np.testing.assert_array_equal(result_true, result)

        result_rv_earth = fkt.equationsmodule.knbp_rv_earth(0, np.array([1, 0, 0, 0, 1, 0]))
        result_rv_moon = fkt.equationsmodule.knbp_rv_moon(0, np.array([1, 0, 0, 0, 1, 0]))
        result_rv_true = np.array([0.,  1.,  0., -1., -0., -0.])
        np.testing.assert_array_equal(result_rv_earth, result_rv_true)
        np.testing.assert_array_equal(result_rv_moon, result_rv_true)

        result_ee_earth = fkt.equationsmodule.knbp_ee_earth(0, np.array([1, 0, 0, 0, 0, 0]))
        result_ee_moon = fkt.equationsmodule.knbp_ee_moon(0, np.array([1, 0, 0, 0, 0, 0]))
        result_ee_true = np.array([0., 0., 0., 0., 0., 1.])
        np.testing.assert_array_equal(result_ee_earth, result_ee_true)
        np.testing.assert_array_equal(result_ee_moon, result_ee_true)

        print('Done test_r2bp_equations.')

    def test_cr3bp_equations(self):

        atol = 1.0e-14
        rtol = 1.0e-14

        for _ in range(100):
            mu = np.random.random()
            fkt.equationsmodule.massparameter = mu
            x = np.random.randn(42)
            fc = fkt.equationsmodule.kcr3bp(0.0, x)
            ff = fkt.equationsmodule.kcr3bp_fb(0.0, x + np.concatenate([[mu], np.zeros(41)]))
            fs = fkt.equationsmodule.kcr3bp_sb(0.0, x + np.concatenate([[mu - 1], np.zeros(41)]))
            np.testing.assert_allclose(fc, ff, rtol=rtol, atol=atol)
            np.testing.assert_allclose(fc, fs, rtol=rtol, atol=atol)
            np.testing.assert_allclose(ff, fs, rtol=rtol, atol=atol)

    def test_khill_kbr4bp(self):

        atol = 1.0e-14
        rtol = 1.0e-14

        for _ in range(10):
            x = np.random.randn(6)
            f = fkt.equationsmodule.khill(0.0, x)

        for _ in range(100):
            mu = np.random.random()
            fkt.equationsmodule.massparameter = mu
            x = np.random.randn(42)
            fc = fkt.equationsmodule.kbr4bp(0.0, x)
            ff = fkt.equationsmodule.kbr4bp_fb(0.0, x + np.concatenate([[mu], np.zeros(41)]))
            fs = fkt.equationsmodule.kbr4bp_sb(0.0, x + np.concatenate([[mu - 1], np.zeros(41)]))
            np.testing.assert_allclose(fc, ff, rtol=rtol, atol=atol)
            np.testing.assert_allclose(fc, fs, rtol=rtol, atol=atol)
            np.testing.assert_allclose(ff, fs, rtol=rtol, atol=atol)

    def test_nbp_rv_equations(self):

        print('Running test_nbp_rv_equations...')

        hard_parameters()

        fkt.equationsmodule.distunit = fkt.equationsmodule.earthdistunit
        fkt.equationsmodule.velunit = fkt.equationsmodule.earthvelunit
        fkt.equationsmodule.timeunit = fkt.equationsmodule.earthtimeunit
        fkt.equationsmodule.accunit = fkt.equationsmodule.earthaccunit
        fkt.equationsmodule.rsun = fkt.equationsmodule.earthrsun
        fkt.equationsmodule.rearth = fkt.equationsmodule.earthrearth
        fkt.equationsmodule.rmoon = fkt.equationsmodule.earthrmoon

        s = np.append(np.array([1.1, 0, 0, 0, 0.9, 0]), np.reshape(np.eye(6), (36,)))
        result_rv_earth = fkt.equationsmodule.knbp_rv_earth(0, s)
        result_rv_true = np.array([0.0, 0.9, 0.0, -0.8275566130556367,
                                   -1.0601077609733587e-07, -4.751264387875831e-06,
                                   0.0, 0.0, 0.0, 1.5066673495447398,
                                   4.3461165073655047e-07, 1.702368614746311e-05,
                                   0.0, 0.0, 0.0, -9.231608247890953e-08,
                                   -0.7523242048572218, 1.026468166409056e-08,
                                   0.0, 0.0, 0.0, 1.7023686147463116e-05,
                                   1.0264681664091692e-08, -0.7543431446875182,
                                   1.0, 0.0, 0.0, -6.952639388061896e-09,
                                   0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                   -1.3905278776123791e-08, 0.0, 0.0, 0.0, 1.0, 0.0,
                                   0.0, -6.952639388061896e-09])
        np.testing.assert_array_equal(result_rv_earth, result_rv_true)

        fkt.equationsmodule.distunit = fkt.equationsmodule.moondistunit
        fkt.equationsmodule.velunit = fkt.equationsmodule.moonvelunit
        fkt.equationsmodule.timeunit = fkt.equationsmodule.moontimeunit
        fkt.equationsmodule.accunit = fkt.equationsmodule.moonaccunit
        fkt.equationsmodule.rsun = fkt.equationsmodule.moonrsun
        fkt.equationsmodule.rearth = fkt.equationsmodule.moonrearth
        fkt.equationsmodule.rmoon = fkt.equationsmodule.moonrmoon

        result_rv_moon = fkt.equationsmodule.knbp_rv_moon(0, s)
        result_rv_true = np.array([0.,  0.9,  0., -0.8265190256349302,
                                   -0.00014936963116617076, -0.00017661107892690553,
                                   0.0, 0.0, 0.0, 1.5019976897487333,
                                   0.0008142002702310721, 0.0014755766187651213,
                                   0.0, 0.0, 0.0, 0.000814200270231209,
                                   -0.7511528073153859, 5.570575324783031e-05,
                                   0.0, 0.0, 0.0, 0.0014755766187652215,
                                   5.57057532479245e-05, -0.7508448824333469,
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0, 0.0])
        np.testing.assert_array_equal(result_rv_moon, result_rv_true)

        print('Done test_nbp_rv_equations.')

    def test_nbp_ee_equations(self):

        print('Running test_nbp_ee_equations...')

        hard_parameters()

        fkt.equationsmodule.stm_required = 0

        fkt.equationsmodule.distunit = fkt.equationsmodule.earthdistunit
        fkt.equationsmodule.velunit = fkt.equationsmodule.earthvelunit
        fkt.equationsmodule.timeunit = fkt.equationsmodule.earthtimeunit
        fkt.equationsmodule.accunit = fkt.equationsmodule.earthaccunit
        fkt.equationsmodule.rsun = fkt.equationsmodule.earthrsun
        fkt.equationsmodule.rearth = fkt.equationsmodule.earthrearth
        fkt.equationsmodule.rmoon = fkt.equationsmodule.earthrmoon

        s = np.array([1.1, 0, 0, 0, 0.0, 0])
        result_ee_earth = fkt.equationsmodule.knbp_ee_earth(0, s)
        result_ee_true = np.array([-1.329782321454122e-07, -2.417786039007494e-07,
                                   0.0008341910544965359, -1.7976458125947337e-06,
                                   0.0, 0.7513148009015775])
        np.testing.assert_array_equal(result_ee_earth, result_ee_true)

        fkt.equationsmodule.distunit = fkt.equationsmodule.moondistunit
        fkt.equationsmodule.velunit = fkt.equationsmodule.moonvelunit
        fkt.equationsmodule.timeunit = fkt.equationsmodule.moontimeunit
        fkt.equationsmodule.accunit = fkt.equationsmodule.moonaccunit
        fkt.equationsmodule.rsun = fkt.equationsmodule.moonrsun
        fkt.equationsmodule.rearth = fkt.equationsmodule.moonrearth
        fkt.equationsmodule.rmoon = fkt.equationsmodule.moonrmoon

        result_ee_moon = fkt.equationsmodule.knbp_ee_moon(0, s)
        result_ee_true = np.array([-0.00011237687363829448, -0.00020432158843326268,
                                   9.778568352127871e-05, -4.363777563548413e-05,
                                   0.0, 0.7513148009015775])
        np.testing.assert_array_equal(result_ee_moon, result_ee_true)

        print('Done test_nbp_ee_equations.')

    def test_r2bp_rv_propagation(self):

        print('Running test_r2bp_rv_propagation...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2*np.pi
        x0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        sources = np.zeros((14,))
        dat = np.array([2459599.5, 1, 100, 1])
        stm = False
        variables = 'rv'

        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([0.9999999999385388, 8.914497825258285e-10, 0.0,
                            -1.126415049339647e-09, 0.9999999998739921, 0.0])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print(finish - start)

        central_body = 'Moon'
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print('Done test_r2bp_rv_propagation.')

    def test_r2bp_rv_propagation_in_parallel(self):

        print('Running test_r2bp_rv_propagation_in_parallel...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        sources = np.zeros((14,))
        dat = np.array([2459599.5, 1, 100, 1])
        stm = False
        variables = 'rv'

        results = []
        exec_list = []
        start = time.perf_counter()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for _ in range(100):
                x0 = np.array([1.0 + np.random.random(), np.random.random(), np.random.random() / 2, np.random.random(),
                               np.random.random(), np.random.random()])
                exec_list.append(executor.submit(propagate, central_body, t0, tf, x0, sources, dat, stm, variables))
            for f in concurrent.futures.as_completed(exec_list):
                results.append(f.result())
        finish = time.perf_counter()
        execution_time = finish - start

        print(f'Parallel execution time: {execution_time} seconds.')

        print('Done test_r2bp_rv_propagation_in_parallel.')

    def test_nbp_rv_propagation(self):

        print('Running test_nbp_rv_propagation...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        x0 = np.array([1.1, 0.0, 0.0, 0.0, 1/np.sqrt(1.1), 0.0])
        sources = np.ones((14,))
        dat = np.array([2459599.5, 1, 100, 50])
        stm = True
        variables = 'rv'

        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([0.749822681194083, -0.8041840415461319, 1.2566453609966705e-05,
                            0.698305751055437, 0.6498103106191498, -8.29360591644147e-06,
                            -12.17261013349227, -13.137142426818757, 0.00012339471911778277,
                            9.932995840165056, -11.067067998376201, 0.00026476346586438277,
                            -0.24610566853663365, 0.7787099782275246, -1.610460636531204e-05,
                            0.660674059784143, 0.22991629883192713, 1.0003309225900676e-06,
                            -0.00012385957765436467, -0.0001117590594340977, 0.6870089560890327,
                            9.313315299294556e-05, -9.470981920066025e-05, 0.6313123050467405,
                            -1.1273628183697881, 0.1119709306108691, -3.7257032247975065e-06,
                            1.443737871056645, -0.4671355198172131, 1.2609583721105881e-05,
                            -15.799487820524556, -15.725035752549203, 0.00014626329165283314,
                            12.443190031932799, -12.83078722666189, 0.00031651446442243134,
                            -3.448747881772827e-06, -1.5622223294322682e-05, -0.8365123950984958,
                            5.76618456705178e-06, 7.414934437045781e-06, 0.6868897073412237])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print(finish - start)

        central_body = 'Moon'
        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([0.738862083492851, -0.8145675540123463, -0.0002750051477969533,
                            0.7066357757835809, 0.640544374036356, -0.0004372173570709712,
                            -12.317309477988225, -12.961547939702546, 0.0072969319903988935,
                            9.744661006226456, -11.17135739808601, 0.0007649559182500601,
                            -0.24657148198216555, 0.7771167164229418, 0.0003880907224057681,
                            0.6691381455825598, 0.2192811582427938, 4.0581979226116456e-05,
                            -0.003927953190322167, -0.005149268602794023, 0.6724183851611704,
                            0.003121382649441994, -0.004129912731157314, 0.6412511651084037,
                            -1.1364923277837544, 0.12430378233330772, 0.0011587591747360749,
                            1.4420460333618805, -0.48595695197695654, -0.00023163499254907213,
                            -15.99829292423572, -15.529941421156733, 0.008393317462763882,
                            12.238530353856325, -12.98751893822642, 0.0011581457241120689,
                            0.0022522122050754666, 0.001519094838514473, -0.8541531640723331,
                            -0.0015322671548077158, 0.0015012186974713165, 0.6726106337158614])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print('Done test_nbp_rv_propagation.')

    def test_nbp_rv_propagation_in_parallel(self):

        print('Running test_nbp_rv_propagation_in_parallel...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        sources = np.ones((14,))
        dat = np.array([2459565.5, 2, 100, 50])
        stm = True
        variables = 'rv'

        results = []
        exec_list = []
        start = time.perf_counter()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for _ in range(5):
                x0 = np.array([1.0 + np.random.random(), np.random.random(), np.random.random() / 2, np.random.random(),
                               np.random.random(), np.random.random()])
                exec_list.append(executor.submit(propagate, central_body, t0, tf, x0, sources, dat, stm, variables))
            for f in concurrent.futures.as_completed(exec_list):
                results.append(f.result())
        finish = time.perf_counter()
        execution_time = finish - start

        print(f'Parallel execution time: {execution_time} seconds.')

        central_body = 'Moon'
        stm = False
        results = []
        exec_list = []
        tf = 2 * np.pi * 5
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for _ in range(4):
                x0 = np.array([1.5, 0.0, 0.0, 0.0, 1/np.sqrt(1.5), 0.0])
                exec_list.append(executor.submit(propagate, central_body, t0, tf, x0, sources, dat, stm, variables))
            for f in concurrent.futures.as_completed(exec_list):
                results.append(f.result())
        finish = time.perf_counter()
        execution_time = finish - start

        print(f'Parallel execution time: {execution_time} seconds.')

        print('Done test_r2bp_rv_propagation_in_parallel.')

    def test_r2bp_ee_propagation(self):

        print('Running test_nbp_rv_propagation...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        x0 = np.array([1.1, 0.0, 0.0, 0.0, 0.0, 0.0])
        sources = np.zeros((14,))
        dat = np.array([2459599.5, 1, 100, 1])
        stm = False
        variables = 'ee'

        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([1.1, 0.0, 0.0, 0.0, 0.0, 4.720650118091347])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print(finish - start)

        central_body = 'Moon'
        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([1.1, 0.0, 0.0, 0.0, 0.0, 4.720650118091347])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print('Done test_nbp_rv_propagation.')

    def test_nbp_ee_propagation(self):

        print('Running test_nbp_rv_propagation...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        x0 = np.array([1.1, 0.0, 0.0, 0.0, 0.0, 0.0])
        sources = np.ones((14,))
        dat = np.array([2459599.5, 1, 100, 50])
        stm = False
        variables = 'ee'

        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([1.0999999853526852, -0.0010879729465415912, -0.0011181220282100541,
                            -5.636709668616053e-06, -1.1281531034928926e-06, 4.733379050591152])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print(finish - start)

        central_body = 'Moon'
        start = time.perf_counter()
        t, y = propagate(central_body, t0, tf, x0, sources, dat, stm, variables)
        finish = time.perf_counter()
        t_true = 6.283185307179586
        yf_true = np.array([1.0999254910284573, -0.00022278776815805909, -0.00022095866374405247,
                            -3.913900279475531e-05, 0.00019286581630369965, 4.722842072934673])
        self.assertEqual(t[1], t_true)
        np.testing.assert_array_equal(y[:, 1], yf_true)

        print('Done test_nbp_rv_propagation.')

    def test_nbp_ee_propagation_in_parallel(self):

        print('Running test_nbp_ee_propagation_in_parallel...')

        central_body = 'Earth'
        t0 = 0.0
        tf = 2 * np.pi
        sources = np.ones((14,))
        dat = np.array([2459565.5, 2, 100, 50])
        stm = False
        variables = 'ee'

        results = []
        exec_list = []
        start = time.perf_counter()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for _ in range(5):
                x0 = np.array([1.0 + np.random.random(), np.random.random()/2, np.random.random()/2, np.random.random()/2,
                               np.random.random()/2, np.random.random()])
                exec_list.append(executor.submit(propagate, central_body, t0, tf, x0, sources, dat, stm, variables))
            for f in concurrent.futures.as_completed(exec_list):
                results.append(f.result())
        finish = time.perf_counter()
        execution_time = finish - start

        print(f'Parallel execution time: {execution_time} seconds.')

        central_body = 'Moon'
        stm = False
        results = []
        exec_list = []
        tf = 2 * np.pi * 5
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for _ in range(4):
                x0 = np.array([1.5, 0.0, 0.0, 0.0, 1/np.sqrt(1.5), 0.0])
                exec_list.append(executor.submit(propagate, central_body, t0, tf, x0, sources, dat, stm, variables))
            for f in concurrent.futures.as_completed(exec_list):
                results.append(f.result())
        finish = time.perf_counter()
        execution_time = finish - start

        print(f'Parallel execution time: {execution_time} seconds.')

        print('Done test_r2bp_rv_propagation_in_parallel.')

    def test_rv2ee_ee2rv(self):

        rtol = 5.0e-8
        atol = 5.0e-8

        rv0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        ee0, _ = fkt.translations.krv2ee(rv0, 1.0, True)
        np.testing.assert_array_equal(ee0, np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        rv1, drv1 = fkt.translations.kee2rv(ee0, 1.0, True)
        np.testing.assert_array_equal(rv1, rv0)

        for _ in range(100):
            rv0 = np.random.random(6)
            mu = 1.0 + np.minimum([np.random.randn()/5.0], [0.7])
            ee0, dee0 = fkt.translations.krv2ee(rv0, mu, True)
            rv1, drv1 = fkt.translations.kee2rv(ee0, mu, True)
            dee0_true = derivative(lambda x: rv2ee(x, mu, False), rv0, 4)
            drv1_true = derivative(lambda x: ee2rv(x, mu, False), ee0, 4)
            np.testing.assert_allclose(rv1, rv0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(dee0, drv1), np.eye(6), rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(drv1, dee0), np.eye(6), rtol=rtol, atol=atol)
            np.testing.assert_allclose(dee0_true, dee0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(drv1_true, drv1, rtol=rtol, atol=atol)
            ee1, dee1 = fkt.translations.krv2ee(rv1, mu, True)
            np.testing.assert_allclose(ee1, ee0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dee1, dee0, rtol=rtol, atol=atol)
            rv2, drv2 = fkt.translations.kee2rv(ee1, mu, True)
            np.testing.assert_allclose(rv2, rv1, rtol=rtol, atol=atol)
            np.testing.assert_allclose(drv2, drv1, rtol=rtol, atol=atol)

    def test_rv2oe_oe2rv(self):

        rtol = 5.0e-7
        atol = 5.0e-6

        rv0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.1])
        oe0, _ = fkt.translations.krv2oe(rv0, 1.0, True)
        np.testing.assert_allclose(oe0, np.array([1.0101010101010102,
                                                  0.009999999999999449,
                                                  0.09966865249116073,
                                                  0.0, 0.0, 0.0]), rtol=1e-13, atol=1e-13)
        rv1, drv1 = fkt.translations.koe2rv(oe0, 1.0, True)
        np.testing.assert_allclose(rv1, rv0, rtol=1e-13, atol=1e-13)

        for _ in range(100):
            rv0 = np.random.random(6)
            mu = 1.0 + np.minimum([np.random.randn() / 5.0], [0.7])
            oe0, doe0 = fkt.translations.krv2oe(rv0, mu, True)
            rv1, drv1 = fkt.translations.koe2rv(oe0, mu, True)
            doe0_true = derivative(lambda x: rv2oe(x, mu, False), rv0, 4)
            drv1_true = derivative(lambda x: oe2rv(x, mu, False), oe0, 4)
            np.testing.assert_allclose(rv1, rv0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(doe0, drv1), np.eye(6), rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(drv1, doe0), np.eye(6), rtol=rtol, atol=atol)
            np.testing.assert_allclose(doe0_true, doe0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(drv1_true, drv1, rtol=rtol, atol=atol)
            oe1, doe1 = fkt.translations.krv2oe(rv1, mu, True)
            np.testing.assert_allclose(oe1, oe0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(doe1, doe0, rtol=rtol, atol=atol)
            rv2, drv2 = fkt.translations.koe2rv(oe1, mu, True)
            np.testing.assert_allclose(rv2, rv1, rtol=rtol, atol=atol)
            np.testing.assert_allclose(drv2, drv1, rtol=rtol, atol=atol)

    def test_cart2sphere_sphere2cart(self):

        x0 = np.array([[0, 0, 1]]).T
        s0 = fkt.translations.kcart2sphere(x0)
        np.testing.assert_equal(s0, np.array([[1, 0, 0]]).T)
        x1 = fkt.translations.ksphere2cart(s0)
        np.testing.assert_equal(x1, x0)
        s1 = fkt.translations.kcart2sphere(x1)
        np.testing.assert_equal(s1, s0)

        x = np.array([[1, 0, 0]]).T
        s = fkt.translations.kcart2sphere(x)
        np.testing.assert_equal(s, np.array([[1, 0, np.pi/2]]).T)

    def test_cart2latlon_latlon2cart(self):

        x0 = np.array([[1, 0, 0]]).T
        s0 = fkt.translations.kcart2latlon(x0)
        np.testing.assert_equal(s0, np.array([[0, 0]]).T)
        x1 = fkt.translations.klatlon2cart(s0)
        np.testing.assert_equal(x1, x0)
        s1 = fkt.translations.kcart2latlon(x1)
        np.testing.assert_equal(s1, s0)

    def test_kitrs2gcrs_kgcrs2itrs(self):

        rtol = 1.0e-13
        atol = 1.0e-13

        for _ in range(100):
            xitrs0 = np.random.randn(3)
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13), np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xgcrs0, dxgcrs0 = fkt.translations.kitrs2gcrs(xitrs0, jd)
            xitrs1, dxitrs1 = fkt.translations.kgcrs2itrs(xgcrs0, jd)
            np.testing.assert_allclose(xitrs1, xitrs0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxgcrs0, np.transpose(dxitrs1), rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(dxgcrs0, dxitrs1), np.eye(3), rtol=rtol, atol=atol)
            xgcrs1, dxgcrs1 = fkt.translations.kitrs2gcrs(xitrs1, jd)
            np.testing.assert_allclose(xgcrs1, xgcrs0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxgcrs1, dxgcrs0, rtol=rtol, atol=atol)
            xitrs2, dxitrs2 = fkt.translations.kgcrs2itrs(xgcrs1, jd)
            np.testing.assert_allclose(xitrs2, xitrs1, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxitrs2, dxitrs1, rtol=rtol, atol=atol)

    def test_kscrs2mer_kmer2scrs(self):

        rtol = 1.0e-13
        atol = 1.0e-13

        for _ in range(100):
            xscrs0 = np.random.randn(3, 1)
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13),
                                          np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xmer0, dxmer0 = fkt.translations.kscrs2mer(xscrs0, jd)
            xscrs1, dxscrs1 = fkt.translations.kmer2scrs(xmer0, jd)
            np.testing.assert_allclose(xscrs1, xscrs0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxmer0[:, :, 0], np.transpose(dxscrs1[:, :, 0]), rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(dxmer0[:, :, 0], dxscrs1[:, :, 0]), np.eye(3), rtol=rtol, atol=atol)
            xmer1, dxmer1 = fkt.translations.kscrs2mer(xscrs1, jd)
            np.testing.assert_allclose(xmer1, xmer0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxmer1, dxmer0, rtol=rtol, atol=atol)

    def test_kscrs2gcrs_kgcrs2scrs(self):

        rtol = 1.0e-10
        atol = 1.0e-10

        for _ in range(100):
            xscrs0 = np.random.randn(6, 1)
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13),
                                          np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xgcrs0 = fkt.translations.kscrs2gcrs(xscrs0, jd, 1.0, 1.0)
            xscrs1 = fkt.translations.kgcrs2scrs(xgcrs0, jd, 1.0, 1.0)
            np.testing.assert_allclose(xscrs1, xscrs0, rtol=rtol, atol=atol)
            xgcrs1 = fkt.translations.kscrs2gcrs(xscrs1, jd, 1.0, 1.0)
            np.testing.assert_allclose(xgcrs1, xgcrs0, rtol=rtol, atol=atol)

    def test_hscrs2gcrs_kgcrs2hcrs(self):

        rtol = 1.0e-8
        atol = 1.0e-8

        for _ in range(100):
            xhcrs0 = np.random.randn(6, 1)*1000.0
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13),
                                          np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xgcrs0 = fkt.translations.khcrs2gcrs(xhcrs0, jd, 1.0, 1.0)
            xhcrs1 = fkt.translations.kgcrs2hcrs(xgcrs0, jd, 1.0, 1.0)
            np.testing.assert_allclose(xhcrs1, xhcrs0, rtol=rtol, atol=atol)
            xgcrs1 = fkt.translations.khcrs2gcrs(xhcrs1, jd, 1.0, 1.0)
            np.testing.assert_allclose(xgcrs1, xgcrs0, rtol=rtol, atol=atol)

    def test_kscrs2sors_ksors2scrs(self):

        rtol = 1.0e-11
        atol = 1.0e-11

        for _ in range(100):
            xscrs0 = np.random.randn(6, 1)
            xscrs0[0:3] = xscrs0[0:3] * 10000.0
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13),
                                          np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xsors0, dxsors0 = fkt.translations.kscrs2sors(xscrs0, jd)
            xscrs1, dxscrs1 = fkt.translations.ksors2scrs(xsors0, jd)
            np.testing.assert_allclose(xscrs1, xscrs0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxsors0[:, :, 0], np.transpose(dxscrs1[:, :, 0]), rtol=rtol, atol=atol)
            np.testing.assert_allclose(np.matmul(dxsors0[:, :, 0], dxscrs1[:, :, 0]), np.eye(6), rtol=rtol, atol=atol)
            xsors1, dxsors1 = fkt.translations.kscrs2sors(xscrs1, jd)
            np.testing.assert_allclose(xsors1, xsors0, rtol=rtol, atol=atol)
            np.testing.assert_allclose(dxsors1, dxsors0, rtol=rtol, atol=atol)

    def test_kine2rotEph_krot2ineEph(self):

        rtol = 1.0e-13
        atol = 1.0e-13

        for _ in range(100):
            xine0 = np.random.randn(6, 1)
            jd = fkt.ephemeris.juliandate(np.random.randint(2020, 2041), np.random.randint(1, 13),
                                          np.random.randint(1, 28),
                                          np.random.randint(0, 24), np.random.randint(0, 60), np.random.randint(0, 60))
            xrot0 = fkt.translations.kine2roteph(xine0, jd, 'Earth', 'Moon', 1.0, 1.0)
            xine1 = fkt.translations.krot2ineeph(xrot0, jd, 'Earth', 'Moon', 1.0, 1.0)
            np.testing.assert_allclose(xine1, xine0, rtol=rtol, atol=atol)
            xrot1 = fkt.translations.kine2roteph(xine1, jd, 'Earth', 'Moon', 1.0, 1.0)
            np.testing.assert_allclose(xrot1, xrot0, rtol=rtol, atol=atol)

    def test_kine2rot_krot2ine(self):

        rtol = 1.0e-13
        atol = 1.0e-13

        for _ in range(100):
            xine0 = np.random.randn(6, 1)
            t = [1.0]
            t0 = [0.0]
            xrot0 = fkt.translations.kine2rot(xine0, t, t0)
            xine1 = fkt.translations.krot2ine(xrot0, t, t0)
            np.testing.assert_allclose(xine1, xine0, rtol=rtol, atol=atol)
            xrot1 = fkt.translations.kine2rot(xine1, t, t0)
            np.testing.assert_allclose(xrot1, xrot0, rtol=rtol, atol=atol)


if __name__ == '__main__':
    unittest.main()
