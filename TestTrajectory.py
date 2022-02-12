import kiam
import Trajectory
import Model
import unittest
import numpy as np
from numpy import pi

class TestTrajectory(unittest.TestCase):

    def test_r2bp(self):

        rtol = 1.0e-8
        atol = 1.0e-8

        t0 = 0.0
        x0 = np.array([1.5, 0.0, 0.0, 0.0, 1.0 / np.sqrt(1.5), 0.1])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        oe0 = kiam.rv2oe(x0, 1.0)
        period = 2*pi*np.sqrt(oe0[0]**3)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
        tr.set_model('rv', 'r2bp', 'earth', [])
        tr.model.data['jd_zero'] = jd0
        tr.propagate(period, 1000)
        np.testing.assert_allclose(tr.states[:, 0], tr.states[:, -1], rtol=rtol, atol=atol)
        tr.show('3d')

    def test_cr3bp(self):

        t0 = 0.0
        x0 = np.array([0.5, 0.0, 0.0, 0.0, 0.5, 0.0])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gsrf_em', 'earth_moon')
        tr.set_model('rv', 'cr3bp_sb', 'earthmoon', [])
        tr.model.data['jd_zero'] = jd0
        tr.propagate(2*pi, 1000)
        tr.show('xy')

    def test_nbp(self):

        t0 = 0.0
        x0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0 / np.sqrt(1.0), 0.0])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
        tr.set_model('rv', 'nbp', 'earth', [])
        tr.model.data['jd_zero'] = jd0
        tr.propagate(2 * pi, 1000)
        tr.show('xy')


if __name__ == '__main__':
    unittest.main()