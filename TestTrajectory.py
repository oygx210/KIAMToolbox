import kiam
import Trajectory
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
        tr.model['data']['jd_zero'] = jd0
        tr.propagate(period, 1000)
        np.testing.assert_allclose(tr.states[:, 0], tr.states[:, -1], rtol=rtol, atol=atol)
        tr.show('3d')

    def test_cr3bp(self):

        t0 = 0.0
        x0 = np.array([0.5, 0.0, 0.0, 0.0, 0.5, 0.0])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gsrf_em', 'earth_moon')
        tr.set_model('rv', 'cr3bp_sb', 'earthmoon', [])
        tr.model['data']['jd_zero'] = jd0
        tr.propagate(2*pi, 1000)
        tr.show('xy')

    def test_nbp(self):

        t0 = 0
        x0 = np.array([1, 0, 0, 0, 1, 0])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
        tr.set_model('rv', 'nbp', 'earth', ['Sun', 'Moon'])
        tr.model['data']['jd_zero'] = jd0
        tr.model['data']['area'] = 2.0
        tr.model['data']['mass'] = 100.0
        tr.model['data']['order'] = 50
        tr.propagate(2 * pi, 1000)
        tr.show('xy')

    def test_nbp_stm(self):

        t0 = 0.0
        x0 = np.append(np.array([1, 0, 0, 0, 1, 0]), kiam.eyevec(6))
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv_stm', 'gcrs', 'earth')
        tr.set_model('rv_stm', 'nbp', 'earth', ['Sun', 'Moon'])
        tr.model['data']['jd_zero'] = jd0
        tr.model['data']['area'] = 2.0
        tr.model['data']['mass'] = 100.0
        tr.model['data']['order'] = 50
        tr.propagate(2 * pi, 1000)
        tr.show('xy')

    def test_to_from_dict(self):

        t0 = 0
        x0 = np.array([1, 0, 0, 0, 1, 0])
        jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
        tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
        tr.set_model('rv', 'nbp', 'earth', ['Sun', 'Moon'])
        tr.model['data']['jd_zero'] = jd0
        tr.model['data']['area'] = 2.0
        tr.model['data']['mass'] = 100.0
        tr.model['data']['order'] = 50
        tr.propagate(2 * pi, 1000)

        d = Trajectory.traj2dict(tr)
        tr2 = Trajectory.dict2traj(d)
        d2 = Trajectory.traj2dict(tr2)

        print(d == d2)

        print(tr.vars == tr2.vars)
        print((tr.states == tr2.states).all())
        print((tr.times == tr2.times).all())
        print(tr.system == tr2.system)
        print(tr.units_name == tr2.units_name)
        print((tr.jds == tr2.jds).all())
        print(tr.initialDate == tr2.initialDate)
        print(tr.finalDate == tr2.finalDate)
        print(tr.units == tr2.units)
        print(tr.parts == tr2.parts)
        print(tr.model == tr2.model)


if __name__ == '__main__':
    unittest.main()
