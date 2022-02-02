import kiam
import Model
import numpy as np

class Trajectory:

    def __init__(self, initial_state, initial_time, initial_jd, variables, system, units_name):

        self.vars = variables

        if variables not in ['rv', 'rvm', 'rv_stm', 'ee', 'eem', 'oe', 'oem']:
            raise Exception('Unknown variables.')

        if variables in ['rv', 'rv_stm', 'ee', 'oe'] and len(initial_state) != 6:
            raise Exception('Wrong number of variables.')

        if variables in ['rvm', 'eem', 'oem'] and len(initial_state) != 7:
            raise Exception('Wrong number of variables.')

        if units_name not in ['earth', 'moon', 'dim', 'earth_moon', 'sun_earth']:
            raise Exception('Unknown units_name.')

        self.states = np.reshape(initial_state, (-1, 1))
        self.times = np.reshape(initial_time, (1,))
        self.system = system
        self.units_name = units_name
        self.jds = np.reshape(initial_jd, (1,))
        self.initialDate = kiam.jd2time(initial_jd)
        self.finalDate = kiam.jd2time(initial_jd)
        self.units = {}
        self.parts = []
        self.model = []

        if units_name == 'earth':
            self.set_earth_units()
        elif units_name == 'dim':
            self.set_dim_units()
        elif units_name == 'earth_moon':
            self.set_earth_moon_units()
        elif units_name == 'sun_earth':
            self.set_sun_earth_units()
        elif units_name == 'moon':
            self.set_moon_units()

    def set_model(self, variables, model_type, primary, sources_cell):
        self.model = Model.Model(variables, model_type, primary, sources_cell)

    # Variables transformations.
    def vars_transform(self, vars1, vars2):
        if vars1 == 'rv' and vars2 == 'ee':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv2ee suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm.')
            for i in range(self.states.shape[1]):
                self.states[0:6, i] = kiam.rv2ee(self.states[0:6, i], 1.0)
            self.vars = 'ee'
        elif vars1 == 'ee' and vars2 == 'rv':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: ee2rv suggests earth or moon, mu = 1.0.')
            elif self.vars != 'ee' and self.vars != 'eem':
                raise Exception('Vars should be ee or eem.')
            for i in range(self.states.shape[1]):
                self.states[0:6, i] = kiam.ee2rv(self.states[0:6, i], 1.0)
            self.vars = 'rv'
        elif vars1 == 'rv' and vars2 == 'oe':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv2oe suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm.')
            for i in range(self.states.shape[1]):
                self.states[0:6, i] = kiam.rv2oe(self.states[0:6, i], 1.0)
            self.vars = 'oe'
        elif vars1 == 'oe' and vars2 == 'rv':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: oe2rv suggests earth or moon, mu = 1.0.')
            elif self.vars != 'oe' and self.vars != 'oem':
                raise Exception('Vars should be oe or oem.')
            for i in range(self.states.shape[1]):
                self.states[0:6, i] = kiam.oe2rv(self.states[0:6, i], 1.0)
            self.vars = 'rv'
        elif vars1 == 'rvm' and vars2 == 'eem':
            self.vars_transform('rv', 'ee')
            self.vars = 'eem'
        elif vars1 == 'eem' and vars2 == 'rvm':
            self.vars_transform('ee', 'rv')
            self.vars = 'rvm'
        elif vars1 == 'rvm' and vars2 == 'oem':
            self.vars_transform('rv', 'oe')
            self.vars = 'oem'
        elif vars1 == 'oem' and vars2 == 'rvm':
            self.vars_transform('oe', 'rv')
            self.vars = 'rvm'
        elif vars1 == 'rv_stm' and vars2 == 'rv':
            np.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'rv'
        elif vars1 == 'oe_stm' and vars2 == 'oe':
            np.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'oe'
        elif vars1 == 'ee_stm' and vars2 == 'ee':
            np.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'ee'
        elif vars1 == 'rvm' and vars2 == 'rv':
            np.delete(self.states, 6, 0)
            self.vars = 'rv'
        elif vars1 == 'oem' and vars2 == 'oe':
            np.delete(self.states, 6, 0)
            self.vars = 'oe'
        elif vars1 == 'eem' and vars2 == 'ee':
            np.delete(self.states, 6, 0)
            self.vars = 'ee'
        elif vars1 == 'rv_stm' and vars2 == 'oe_stm':
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv_stm2oe_stm suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv_stm':
                raise Exception('Vars should be rv_stm.')
            _, doe0 = kiam.rv2oe(self.states[0:6, 0], 1.0, True)
            for i in range(self.states.shape[1]):
                oe, doe = kiam.rv2oe(self.states[0:6, i], 1.0, True)
                self.states[0:6, i] = oe
                phi_rv = np.reshape(self.states[6:42, i], (6, 6))
                phi_oe = kiam.dotAinvB(np.matmul(doe, phi_rv), doe0)
                self.states[6:42, i] = np.reshape(phi_oe, (36, 1))
            self.vars = 'oe_stm'
        elif vars1 == 'oe_stm' and vars2 == 'rv_stm':
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: oe_stm2rv_stm suggests earth or moon, mu = 1.0.')
            elif self.vars != 'oe_stm':
                raise Exception('Vars should be oe_stm.')
            _, drv0 = kiam.oe2rv(self.states[0:6, 0], 1.0, True)
            for i in range(self.states.shape[1]):
                rv, drv = kiam.oe2rv(self.states[0:6, i], 1.0, True)
                self.states[0:6, i] = rv
                phi_oe = np.reshape(self.states[6:42, i], (6, 6))
                phi_rv = kiam.dotAinvB(np.matmul(drv, phi_oe), drv0)
                self.states[6:42, i] = np.reshape(phi_rv, (36, 1))
            self.vars = 'rv_stm'
        elif vars1 == 'rv_stm' and vars2 == 'ee_stm':
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv_stm2ee_stm suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv_stm':
                raise Exception('Vars should be rv_stm.')
            _, dee0 = kiam.rv2ee(self.states[0:6, 0], 1.0, True)
            for i in range(self.states.shape[1]):
                ee, dee = kiam.rv2ee(self.states[0:6, i], 1.0, True)
                self.states[0:6, i] = ee
                phi_rv = np.reshape(self.states[6:42, i], (6, 6))
                phi_ee = kiam.dotAinvB(np.matmul(dee, phi_rv), dee0)
                self.states[6:42, i] = np.reshape(phi_ee, (36, 1))
            self.vars = 'ee_stm'
        elif vars1 == 'ee_stm' and vars2 == 'rv_stm':
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: ee_stm2rv_stm suggests earth or moon, mu = 1.0.')
            elif self.vars != 'ee_stm':
                raise Exception('Vars should be ee_stm.')
            _, drv0 = kiam.ee2rv(self.states[0:6, 0], 1.0, True)
            for i in range(self.states.shape[1]):
                rv, drv = kiam.ee2rv(self.states[0:6, i], 1.0, True)
                self.states[0:6, i] = rv
                phi_ee = np.reshape(self.states[6:42, i], (6, 6))
                phi_rv = kiam.dotAinvB(np.matmul(drv, phi_ee), drv0)
                self.states[6:42, i] = np.reshape(phi_rv, (36, 1))
            self.vars = 'rv_stm'
        else:
            raise Exception('Unknown variable transformaton.')

    # System transformations.
    def system_transform(self, system1, system2):
        if system1 == 'itrs' and system2 == 'gcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'itrs':
                raise Exception('System should be itrs.')
            for i in range(self.states.shape[1]):
                self.states[0:3, i] = kiam.itrs2gcrs(self.states[0:3, i], self.jds[i])
                self.states[3:6, i] = kiam.itrs2gcrs(self.states[3:6, i], self.jds[i])
            self.system = 'gcrs'
        elif system1 == 'gcrs' and system2 == 'itrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gcrs':
                raise Exception('System should be gcrs.')
            for i in range(self.states.shape[1]):
                self.states[0:3, i] = kiam.gcrs2itrs(self.states[0:3, i], self.jds[i])
                self.states[3:6, i] = kiam.gcrs2itrs(self.states[3:6, i], self.jds[i])
            self.system = 'itrs'
        elif system1 == 'gcrs' and system2 == 'gsrf_em':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gcrs':
                raise Exception('System should be gcrs.')
            self.states[0:6, :] = kiam.ine2rotEph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gsrf_em'
        elif system1 == 'gsrf_em' and system2 == 'gcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gsrf_em':
                raise Exception('System should be gsrf_em.')
            self.states[0:6, :] = kiam.rot2ineEph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gcrs'
        elif system1 == 'gcrs' and system2 == 'gsrf_se':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gcrs':
                raise Exception('System should be gcrs.')
            self.states[0:6, :] = kiam.ine2rotEph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gsrf_se'
        elif system1 == 'gsrf_se' and system2 == 'gcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gsrf_se':
                raise Exception('System should be gsrf_se.')
            self.states[0:6, :] = kiam.rot2ineEph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gcrs'
        elif system1 == 'scrs' and system2 == 'ssrf_em':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            self.states[0:6, :] = kiam.ine2rotEph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'ssrf_em'
        elif system1 == 'ssrf_em' and system2 == 'scrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'ssrf_em':
                raise Exception('System should be ssrf_em.')
            self.states[0:6, :] = kiam.rot2ineEph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'scrs'
        elif system1 == 'gcrs' and system2 == 'scrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gcrs':
                raise Exception('System should be gcrs.')
            self.states[0:6, :] = kiam.gcrs2scrs(self.states[0:6, :], self.jds,
                                                 self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'scrs'
        elif system1 == 'scrs' and system2 == 'gcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            self.states[0:6, :] = kiam.scrs2gcrs(self.states[0:6, :], self.jds,
                                                 self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gcrs'
        elif system1 == 'scrs' and system2 == 'mer':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            self.states[0:6, :] = kiam.scrs2mer(self.states[0:6, :], self.jds)
            self.system = 'mer'
        elif system1 == 'mer' and system2 == 'scrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'mer':
                raise Exception('System should be mer.')
            self.states[0:6, :] = kiam.mer2scrs(self.states[0:6, :], self.jds)
            self.system = 'scrs'
        elif system1 == 'gcrs' and system2 == 'hcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'gcrs':
                raise Exception('System should be gcrs.')
            self.states[0:6, :] = kiam.gcrs2hcrs(self.states[0:6, :], self.jds,
                                                 self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'hcrs'
        elif system1 == 'hcrs' and system2 == 'gcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'hcrs':
                raise Exception('System should be hcrs.')
            self.states[0:6, :] = kiam.hcrs2gcrs(self.states[0:6, :], self.jds,
                                                 self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'gcrs'
        elif system1 == 'hcrs' and system2 == 'hsrf_se':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'hcrs':
                raise Exception('System should be hcrs.')
            self.states[0:6, :] = kiam.ine2rotEph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'hsrf_se'
        elif system1 == 'hsrf_se' and system2 == 'hcrs':
            if self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'hsrf_se':
                raise Exception('System should be hsrf_se.')
            self.states[0:6, :] = kiam.rot2ineEph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                  self.units['DistUnit'], self.units['VelUnit'])
            self.system = 'hcrs'
        elif system1 == 'scrs' and system2 == 'sors':
            if self.vars != 'rv' and self.vars != 'rvm' and self.vars != 'rv_stm':
                raise Exception('Vars should be rv, rvm or rv_stm.')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            if self.vars != 'rv_stm':
                self.states[0:6, :] = kiam.scrs2sors(self.states[0:6, :], self.jds, False)
            else:
                xsors, dxsors = kiam.scrs2sors(self.states[0:6, :], self.jds, True)
                self.states[0:6, :] = xsors
                for i in range(dxsors.shape[2]):
                    phi_scrs = np.reshape(self.states[6:42, i], (6, 6))
                    phi_sors = kiam.dotAinvB(np.matmul(dxsors[:, :, i], phi_scrs), dxsors[:, :, 0])
                    self.states[6:42, i] = np.reshape(phi_sors, (36, 1))
            self.system = 'sors'
        elif system1 == 'sors' and system2 == 'scrs':
            if self.vars != 'rv' and self.vars != 'rvm' and self.vars != 'rv_stm':
                raise Exception('Vars should be rv, rvm or rv_stm.')
            elif self.system != 'sors':
                raise Exception('System should be sors.')
            if self.vars != 'rv_stm':
                self.states[0:6, :] = kiam.sors2scrs(self.states[0:6, :], self.jds, False)
            else:
                xscrs, dxscrs = kiam.sors2scrs(self.states[0:6, :], self.jds, True)
                self.states[0:6, :] = xscrs
                for i in range(dxscrs.shape[2]):
                    phi_sors = np.reshape(self.states[6:42, i], (6, 6))
                    phi_scrs = kiam.dotAinvB(np.matmul(dxscrs[:, :, i], phi_sors), dxscrs[:, :, 0])
                    self.states[6:42, i] = np.reshape(phi_scrs, (36, 1))
            self.system = 'scrs'


    # Units transformations and settings.
    def set_earth_units(self):
        ku = kiam.units('Earth')
        self.units['mu'] = 1.0
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']    # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']    # m/s^2
    def set_moon_units(self):
        ku = kiam.units('Moon')
        self.units['mu'] = 1.0
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']    # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']    # m/s^2
    def set_earth_moon_units(self):
        ku = kiam.units('Earth', 'Moon')
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']    # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']    # m/s^2
    def set_sun_earth_units(self):
        ku = kiam.units('Sun', 'Earth')
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']    # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']    # m/s^2
    def set_dim_units(self):
        ku = kiam.units('Earth')
        self.units['mu'] = ku['GM']    # km^3/s^2
        self.units['DistUnit'] = 1.0   # km
        self.units['VelUnit'] = 1.0    # km/s
        self.units['TimeUnit'] = 1.0   # days
        self.units['AccUnit'] = 1.0    # m/s^2
