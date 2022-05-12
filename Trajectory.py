import kiam
import numpy as np
import networkx as nx
import math
import copy

def traj2dict(tr):

    d = {}
    d['vars'] = tr.vars
    d['states'] = tr.states
    d['times'] = tr.times
    d['system'] = tr.system
    d['units_name'] = tr.units_name
    d['jds'] = tr.jds
    d['initialDate'] = tr.initialDate
    d['finalDate'] = tr.finalDate
    d['units'] = tr.units
    d['parts'] = tr.parts
    d['model'] = tr.model
    return d
def dict2traj(d):

    initial_state = d['states'][:, 0]
    initial_time = d['times'][0]
    initial_jd = d['jds'][0]
    variables = d['vars']
    system = d['system']
    units_name = d['units_name']

    tr = Trajectory(initial_state, initial_time, initial_jd, variables, system, units_name)

    tr.model = d['model']
    tr.states = d['states']
    tr.times = d['times']
    tr.jds = d['jds']
    tr.initialDate = d['initialDate']
    tr.finalDate = d['finalDate']
    tr.units = d['units']
    tr.parts = d['parts']

    return tr

class Trajectory:

    def __init__(self, initial_state, initial_time, initial_jd, variables, system, units_name):

        variables = variables.lower()
        system = system.lower()
        units_name = units_name.lower()

        if variables not in ['rv', 'rvm', 'rv_stm', 'ee', 'eem', 'ee_stm', 'oe', 'oem', 'oe_stm']:
            raise Exception('Unknown variables.')

        if variables in ['rv', 'ee', 'oe'] and len(initial_state) != 6:
            raise Exception('Wrong number of variables.')

        if variables in ['rvm', 'eem', 'oem'] and len(initial_state) != 7:
            raise Exception('Wrong number of variables.')

        if variables in ['rv_stm', 'ee_stm', 'oe_stm'] and len(initial_state) != 42:
            raise Exception('Wrong number of variables.')

        if units_name not in ['earth', 'moon', 'dim', 'earth_moon', 'sun_earth']:
            raise Exception('Unknown units_name.')

        initial_state, initial_time, initial_jd = kiam.to_float(initial_state, initial_time, initial_jd)

        self.vars = variables
        self.states = np.reshape(initial_state, (-1, 1))
        self.times = np.reshape(initial_time, (1,))
        self.system = system
        self.units_name = units_name
        self.jds = np.reshape(initial_jd, (1,))
        self.initialDate = kiam.jd2time(initial_jd)
        self.finalDate = kiam.jd2time(initial_jd)
        self.units = {}
        self.parts = []
        self.model = {}
        self.vars_graph = nx.DiGraph()
        self.systems_graph = nx.Graph()
        self.units_graph = nx.Graph()

        self.allocate_vars_graph()
        self.allocate_systems_graph()
        self.allocate_units_graph()

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
    def set_model(self, variables, model_type, primary, sources_list):

        self.model = {}

        variables = variables.lower()
        model_type = model_type.lower()
        primary = primary.lower()
        sources_list = [source.lower() for source in sources_list]

        self.model['vars'] = variables
        self.model['type'] = model_type
        self.model['primary'] = primary
        self.model['sources_list'] = sources_list
        self.model['data'] = {}
        self.model['units'] = {}

        if model_type == 'r2bp':

            if (primary != 'earth') and (primary != 'moon'):
                raise Exception('Earth or Moon as primary are required for r2bp type of model.')

            self.model['data']['grav_parameter'] = 1.0
            self.set_model_units(primary)
            if primary == 'earth':
                self.model['system'] = 'gcrs'
            elif primary == 'moon':
                self.model['system'] = 'scrs'
            else:
                raise Exception('Unknown primary')

        elif model_type == 'cr3bp_fb':

            if primary == 'earthmoon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self.set_model_units('earth_moon')
                self.model['system'] = 'gsrf_em'
            elif primary == 'sunearth':
                ku = kiam.units('Sun', 'Earth')
                self.model['data']['mass_parameter'] = ku['mu']
                self.set_model_units('sun_earth')
                self.model['system'] = 'hsrf_se'
            else:
                raise Exception('Unknown primary.')

        elif model_type == 'cr3bp_sb':

            if primary == 'earthmoon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self.set_model_units('earth_moon')
                self.model['system'] = 'ssrf_em'
            elif primary == 'sunearth':
                ku = kiam.units('Sun', 'Earth')
                self.model['data']['mass_parameter'] = ku['mu']  # Sun - (Earth + Moon)
                self.set_model_units('sun_earth')
                self.model['system'] = 'gsrf_se'
            else:
                raise Exception('Unknown primary.')

        elif model_type == 'nbp':

            self.set_model_sources()
            if variables == 'rv' and primary == 'earth':
                self.set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rv_stm' and primary == 'earth':
                self.set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rv' and primary == 'moon':
                self.set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'rv_stm' and primary == 'moon':
                self.set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'rvm' and primary == 'earth':
                self.set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rvm' and primary == 'moon':
                self.set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'ee' and primary == 'earth':
                self.set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'ee' and primary == 'moon':
                self.set_model_units('moon')
                self.model['system'] = 'scrs'
            else:
                raise Exception('Unknown model.')

    def propagate(self, tof, npoints=2):

        self.change_units(self.model['units']['name'])
        self.change_vars(self.model['vars'])
        self.change_system(self.model['system'])

        if self.vars in ['rv_stm', 'ee_stm', 'oe_stm']:
            stm = True
        else:
            stm = False

        tspan = np.linspace(self.times[-1], self.times[-1] + tof, npoints)

        if self.model['type'] == 'nbp':
            T, X = kiam.propagate_nbp(self.model['primary'], tspan, self.states[0:, -1],
                                      self.model['sources'], self.model['data'], stm, self.vars)
        elif self.model['type'] == 'r2bp':
            T, X = kiam.propagate_r2bp(tspan, self.states[0:, -1])
        elif self.model['type'] == 'cr3bp_fb':
            T, X = kiam.propagate_cr3bp(central_body='First', tspan=tspan, x0=self.states[0:, -1],
                                        mu=self.model['data']['mass_parameter'], stm=stm)
        elif self.model['type'] == 'cr3bp_sb':
            T, X = kiam.propagate_cr3bp(central_body='Secondary', tspan=tspan, x0=self.states[0:, -1],
                                        mu=self.model['data']['mass_parameter'], stm=stm)
        else:
            raise Exception('Unknown model_type.')


        if len(self.parts) == 0:
            self.parts = [0, len(T)]
        else:
            self.parts.append(self.parts[-1] + len(T))

        self.jds = np.append(self.jds[0:-1], self.jds[-1] + (T - self.times[-1]) * self.units['TimeUnit'])
        self.times = np.append(self.times[0:-1], T)
        self.states = np.append(self.states[:, 0:-1], X, axis=1)
        self.finalDate = kiam.jd2time(self.jds[-1])
    def show(self, variables):
        if self.units_name == 'dim':
            tlabel = 'Time of flight, days'
        else:
            tlabel = 'Time of flight, nondimensional'
        if self.vars in ['rv', 'rvm', 'rv_stm']:
            if variables == 'xy':
                if self.units_name == 'earth':
                    xlabel = '$x$, Earth radii'
                    ylabel = '$y$, Earth radii'
                elif self.units_name == 'moon':
                    xlabel = '$x$, Moon radii'
                    ylabel = '$y$, Moon radii'
                elif self.units_name == 'dim':
                    xlabel = '$x$, km'
                    ylabel = '$y$, km'
                elif self.units_name in ['earth_moon', 'sun_earth']:
                    xlabel = '$x$, nondimensional'
                    ylabel = '$y$, nondimensional'
                else:
                    raise Exception('Unknown units.')
                kiam.plotcol(self.states[0:2, :], xlabel=xlabel, ylabel=ylabel)
            elif variables == '3d':
                if self.units_name == 'earth':
                    xlabel = '$x$, Earth radii'
                    ylabel = '$y$, Earth radii'
                    zlabel = '$z$, Earth radii'
                elif self.units_name == 'moon':
                    xlabel = '$x$, Moon radii'
                    ylabel = '$y$, Moon radii'
                    zlabel = '$z$, Moon radii'
                elif self.units_name == 'dim':
                    xlabel = '$x$, km'
                    ylabel = '$y$, km'
                    zlabel = '$y$, km'
                elif self.units_name in ['earth_moon', 'sun_earth']:
                    xlabel = '$x$, nondimensional'
                    ylabel = '$y$, nondimensional'
                    zlabel = '$z$, nondimensional'
                else:
                    raise Exception('Unknown units.')
                kiam.plotcol(self.states[0:3, :], xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
        elif self.vars in ['oe', 'oem', 'oe_stm']:
            if variables == 'a':
                if self.units_name == 'dim':
                    ylabel = 'Semi-major axis, km'
                elif self.units_name == 'earth':
                    ylabel = 'Semi-major axis, Earth''s radii'
                else:
                    ylabel = ''
                kiam.plot(self.times, self.states[0, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'e':
                ylabel = 'Eccentricity'
                kiam.plot(self.times, self.states[1, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'inc':
                ylabel = 'Inclination, degrees'
                kiam.plot(self.times, self.states[2, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'Om':
                ylabel = 'Longitude of the ascending node, degrees'
                kiam.plot(self.times, self.states[3, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'w':
                ylabel = 'Argument of pericenter, degrees'
                kiam.plot(self.times, self.states[4, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'th':
                ylabel = 'True anomaly, degrees'
                kiam.plot(self.times, self.states[5, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
        elif self.vars in ['ee', 'eem', 'ee_stm']:
            if variables == 'h':
                if self.units_name == 'dim':
                    ylabel = r'$h$, (km/s)^{-1}'
                elif self.units_name == 'earth':
                    ylabel = r'$h$, nondimensional'
                else:
                    ylabel = ''
                kiam.plot(self.times, self.states[0, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ex':
                ylabel = '$e_x$'
                kiam.plot(self.times, self.states[1, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ey':
                ylabel = '$e_y$'
                kiam.plot(self.times, self.states[2, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ix':
                ylabel = '$i_x$'
                kiam.plot(self.times, self.states[3, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'iy':
                ylabel = '$i_y$'
                kiam.plot(self.times, self.states[4, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'L':
                ylabel = 'True longitude, degrees'
                kiam.plot(self.times, self.states[5, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
    def copy(self):
        return copy.deepcopy(self)
    def clear(self):
        self.states = np.reshape(self.states[:, 0], (6, 1))
        self.times = np.reshape(self.times[0], (1,))
        self.jds = np.reshape(self.jds[0], (1,))
        self.parts = []
        self.finalDate = self.initialDate

    def change_vars(self, new_vars):
        if self.vars == new_vars:
            return
        if new_vars not in self.vars_graph.nodes:
            raise Exception('Unknown new_vars.')
        p = nx.shortest_path(self.vars_graph, self.vars, new_vars)
        for i in range(len(p)-1):
            self.vars_transform(p[i], p[i+1])
    def change_system(self, new_system):
        if self.system == new_system:
            return
        if new_system not in self.systems_graph.nodes:
            raise Exception('Unknown new_system.')
        p = nx.shortest_path(self.systems_graph, self.system, new_system)
        for i in range(len(p)-1):
            self.system_transform(p[i], p[i+1])
    def change_units(self, new_units):
        if self.units_name == new_units:
            return
        if new_units not in self.units_graph.nodes:
            raise Exception('Unknown new_units.')
        p = nx.shortest_path(self.units_graph, self.units_name, new_units)
        for i in range(len(p)-1):
            self.units_transform(p[i], p[i+1])

    # Variables transformations.
    def allocate_vars_graph(self):

        self.vars_graph.add_edge('rv', 'ee')
        self.vars_graph.add_edge('ee', 'rv')

        self.vars_graph.add_edge('rv', 'oe')
        self.vars_graph.add_edge('oe', 'rv')

        self.vars_graph.add_edge('rvm', 'eem')
        self.vars_graph.add_edge('eem', 'rvm')

        self.vars_graph.add_edge('rvm', 'oem')
        self.vars_graph.add_edge('oem', 'rvm')

        self.vars_graph.add_edge('rv_stm', 'rv')
        self.vars_graph.add_edge('oe_stm', 'oe')
        self.vars_graph.add_edge('ee_stm', 'ee')

        self.vars_graph.add_edge('rvm', 'rv')
        self.vars_graph.add_edge('oem', 'oe')
        self.vars_graph.add_edge('eem', 'ee')

        self.vars_graph.add_edge('rv_stm', 'oe_stm')
        self.vars_graph.add_edge('oe_stm', 'rv_stm')

        self.vars_graph.add_edge('rv_stm', 'ee_stm')
        self.vars_graph.add_edge('ee_stm', 'rv_stm')
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
                phi_rv = np.reshape(self.states[6:42, i], (6, 6)).T
                phi_oe = kiam.dotAinvB(np.matmul(doe, phi_rv), doe0)
                self.states[6:42, i] = np.reshape(phi_oe.T, (36,))
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
                phi_oe = np.reshape(self.states[6:42, i], (6, 6)).T
                phi_rv = kiam.dotAinvB(np.matmul(drv, phi_oe), drv0)
                self.states[6:42, i] = np.reshape(phi_rv.T, (36,))
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
                phi_rv = np.reshape(self.states[6:42, i], (6, 6)).T
                phi_ee = kiam.dotAinvB(np.matmul(dee, phi_rv), dee0)
                self.states[6:42, i] = np.reshape(phi_ee.T, (36,))
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
                phi_ee = np.reshape(self.states[6:42, i], (6, 6)).T
                phi_rv = kiam.dotAinvB(np.matmul(drv, phi_ee), drv0)
                self.states[6:42, i] = np.reshape(phi_rv.T, (36,))
            self.vars = 'rv_stm'
        else:
            raise Exception('Unknown variable transformaton.')

    # System transformations.
    def allocate_systems_graph(self):
        self.systems_graph.add_edge('itrs', 'gcrs')
        self.systems_graph.add_edge('gcrs', 'gsrf_em')
        self.systems_graph.add_edge('gcrs', 'gsrf_se')
        self.systems_graph.add_edge('scrs', 'ssrf_em')
        self.systems_graph.add_edge('gcrs', 'scrs')
        self.systems_graph.add_edge('scrs', 'mer')
        self.systems_graph.add_edge('gcrs', 'hcrs')
        self.systems_graph.add_edge('hcrs', 'hsrf_se')
        self.systems_graph.add_edge('scrs', 'sors')
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
            if self.vars not in ['rv', 'rvm', 'rv_stm']:
                raise Exception('Vars should be rv or rvm')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            self.states[0:6, :] = kiam.scrs2mer(self.states[0:6, :], self.jds)
            self.system = 'mer'
        elif system1 == 'mer' and system2 == 'scrs':
            if self.vars not in ['rv', 'rvm', 'rv_stm']:
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
            if self.vars not in ['rv', 'rvm', 'rv_stm']:
                raise Exception('Vars should be rv, rvm or rv_stm.')
            elif self.system != 'scrs':
                raise Exception('System should be scrs.')
            if self.vars != 'rv_stm':
                self.states[0:6, :] = kiam.scrs2sors(self.states[0:6, :], self.jds, False)
            else:
                xsors, dxsors = kiam.scrs2sors(self.states[0:6, :], self.jds, True)
                self.states[0:6, :] = xsors
                for i in range(dxsors.shape[2]):
                    phi_scrs = np.reshape(self.states[6:42, i], (6, 6)).T
                    phi_sors = kiam.dotAinvB(np.matmul(dxsors[:, :, i], phi_scrs), dxsors[:, :, 0])
                    self.states[6:42, i] = np.reshape(phi_sors.T, (36,))
            self.system = 'sors'
        elif system1 == 'sors' and system2 == 'scrs':
            if self.vars not in ['rv', 'rvm', 'rv_stm']:
                raise Exception('Vars should be rv, rvm or rv_stm.')
            elif self.system != 'sors':
                raise Exception('System should be sors.')
            if self.vars != 'rv_stm':
                self.states[0:6, :] = kiam.sors2scrs(self.states[0:6, :], self.jds, False)
            else:
                xscrs, dxscrs = kiam.sors2scrs(self.states[0:6, :], self.jds, True)
                self.states[0:6, :] = xscrs
                for i in range(dxscrs.shape[2]):
                    phi_sors = np.reshape(self.states[6:42, i], (6, 6)).T
                    phi_scrs = kiam.dotAinvB(np.matmul(dxscrs[:, :, i], phi_sors), dxscrs[:, :, 0])
                    self.states[6:42, i] = np.reshape(phi_scrs.T, (36,))
            self.system = 'scrs'

    # Units transformations and settings.
    def allocate_units_graph(self):
        self.units_graph.add_edge('dim', 'earth')
        self.units_graph.add_edge('dim', 'moon')
        self.units_graph.add_edge('dim', 'earth_moon')
        self.units_graph.add_edge('dim', 'sun_earth')
    def units_transform(self, units1, units2):
        if units1 == 'dim' and units2 == 'earth':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self.set_earth_units()
                self.undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self.set_earth_units()
                self.undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self.set_earth_units()
                self.undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'earth'
        elif units1 == 'dim' and units2 == 'moon':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self.set_moon_units()
                self.undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self.set_moon_units()
                self.undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self.set_moon_units()
                self.undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'moon'
        elif units1 == 'dim' and units2 == 'earth_moon':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self.set_earth_moon_units()
                self.undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self.set_earth_moon_units()
                self.undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self.set_earth_moon_units()
                self.undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'earth_moon'
        elif units1 == 'dim' and units2 == 'sun_earth':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self.set_sun_earth_units()
                self.undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self.set_sun_earth_units()
                self.undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self.set_sun_earth_units()
                self.undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'sun_earth'
        elif units2 == 'dim':
            if units1 != self.units_name:
                raise Exception(f'Not {units1} units.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self.dim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self.dim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self.dim_oe()
            else:
                raise Exception('Unknown vars.')
            self.set_dim_units()
            self.units_name = 'dim'
    def undim_rv(self):
        self.times = self.times/self.units['TimeUnit']
        self.states[0:3, :] = self.states[0:3, :]/self.units['DistUnit']
        self.states[3:6, :] = self.states[3:6, :]/self.units['VelUnit']
    def dim_rv(self):
        self.times = self.times*self.units['TimeUnit']
        self.states[0:3, :] = self.states[0:3, :]*self.units['DistUnit']
        self.states[3:6, :] = self.states[3:6, :]*self.units['VelUnit']
    def undim_ee(self):
        self.times = self.times/self.units['TimeUnit']
        self.states[0, :] = self.states[0, :]*self.units['VelUnit']
    def dim_ee(self):
        self.times = self.times*self.units['TimeUnit']
        self.states[0, :] = self.states[0, :]/self.units['VelUnit']
    def undim_oe(self):
        self.times = self.times/self.units['TimeUnit']
        self.states[0, :] = self.states[0, :]/self.units['DistUnit']
    def dim_oe(self):
        self.times = self.times*self.units['TimeUnit']
        self.states[0, :] = self.states[0, :]*self.units['DistUnit']
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

    # Auxilary model routines.
    def set_model_units(self, units_name):

        if units_name == 'earth':

            self.model['units']['name'] = 'earth'
            units = kiam.units('Earth')
            self.model['units']['mu'] = units['GM']

        elif units_name == 'moon':

            self.model['units']['name'] = 'moon'
            units = kiam.units('Moon')
            self.model['units']['mu'] = units['GM']

        elif units_name == 'earth_moon':

            self.model['units']['name'] = 'earth_moon'
            units = kiam.units('Earth', 'Moon')

        elif units_name == 'sun_earth':

            self.model['units']['name'] = 'sun_earth'
            units = kiam.units('Sun', 'Earth')

        else:

            raise Exception('Unknown units_name.')

        self.model['units']['DistUnit'] = units['DistUnit']  # km
        self.model['units']['VelUnit'] = units['VelUnit']  # km/s
        self.model['units']['TimeUnit'] = units['TimeUnit']  # days
        self.model['units']['AccUnit'] = units['AccUnit']  # m/s^2

        _, star, planet, moon, _ = kiam.astro_const()

        self.model['units']['RSun'] = star['Sun']['MeanRadius'] / units['DistUnit']  # 695700 km
        self.model['units']['REarth'] = planet['Earth']['EquatorRadius'] / units['DistUnit']  # 6378.1366 km
        self.model['units']['RMoon'] = moon['Moon']['MeanRadius'] / units['DistUnit']  # 1737.4 km
    def set_model_sources(self):

        self.model['sources'] = kiam.prepare_sources_dict()

        for source in self.model['sources_list']:
            self.model['sources'][source.lower()] = True
