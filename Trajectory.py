"""
This Python module is a part of the KIAM Astrodynamics Toolbox developed in
Keldysh Institute of Applied Mathematics (KIAM), Moscow, Russia.

The module provides a useful and convenient a class Trajectory for the design of
space trajectories. With this class, the users can propagate trajectories
within the selected model of motion, perform change of variables,
coordinate systems, and units systems, and visualize the trajectories and
their characteristics.

The toolbox is licensed under the MIT License.

For more information see GitHub page of the project:
https://github.com/shmaxg/KIAMToolbox
"""

import kiam
import numpy
import networkx as nx
import math
import copy

variables_all = {'rv', 'rvm', 'rv_stm',
                 'ee', 'eem', 'ee_stm',
                 'oe', 'oem', 'oe_stm'}

systems_ephemeris = {'itrs', 'gcrs', 'gsrf_em', 'gsrf_se',  # Earth-cent. systems
                     'scrs', 'ssrf_em', 'mer', 'sors',  # Moon-cent. systems
                     'hcrs', 'hsrf_se'}  # Sun-cent. systems

systems_cr3bp = {'ine_fb', 'ine_sb', 'ine_cm', 'rot_fb', 'rot_sb', 'rot_cm'}

systems_br4bp = {'ine_fb', 'ine_sb', 'ine_cm', 'rot_fb', 'rot_sb', 'rot_cm'}

systems_hill = {'rot_sb', 'ine_sb'}

systems_all = {*systems_cr3bp, *systems_br4bp, *systems_hill, *systems_ephemeris}

units_all = {'earth', 'moon', 'dim', 'earth_moon', 'sun_earth'}

models_nbp = {'nbp'}

models_cr3bp = {'cr3bp_fb', 'cr3bp_sb'}

models_br4bp = {'br4bp_fb', 'br4bp_sb'}

models_hill = {'hill'}

models_all = {*models_cr3bp, *models_br4bp, *models_hill, *models_nbp}


class Trajectory:

    def __init__(self, initial_state: numpy.ndarray, initial_time: float, initial_jd: float,
                 variables: str, system: str, units_name: str) -> None:
        """
        Initialize the Trajectory object with given initial conditions and properties.

        Parameters:
        -----------
        `initial_state` : numpy.ndarray, shape (6,), (7,), (42,)

        Initial state of the spacecraft.

        It can be a

        1. position-velocity vector ('rv'), shape (6,)

        2. position-velocity-mass vector ('rvm'), shape (7,)

        3. position-velocity-stm vector('rv_stm'), shape (42,)

        4. equinoctial orbital elements vector ('ee'), shape (6,)

        5. equinoctial orbital elements extended by mass vector ('eem'), shape (7,)

        6. equinoctial orbital elements extended by stm vector ('ee_stm'), shape (42,)

        7. classical orbital elements vector ('oe'), shape (6,)

        8. classical orbital elements extended by mass vector ('oem'), shape (7,)

        9. classical orbital elements extended by stm vector ('oe_stm'), shape (42,)

        stm means state-transition matrix

        The classical orbital elements:

        a (semi-major axis),

        e (eccentricity),

        i (inclination),

        Omega (right ascension of the ascending node),

        omega (argument of pericenter),

        theta (true anomaly)

        The equinoctial orbital elements:

        h, ex, ey, ix, iy, L

        h = sqrt(p/mu),

        ex = e*cos(Omega+omega),

        ey = e*sin(Omega+omega),

        ix = tan(i/2)*cos(Omega),

        iy = tan(i/2)*sin(Omega),

        L = theta + omega + Omega,

        where

        mu - gravitational parameter,

        p - semi-latus rectum (focal parameter)

        `initial_time` : float

        Initial time.

        `initial_jd` : float

        Julian date corresponding to initial_time.

        `variables` : str

        Variables in terms of which initial_state is given.

        It can be 'rv', 'rvm', 'rv_stm', 'ee', 'eem', 'ee_stm', 'oe', 'oem', 'oe_stm'.

        `system` : str

        The coordinate system in which initial_state is given.

        It can be

        'scrs', 'mer', 'sors', 'ssrm_em' (the Moon-centered systems)

        'gcrs', 'itrs', 'gsrf_em', 'gsrf_se' (the Earth-centered systems)

        'hcrs', 'hsrf_se' (the Sun-centered systems)

        'ine_fb', 'rot_fb' (first primary-centered systems in CR3BP)

        'ine_sb', 'rot_sb' (second primary-centered systems in CR3BP)

        'ine_cm', 'rot_cm' (baricenter-centered systems in CR3BP)

        `units_name` : str

        The name of the units in which initial_state is given.

        It can be 'earth', 'moon', 'dim', 'earth_moon', 'sun_earth'.

        Properties:
        -----------
        `finalDate` : datetime.datetime

        The date and time of the last phase state in trajectory.

        `initialDate` : datetime.datetime

        The date and time of the first phase state in trajectory.

        `jds` : numpy.ndarray, shape(n,)

        Julian dates corresponding to phase vectors in states property.

        `model` : dict

        The model used for trajectory propagation.

        `parts` : list[int]

        The numbers of nodes between which the trajectory was sequentially propagated.

        `states` : numpy.array, shape(m,n)

        Phase vectors along the trajectory.

        `system` : str

        The current coordinate system.

        `system_graph` : networkx.classes.graph.Graph

        The graph of transormations between coordinate systems.

        `times` : numpy.array(n,)

        Times along the trajectory.

        `units` : dict

        The dictionary of the current units.

        `units_graph` : networkx.classes.graph.Graph

        The graph of transormations between units.

        `units_name` : str

        The name of the current units.

        `vars` : str

        The current variables.

        `vars_graph` : networkx.classes.digraph.DiGraph

        The graph of transormations between variables.

        Examples:
        ---------

        Examples can be found in kiam_examples.py.

        """

        variables = variables.lower()
        system = system.lower()
        units_name = units_name.lower()

        if variables not in variables_all:
            raise Exception('Unknown variables.')

        if system not in systems_all:
            raise Exception('Unknown system.')

        if units_name not in units_all:
            raise Exception('Unknown units_name.')

        if variables in ['rv', 'ee', 'oe'] and len(initial_state) != 6:
            raise Exception('Wrong number of variables.')

        if variables in ['rvm', 'eem', 'oem'] and len(initial_state) != 7:
            raise Exception('Wrong number of variables.')

        if variables in ['rv_stm', 'ee_stm', 'oe_stm'] and len(initial_state) != 42:
            raise Exception('Wrong number of variables.')

        initial_state, initial_time, initial_jd = kiam.to_float(initial_state, initial_time, initial_jd)

        self.vars = variables
        self.states = numpy.reshape(initial_state, (-1, 1))
        self.times = numpy.reshape(initial_time, (1,))
        self.system = system
        self.units_name = units_name
        self.jds = numpy.reshape(initial_jd, (1,))
        self.initialDate = kiam.jd2time(initial_jd)
        self.finalDate = kiam.jd2time(initial_jd)
        self.units = {}
        self.parts = []
        self.model = {}
        self.vars_graph = nx.DiGraph()
        self.systems_graph = nx.Graph()
        self.units_graph = nx.Graph()

        self._allocate_vars_graph()
        self._allocate_systems_graph()
        self._allocate_units_graph()

        if units_name == 'earth':
            self._set_earth_units()
        elif units_name == 'dim':
            self._set_dim_units()
        elif units_name == 'earth_moon':
            self._set_earth_moon_units()
        elif units_name == 'sun_earth':
            self._set_sun_earth_units()
        elif units_name == 'moon':
            self._set_moon_units()

    def set_model(self, variables: str, model_type, primary: str, sources_list: list[str]) -> None:
        """
        Set the model used for propagating the trajectory.
        It also initializes tr.model['data'] property.

        Parameters:
        -----------
        `variables` : str

        Variables in terms of which the equations will be chosen.

        `model_type`: str

        The model type of the dynamics.

        Options:

        1. r2bp (restricted two-body problem)

        2. cr3bp_fb (circular restricted three-body problem with center at the first primary body)

        3. cr3bp_sb (circular restricted three-body problem with center at the secondary primary body)

        4. nbp (generic n-body problem)

        `primary` : str

        The center body or the primaries wrt which the equations will be chosen.

        Options:

        1. Any string if model_type = 'r2bp'

        2. 'earth_moon', 'sun_earth' if model_type = 'cr3bp_fb' or 'cr3bp_sb'

        3. 'earth', 'moon' if model_type = 'nbp'

        `sources_lits` : list of str

        Listed sources of perturbations to be taken into account if model_type = 'nbp'.

        The full list of sources:

        'atm'       (Earth's atmosphere)

        'j2'        (Earth's J2)

        'srp'       (Solar radiation pressure)

        'sun'       (Gravitational acceleration of the Sun)

        'mercury'   (Gravitational acceleration of Mercury)

        'venus'     (Gravitational acceleration of Venus)

        'earth'     (Gravitational acceleration of the Earth)

        'mars'      (Gravitational acceleration of Mars)

        'jupiter'   (Gravitational acceleration of Jupiter)

        'saturn'    (Gravitational acceleration of Saturn)

        'uranus'    (Gravitational acceleration of Uranus)

        'neptune'   (Gravitational acceleration of Neptune)

        'moon'      (Gravitational acceleration of the Moon)

        'cmplxmoon' (Complex gravitational acceleration of the Moon)

        If `model_type` is 'cr3bp', then `model_specifics` is a dictionary.

        The dictionary contains t0 - the time at which the rotating and non-rotating coordinate systems coincide.

        """

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
            self._set_model_units(primary)
            if primary == 'earth':
                self.model['system'] = 'gcrs'
            elif primary == 'moon':
                self.model['system'] = 'scrs'
            else:
                raise Exception('Unknown primary')

        elif model_type == 'cr3bp_fb':

            if variables not in ['rv', 'rv_stm']:
                raise Exception('CR3BP model assumes rv or rv_stm variables.')

            if primary == 'earth_moon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self._set_model_units('earth_moon')
                self.model['system'] = 'rot_fb'
            elif primary == 'sun_earth':
                ku = kiam.units('Sun', 'Earth')
                self.model['data']['mass_parameter'] = ku['mu']
                self._set_model_units('sun_earth')
                self.model['system'] = 'rot_fb'
            else:
                raise Exception('Unknown primary.')

        elif model_type == 'cr3bp_sb':

            if variables not in ['rv', 'rv_stm']:
                raise Exception('CR3BP model assumes rv or rv_stm variables.')

            if primary == 'earth_moon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self._set_model_units('earth_moon')
                self.model['system'] = 'rot_sb'
            elif primary == 'sun_earth':
                ku = kiam.units('Sun', 'Earth')
                self.model['data']['mass_parameter'] = ku['mu']  # Sun - (Earth + Moon)
                self._set_model_units('sun_earth')
                self.model['system'] = 'rot_sb'
            else:
                raise Exception('Unknown primary.')

        elif model_type == 'br4bp_fb':

            if variables not in ['rv', 'rv_stm']:
                raise Exception('BR4BP model assumes rv or rv_stm variables.')

            if primary == 'earth_moon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self.model['data']['gm4b'] = 3.289005596e+05
                self.model['data']['a4b'] = 389.170375544352
                self.model['data']['theta0'] = 0.0
                self._set_model_units('earth_moon')
                self.model['system'] = 'rot_fb'
            else:
                raise Exception('Unknown primary.')

        elif model_type == 'br4bp_sb':

            if variables not in ['rv', 'rv_stm']:
                raise Exception('BR4BP model assumes rv or rv_stm variables.')

            if primary == 'earth_moon':
                ku = kiam.units('Earth', 'Moon')
                self.model['data']['mass_parameter'] = ku['mu']
                self.model['data']['gm4b'] = 3.289005596e+05
                self.model['data']['a4b'] = 389.170375544352  # 328 900·5596
                self.model['data']['theta0'] = 0.0
                self._set_model_units('earth_moon')
                self.model['system'] = 'rot_sb'

        elif model_type == 'hill':

            if variables not in ['rv']:
                raise Exception('Hill model assumes rv variables.')

            if primary == 'earth_moon':
                self._set_model_units('earth_moon')
                self.model['system'] = 'rot_sb'
            elif primary == 'sun_earth':
                self._set_model_units('sun_earth')
                self.model['system'] = 'rot_sb'

        elif model_type == 'nbp':

            self._set_model_sources()
            if variables == 'rv' and primary == 'earth':
                self._set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rv_stm' and primary == 'earth':
                self._set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rv' and primary == 'moon':
                self._set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'rv_stm' and primary == 'moon':
                self._set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'rvm' and primary == 'earth':
                self._set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'rvm' and primary == 'moon':
                self._set_model_units('moon')
                self.model['system'] = 'scrs'
            elif variables == 'ee' and primary == 'earth':
                self._set_model_units('earth')
                self.model['system'] = 'gcrs'
            elif variables == 'ee' and primary == 'moon':
                self._set_model_units('moon')
                self.model['system'] = 'scrs'
            else:
                raise Exception('Unknown model for the current variables and primary.')

        else:
            raise Exception('Unknown model type.')

    def propagate(self, tof: float, npoints: int = 2) -> None:
        """
        Propagate the trajectory.

        Parameters:
        -----------
        `tof` : float

        The time of flight in units specified by the chosen by set_model model.

        The time interval will be [t0, t0 + tof].

        `npoints` : int

        The number of nodes in the propagation time interval.
        """

        self.change_units(self.model['units']['name'])
        self.change_vars(self.model['vars'])
        self.change_system(self.model['system'])

        if self.vars in ['rv_stm', 'ee_stm', 'oe_stm']:
            stm = True
        else:
            stm = False

        tspan = numpy.linspace(self.times[-1], self.times[-1] + tof, npoints)

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
        elif self.model['type'] == 'br4bp_fb':
            T, X = kiam.propagate_br4bp(central_body='First', tspan=tspan, x0=self.states[0:, -1],
                                        mu=self.model['data']['mass_parameter'],
                                        gm4b=self.model['data']['gm4b'],
                                        a4b=self.model['data']['a4b'],
                                        theta0=self.model['data']['theta0'], stm=stm)
        elif self.model['type'] == 'br4bp_sb':
            T, X = kiam.propagate_br4bp(central_body='Secondary', tspan=tspan, x0=self.states[0:, -1],
                                        mu=self.model['data']['mass_parameter'],
                                        gm4b=self.model['data']['gm4b'],
                                        a4b=self.model['data']['a4b'],
                                        theta0=self.model['data']['theta0'], stm=stm)
        elif self.model['type'] == 'hill':
            T, X = kiam.propagate_hill(tspan=tspan, x0=self.states[0:, -1])
        else:
            raise Exception('Unknown model_type.')

        if len(self.parts) == 0:
            self.parts = [0, len(T)]
        else:
            self.parts.append(self.parts[-1] + len(T))

        self.jds = numpy.append(self.jds[0:-1], self.jds[-1] + (T - self.times[-1]) * self.units['TimeUnit'])
        self.times = numpy.append(self.times[0:-1], T)
        self.states = numpy.append(self.states[:, 0:-1], X, axis=1)
        self.finalDate = kiam.jd2time(self.jds[-1])

    def repropagate(self, tof: float, npoints: int = 2, start_index: int = 0) -> None:
        """
        Clears the calculated data in trajectory and propagate the trajectory from the beginning.

        Parameters:
        -----------
        `tof` : float

        The time of flight in units specified by the chosen by set_model model.

        The time interval will be [t0, t0 + tof].

        `npoints` : int

        The number of nodes in the propagation time interval.

        `start_index` : int

        The number of a time node from which the trajectory should be repropagated.

        At the moment, in can be only 0.
        """
        if start_index != 0:
            raise 'TBD.'
        self.clear()
        self.propagate(tof, npoints)

    def show(self, variables: str, draw=True):
        """
        Plots the specified characteristics of the trajectory.

        Parameters:
        -----------
        `variables` : str

        Variables that will be plotted.

        Options:

        1. 'xy', '3d' if self.vars in ['rv', 'rvm', 'rv_stm']

            'xy' plots the trajectory in x-y axes.

            '3d' plots the trajectory in 3d.

        2. 'a', 'e', 'inc', 'Om', 'w', 'th' if self.vars in ['oe', 'oem', 'oe_stm']

            'a' plots the semi-major axis wrt time

            'e' plots the eccentricity wrt time

            'inc' plots the inclination wrt time

            'Om' plots the right ascension of the ascending node wrt time

            'w' plots the argument of pericenter wrt time

            'th' plots the true anomaly wrt time

        3. 'h', 'ex', 'ey', 'ix', 'iy', 'L' if self.vars in ['ee', 'eem', 'ee_stm']

            'h' plots h = sqrt(p/mu) wrt time

            'ex' plots ex = e*cos(Omega+omega) wrt time

            'ey' plots ey = e*sin(Omega+omega) wrt time

            'ix' plots ix = tan(i/2)*cos(Omega) wrt time

            'iy' plots iy = tan(i/2)*sin(Omega) wrt time

            'L' plots L = theta + omega + Omega wrt time

        `draw` : bool

        If True (by default), the fig plotly object will be returned and figure will be showed.

        If False, the fig plot object will be returned and the figure will not be showed.

        Returns:
        --------
        `ax` : matplotlib axis object

        The matplotlib axis object for further work.
        """
        if self.units_name == 'dim':
            tlabel = 'Time of flight, days'
        else:
            tlabel = 'Time of flight, nondimensional'
        if self.vars in ['rv', 'rvm', 'rv_stm']:
            if variables == 'xy':
                if self.units_name == 'earth':
                    xlabel = r'$x\text{, Earth radii}$'
                    ylabel = r'$y\text{, Earth radii}$'
                elif self.units_name == 'moon':
                    xlabel = r'$x\text{, Moon radii}$'
                    ylabel = r'$y\text{, Moon radii}$'
                elif self.units_name == 'dim':
                    xlabel = r'$x\text{, km}$'
                    ylabel = r'$y\text{, km}$'
                elif self.units_name in ['earth_moon', 'sun_earth']:
                    xlabel = r'$x\text{, nondimensional}$'
                    ylabel = r'$y\text{, nondimensional}$'
                else:
                    raise Exception('Unknown units.')
                fig = kiam.plot(self.states[0, :], self.states[1, :], xlabel=xlabel, ylabel=ylabel)
            elif variables == '3d':
                if self.units_name == 'earth':
                    xlabel = r'x, Earth radii'
                    ylabel = r'y, Earth radii'
                    zlabel = r'z, Earth radii'
                elif self.units_name == 'moon':
                    xlabel = r'x, Moon radii'
                    ylabel = r'y, Moon radii'
                    zlabel = r'z, Moon radii'
                elif self.units_name == 'dim':
                    xlabel = r'x, km}'
                    ylabel = r'y, km}'
                    zlabel = r'y, km}'
                elif self.units_name in ['earth_moon', 'sun_earth']:
                    xlabel = r'x, nondimensional'
                    ylabel = r'y, nondimensional'
                    zlabel = r'z, nondimensional'
                else:
                    raise Exception('Unknown units.')
                fig = kiam.plot3(self.states[0, :], self.states[1, :], self.states[2, :], xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
            else:
                raise 'Unknown variables to show.'
        elif self.vars in ['oe', 'oem', 'oe_stm']:
            if variables == 'a':
                if self.units_name == 'dim':
                    ylabel = 'Semi-major axis, km'
                elif self.units_name == 'earth':
                    ylabel = 'Semi-major axis, Earth''s radii'
                else:
                    ylabel = ''
                fig = kiam.plot(self.times, self.states[0, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'e':
                ylabel = 'Eccentricity'
                fig = kiam.plot(self.times, self.states[1, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'inc':
                ylabel = 'Inclination, degrees'
                fig = kiam.plot(self.times, self.states[2, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'Om':
                ylabel = 'Right ascension of the ascending node, degrees'
                fig = kiam.plot(self.times, self.states[3, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'w':
                ylabel = 'Argument of pericenter, degrees'
                fig = kiam.plot(self.times, self.states[4, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            elif variables == 'th':
                ylabel = 'True anomaly, degrees'
                fig = kiam.plot(self.times, self.states[5, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            else:
                raise 'Unknown classical orbital element. Elements: a, e, inc, Om, w, th.'
        elif self.vars in ['ee', 'eem', 'ee_stm']:
            if variables == 'h':
                if self.units_name == 'dim':
                    ylabel = r'$h\text{, (km/s)}^{-1}$'
                elif self.units_name == 'earth':
                    ylabel = r'$h\text{, nondimensional}$'
                else:
                    ylabel = ''
                fig = kiam.plot(self.times, self.states[0, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ex':
                ylabel = '$e_x$'
                fig = kiam.plot(self.times, self.states[1, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ey':
                ylabel = '$e_y$'
                fig = kiam.plot(self.times, self.states[2, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'ix':
                ylabel = '$i_x$'
                fig = kiam.plot(self.times, self.states[3, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'iy':
                ylabel = '$i_y$'
                fig = kiam.plot(self.times, self.states[4, :], xlabel=tlabel, ylabel=ylabel)
            elif variables == 'L':
                ylabel = 'True longitude, degrees'
                fig = kiam.plot(self.times, self.states[5, :] / math.pi * 180, xlabel=tlabel, ylabel=ylabel)
            else:
                raise 'Unknown equinoctial orbital element. Elements: h, ex, ey, ix, iy, L.'
        else:
            raise 'Unknown variables.'
        if draw:
            fig.show()
        return fig

    def copy(self):
        """
        Returns independent copy of the trajectory object.
        """
        return copy.deepcopy(self)

    def clear(self) -> None:
        """
        Clears states, times, julian dates, parts, and resets the final date.
        """
        self.states = numpy.reshape(self.states[:, 0], (6, 1))
        self.times = numpy.reshape(self.times[0], (1,))
        self.jds = numpy.reshape(self.jds[0], (1,))
        self.parts = []
        self.finalDate = self.initialDate

    def change_vars(self, new_vars: str) -> None:
        """
        Change variables to the specified ones.

        Parameters:
        -----------
        `new_vars` : str

        The name of the new variables.

        It can be 'rv', 'rvm', 'rv_stm', 'ee', 'eem', 'ee_stm', 'oe', 'oem', 'oe_stm'.

        All transformations are possible with given rules:

        1. units_name are 'earth' or 'moon'

        2. transformations that increases the number of variables are impossible.

        E.g. rv -> rvm or rv -> rv_stm transformations are impossible.

        The routine automatically find the chain of transformations from the current variables to
        the specified variables.
        """
        if self.vars == new_vars:
            return
        if new_vars not in self.vars_graph.nodes:
            raise Exception('Unknown new_vars.')
        p = nx.shortest_path(self.vars_graph, self.vars, new_vars)
        for i in range(len(p) - 1):
            self._vars_transform(p[i], p[i + 1])

    def change_system(self, new_system: str) -> None:
        """
        Change coordinate system to the specified one.

        Parameters:
        -----------
        `new_system` : str

        The name of the new coordinate system.

        Options:

        'scrs', 'mer', 'sors', 'ssrm_em' (the Moon-centered systems)

        'gcrs', 'itrs', 'gsrf_em', 'gsrf_se' (the Earth-centered systems)

        'hcrs', 'hsrf_se' (the Sun-centered systems)

        'ine_fb' (non-rotating frame in CR3BP, BR4BP at first primary)

        'ine_cm' (non-rotating frame in CR3BP, BR4BP at barycenter)

        'ine_sb' (non-rotating frame in CR3BP, BR4BP, Hill at secondary primary)

        'rot_fb' (rotating frame in CR3BP, BR4BP at first primary)

        'rot_cm' (rotating frame in CR3BP, BR4BP at barycenter)

        'rot_sb' (rotating frame in CR3BP, BR4BP, Hill at secondary primary)

        All transformations are possible within given rules:

        1. The variables are 'rv', 'rvm', or 'rv_stm'.

        2. Changes with variables 'rv_stm' are possible only in scrs <-> sors and scrs <-> mer translations.

        Implemented transformations:

        Model nbp, both systems should be from: systems_ephemeris.

        Model CR3BP, both systems should be from: systems_cr3bp.

        Model BR4BP, both systems should be from: systems_br4bp.

        Model Hill, both systems should be from: systems_hill.

        The routine automatically find the chain of transformations from the current coordinate system to
        the specified coordinate system.
        """
        if self.system == new_system:
            return
        if new_system not in self.systems_graph.nodes:
            raise Exception('Unknown new_system.')
        p = nx.shortest_path(self.systems_graph, self.system, new_system)
        for i in range(len(p) - 1):
            self._system_transform(p[i], p[i + 1])

    def change_units(self, new_units: str) -> None:
        """
        Change units to the specified ones.

        Parameters:
        -----------
        `new_units` : str

        The name of the new units.

        Options:

        'earth', 'moon', 'dim', 'earth_moon', 'sun_earth'.

        All transformations are possible.

        The routine automatically find the chain of transformations from the current units to
        the specified units.
        """
        if self.units_name == new_units:
            return
        if new_units not in self.units_graph.nodes:
            raise Exception('Unknown new_units.')
        p = nx.shortest_path(self.units_graph, self.units_name, new_units)
        for i in range(len(p) - 1):
            self._units_transform(p[i], p[i + 1])

    # Variables transformations.
    def _allocate_vars_graph(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Allocates the graph of variables for further automating transformations.
        """

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

    def _vars_transform(self, vars1: str, vars2: str) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Call the routine that transforms variables from one to another.
        Routine checks the rules for the transformation.

        Parameters:
        -----------
        `vars1` : str

        Variables before transformation.

        `vars2` : str

        Variables after transformation.
        """
        if vars1 == 'rv' and vars2 == 'ee':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv2ee suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm.')
            self.states = kiam.rv2ee(self.states, 1.0)
            # for i in range(self.states.shape[1]):
            #    self.states[0:6, i] = kiam.rv2ee(self.states[0:6, i], 1.0)
            self.vars = 'ee'
        elif vars1 == 'ee' and vars2 == 'rv':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: ee2rv suggests earth or moon, mu = 1.0.')
            elif self.vars != 'ee' and self.vars != 'eem':
                raise Exception('Vars should be ee or eem.')
            self.states = kiam.ee2rv(self.states, 1.0)
            # for i in range(self.states.shape[1]):
            #    self.states[0:6, i] = kiam.ee2rv(self.states[0:6, i], 1.0)
            self.vars = 'rv'
        elif vars1 == 'rv' and vars2 == 'oe':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: rv2oe suggests earth or moon, mu = 1.0.')
            elif self.vars != 'rv' and self.vars != 'rvm':
                raise Exception('Vars should be rv or rvm.')
            self.states = kiam.rv2oe(self.states, 1.0)
            # for i in range(self.states.shape[1]):
            #    self.states[0:6, i] = kiam.rv2oe(self.states[0:6, i], 1.0)
            self.vars = 'oe'
        elif vars1 == 'oe' and vars2 == 'rv':  # mu = 1.0
            if self.units_name != 'earth' and self.units_name != 'moon':
                raise Exception('Wrong units: oe2rv suggests earth or moon, mu = 1.0.')
            elif self.vars != 'oe' and self.vars != 'oem':
                raise Exception('Vars should be oe or oem.')
            self.states = kiam.oe2rv(self.states, 1.0)
            # for i in range(self.states.shape[1]):
            #    self.states[0:6, i] = kiam.oe2rv(self.states[0:6, i], 1.0)
            self.vars = 'rv'
        elif vars1 == 'rvm' and vars2 == 'eem':
            self._vars_transform('rv', 'ee')
            self.vars = 'eem'
        elif vars1 == 'eem' and vars2 == 'rvm':
            self._vars_transform('ee', 'rv')
            self.vars = 'rvm'
        elif vars1 == 'rvm' and vars2 == 'oem':
            self._vars_transform('rv', 'oe')
            self.vars = 'oem'
        elif vars1 == 'oem' and vars2 == 'rvm':
            self._vars_transform('oe', 'rv')
            self.vars = 'rvm'
        elif vars1 == 'rv_stm' and vars2 == 'rv':
            numpy.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'rv'
        elif vars1 == 'oe_stm' and vars2 == 'oe':
            numpy.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'oe'
        elif vars1 == 'ee_stm' and vars2 == 'ee':
            numpy.delete(self.states, [i for i in range(6, 42)], 0)
            self.vars = 'ee'
        elif vars1 == 'rvm' and vars2 == 'rv':
            numpy.delete(self.states, 6, 0)
            self.vars = 'rv'
        elif vars1 == 'oem' and vars2 == 'oe':
            numpy.delete(self.states, 6, 0)
            self.vars = 'oe'
        elif vars1 == 'eem' and vars2 == 'ee':
            numpy.delete(self.states, 6, 0)
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
                phi_rv = numpy.reshape(self.states[6:42, i], (6, 6)).T
                phi_oe = kiam.dotainvb(numpy.matmul(doe, phi_rv), doe0)
                self.states[6:42, i] = numpy.reshape(phi_oe.T, (36,))
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
                phi_oe = numpy.reshape(self.states[6:42, i], (6, 6)).T
                phi_rv = kiam.dotainvb(numpy.matmul(drv, phi_oe), drv0)
                self.states[6:42, i] = numpy.reshape(phi_rv.T, (36,))
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
                phi_rv = numpy.reshape(self.states[6:42, i], (6, 6)).T
                phi_ee = kiam.dotainvb(numpy.matmul(dee, phi_rv), dee0)
                self.states[6:42, i] = numpy.reshape(phi_ee.T, (36,))
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
                phi_ee = numpy.reshape(self.states[6:42, i], (6, 6)).T
                phi_rv = kiam.dotainvb(numpy.matmul(drv, phi_ee), drv0)
                self.states[6:42, i] = numpy.reshape(phi_rv.T, (36,))
            self.vars = 'rv_stm'
        else:
            raise Exception('Unknown variable transformaton.')

    # System transformations.
    def _allocate_systems_graph(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Allocates the graph of coordinate systems for further automating transformations.
        """

        # Ephemeris coordinate systems
        self.systems_graph.add_edge('itrs', 'gcrs')
        self.systems_graph.add_edge('gcrs', 'gsrf_em')
        self.systems_graph.add_edge('gcrs', 'gsrf_se')
        self.systems_graph.add_edge('scrs', 'ssrf_em')
        self.systems_graph.add_edge('gcrs', 'scrs')
        self.systems_graph.add_edge('scrs', 'mer')
        self.systems_graph.add_edge('gcrs', 'hcrs')
        self.systems_graph.add_edge('hcrs', 'hsrf_se')
        self.systems_graph.add_edge('scrs', 'sors')

        # CR3BP and BR4BP coordinate systems
        self.systems_graph.add_edge('ine_fb', 'rot_fb')
        self.systems_graph.add_edge('ine_sb', 'rot_sb')
        self.systems_graph.add_edge('ine_cm', 'rot_cm')
        self.systems_graph.add_edge('rot_fb', 'ine_fb')
        self.systems_graph.add_edge('rot_sb', 'ine_sb')
        self.systems_graph.add_edge('rot_cm', 'ine_cm')
        self.systems_graph.add_edge('rot_fb', 'rot_cm')
        self.systems_graph.add_edge('rot_cm', 'rot_fb')
        self.systems_graph.add_edge('rot_fb', 'rot_sb')
        self.systems_graph.add_edge('rot_sb', 'rot_fb')
        self.systems_graph.add_edge('rot_sb', 'rot_cm')
        self.systems_graph.add_edge('rot_cm', 'rot_sb')

    def _system_transform(self, system1: str, system2: str) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Call the routine that transforms coordinate systems from one to another.
        Routine checks the rules for the transformation.

        Parameters:
        -----------
        `system1` : str

        Coordinate system before transformation.

        `system2` : str

        Coordinate system after transformation.
        """

        case_nbp = self.model['type'] == 'nbp' and system1 in systems_ephemeris and system2 in systems_ephemeris
        if case_nbp:
            if system1 == 'itrs' and system2 == 'gcrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'itrs':
                    raise Exception('System should be itrs.')
                self.states[0:6] = kiam.itrs2gcrs(self.states[0:6], self.jds)
                # self.states[0:3, :] = kiam.itrs2gcrs(self.states[0:3, :], self.jds)
                # self.states[3:6, :] = kiam.itrs2gcrs(self.states[3:6, :], self.jds)
                # for i in range(self.states.shape[1]):
                #    self.states[0:3, i] = kiam.itrs2gcrs(self.states[0:3, i], self.jds[i])
                #    self.states[3:6, i] = kiam.itrs2gcrs(self.states[3:6, i], self.jds[i])
                self.system = 'gcrs'
            elif system1 == 'gcrs' and system2 == 'itrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'gcrs':
                    raise Exception('System should be gcrs.')
                self.states[0:6] = kiam.gcrs2itrs(self.states[0:6], self.jds)
                # self.states[0:3, :] = kiam.gcrs2itrs(self.states[0:3, :], self.jds)
                # self.states[3:6, :] = kiam.gcrs2itrs(self.states[3:6, :], self.jds)
                # for i in range(self.states.shape[1]):
                #    self.states[0:3, i] = kiam.gcrs2itrs(self.states[0:3, i], self.jds[i])
                #    self.states[3:6, i] = kiam.gcrs2itrs(self.states[3:6, i], self.jds[i])
                self.system = 'itrs'
            elif system1 == 'gcrs' and system2 == 'gsrf_em':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'gcrs':
                    raise Exception('System should be gcrs.')
                self.states[0:6, :] = kiam.ine2rot_eph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'gsrf_em'
            elif system1 == 'gsrf_em' and system2 == 'gcrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'gsrf_em':
                    raise Exception('System should be gsrf_em.')
                self.states[0:6, :] = kiam.rot2ine_eph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'gcrs'
            elif system1 == 'gcrs' and system2 == 'gsrf_se':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'gcrs':
                    raise Exception('System should be gcrs.')
                self.states[0:6, :] = kiam.ine2rot_eph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'gsrf_se'
            elif system1 == 'gsrf_se' and system2 == 'gcrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'gsrf_se':
                    raise Exception('System should be gsrf_se.')
                self.states[0:6, :] = kiam.rot2ine_eph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'gcrs'
            elif system1 == 'scrs' and system2 == 'ssrf_em':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'scrs':
                    raise Exception('System should be scrs.')
                self.states[0:6, :] = kiam.ine2rot_eph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'ssrf_em'
            elif system1 == 'ssrf_em' and system2 == 'scrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'ssrf_em':
                    raise Exception('System should be ssrf_em.')
                self.states[0:6, :] = kiam.rot2ine_eph(self.states[0:6, :], self.jds, 'Earth', 'Moon',
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
                self.states[0:6, :] = kiam.ine2rot_eph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
                                                       self.units['DistUnit'], self.units['VelUnit'])
                self.system = 'hsrf_se'
            elif system1 == 'hsrf_se' and system2 == 'hcrs':
                if self.vars != 'rv' and self.vars != 'rvm':
                    raise Exception('Vars should be rv or rvm')
                elif self.system != 'hsrf_se':
                    raise Exception('System should be hsrf_se.')
                self.states[0:6, :] = kiam.rot2ine_eph(self.states[0:6, :], self.jds, 'Sun', 'Earth',
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
                        phi_scrs = numpy.reshape(self.states[6:42, i], (6, 6)).T
                        phi_sors = kiam.dotainvb(numpy.matmul(dxsors[:, :, i], phi_scrs), dxsors[:, :, 0])
                        self.states[6:42, i] = numpy.reshape(phi_sors.T, (36,))
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
                        phi_sors = numpy.reshape(self.states[6:42, i], (6, 6)).T
                        phi_scrs = kiam.dotainvb(numpy.matmul(dxscrs[:, :, i], phi_sors), dxscrs[:, :, 0])
                        self.states[6:42, i] = numpy.reshape(phi_scrs.T, (36,))
                self.system = 'scrs'
            return

        case_cr3bp = self.model['type'] in models_cr3bp and system1 in systems_cr3bp and system2 in systems_cr3bp
        case_br4bp = self.model['type'] in models_br4bp and system1 in systems_br4bp and system2 in systems_br4bp
        if case_cr3bp or case_br4bp:
            if case_cr3bp:
                if self.vars not in ['rv', 'rvm']:
                    raise Exception('Vars should be rv or rvm.')
                if self.model['primary'] not in ['sun_earth', 'earth_moon']:
                    raise Exception("Model's primary should be sun_earth or earth_moon.")
                if self.units_name != self.model['primary']:
                    raise Exception('Units_name should equal primary.')
            elif case_br4bp:
                if self.vars not in ['rv', 'rvm']:
                    raise Exception('Vars should be rv or rvm.')
                if self.model['primary'] not in ['earth_moon']:
                    raise Exception("Model's primary should be earth_moon.")
                if self.units_name != self.model['primary']:
                    raise Exception('Units_name should equal primary.')
            if (system1 == 'ine_fb' and system2 == 'rot_fb') or \
                    (system1 == 'ine_sb' and system2 == 'rot_sb') or \
                    (system1 == 'ine_cm' and system2 == 'rot_cm'):
                t0 = self.model['data'].get('t0')
                if t0 is None:
                    raise Exception('Please set t0 to self.model["data"]')
                self.states[0:6, :] = kiam.ine2rot(self.states[0:6, :], self.times, t0)
                self.system = system2
            elif (system1 == 'rot_fb' and system2 == 'ine_fb') or \
                    (system1 == 'rot_sb' and system2 == 'ine_sb') or \
                    (system1 == 'rot_cm' and system2 == 'ine_cm'):
                t0 = self.model['data'].get('t0')
                if t0 is None:
                    raise Exception('Please set t0 to self.model["data"]')
                self.states[0:6, :] = kiam.rot2ine(self.states[0:6, :], self.times, t0)
                self.system = system2
            elif system1 == 'rot_fb' and system2 == 'rot_sb':
                self.states[0] = self.states[0] - 1.0
                self.system = system2
            elif system1 == 'rot_fb' and system2 == 'rot_cm':
                self.states[0] = self.states[0] - self.model['data']['mass_parameter']
                self.system = system2
            elif system1 == 'rot_sb' and system2 == 'rot_fb':
                self.states[0] = self.states[0] + 1.0
                self.system = system2
            elif system1 == 'rot_sb' and system2 == 'rot_cm':
                self.states[0] = self.states[0] + 1.0 - self.model['data']['mass_parameter']
                self.system = system2
            elif system1 == 'rot_cm' and system2 == 'rot_fb':
                self.states[0] = self.states[0] + self.model['data']['mass_parameter']
                self.system = system2
            elif system1 == 'rot_cm' and system2 == 'rot_sb':
                self.states[0] = self.states[0] - 1.0 + self.model['data']['mass_parameter']
                self.system = system2
            return

        case_hill = self.model['type'] in models_hill and system1 in systems_hill and system2 in systems_hill
        if case_hill:
            if self.vars not in ['rv', 'rvm']:
                raise Exception('Vars should be rv or rvm.')
            if self.model['primary'] not in ['earth_moon', 'sun_earth']:
                raise Exception("Model's primary should be earth_moon.")
            if self.units_name != self.model['primary']:
                raise Exception('Units_name should equal primary.')
            if system1 == 'ine_sb' and system2 == 'rot_sb':
                t0 = self.model['data'].get('t0')
                if t0 is None:
                    raise Exception('Please set t0 to self.model["data"]')
                self.states[0:6, :] = kiam.ine2rot(self.states[0:6, :], self.times, t0)
                self.system = system2
            elif system1 == 'rot_sb' and system2 == 'ine_sb':
                t0 = self.model['data'].get('t0')
                if t0 is None:
                    raise Exception('Please set t0 to self.model["data"]')
                self.states[0:6, :] = kiam.rot2ine(self.states[0:6, :], self.times, t0)
                self.system = system2
            return

        raise Exception('Possible relations:\n'
                        f'Model: nbp, both systems should be from: {", ".join(systems_ephemeris)}.\n'
                        f'Model: CR3BP, both systems should be from: {", ".join(systems_cr3bp)}.\n'
                        f'Model: BR4BP, both systems should be from: {", ".join(systems_br4bp)}.\n'
                        f'Model: Hill, both systems should be from: {", ".join(systems_hill)}.\n')

    # Units transformations and settings.
    def _allocate_units_graph(self):
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Allocates the graph of units for further automating transformations.
        """
        self.units_graph.add_edge('dim', 'earth')
        self.units_graph.add_edge('dim', 'moon')
        self.units_graph.add_edge('dim', 'earth_moon')
        self.units_graph.add_edge('dim', 'sun_earth')

    def _units_transform(self, units1: str, units2: str) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Call the routine that transforms units from one to another.
        Routine checks the rules for the transformation.

        Parameters:
        -----------
        `units1` : str

        Units before transformation.

        `units2` : str

        Units after transformation.
        """
        if units1 == 'dim' and units2 == 'earth':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self._set_earth_units()
                self._undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self._set_earth_units()
                self._undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self._set_earth_units()
                self._undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'earth'
        elif units1 == 'dim' and units2 == 'moon':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self._set_moon_units()
                self._undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self._set_moon_units()
                self._undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self._set_moon_units()
                self._undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'moon'
        elif units1 == 'dim' and units2 == 'earth_moon':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self._set_earth_moon_units()
                self._undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self._set_earth_moon_units()
                self._undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self._set_earth_moon_units()
                self._undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'earth_moon'
        elif units1 == 'dim' and units2 == 'sun_earth':
            if self.units_name != 'dim':
                raise Exception('Units should be dim.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self._set_sun_earth_units()
                self._undim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self._set_sun_earth_units()
                self._undim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self._set_sun_earth_units()
                self._undim_oe()
            else:
                raise Exception('Unknown vars.')
            self.units_name = 'sun_earth'
        elif units2 == 'dim':
            if units1 != self.units_name:
                raise Exception(f'Not {units1} units.')
            if self.vars in ['rv', 'rvm', 'rv_stm']:
                self._dim_rv()
            elif self.vars in ['ee', 'eem', 'ee_stm']:
                self._dim_ee()
            elif self.vars in ['oe', 'oem', 'oe_stm']:
                self._dim_oe()
            else:
                raise Exception('Unknown vars.')
            self._set_dim_units()
            self.units_name = 'dim'

    def _undim_rv(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Undimensionalize the time and position-velocity states
        """
        self.times = self.times / self.units['TimeUnit']
        self.states[0:3, :] = self.states[0:3, :] / self.units['DistUnit']
        self.states[3:6, :] = self.states[3:6, :] / self.units['VelUnit']

    def _dim_rv(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Dimensionalize the time and position-velocity states
        """
        self.times = self.times * self.units['TimeUnit']
        self.states[0:3, :] = self.states[0:3, :] * self.units['DistUnit']
        self.states[3:6, :] = self.states[3:6, :] * self.units['VelUnit']

    def _undim_ee(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Undimensionalize the time and equinoctial orbits elements in states.
        """
        self.times = self.times / self.units['TimeUnit']
        self.states[0, :] = self.states[0, :] * self.units['VelUnit']

    def _dim_ee(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Dimensionalize the time and equinoctial orbits elements in states.
        """
        self.times = self.times * self.units['TimeUnit']
        self.states[0, :] = self.states[0, :] / self.units['VelUnit']

    def _undim_oe(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Undimensionalize the time and classical orbits elements in states.
        """
        self.times = self.times / self.units['TimeUnit']
        self.states[0, :] = self.states[0, :] / self.units['DistUnit']

    def _dim_oe(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Dimensionalize the time and classical orbits elements in states.
        """
        self.times = self.times * self.units['TimeUnit']
        self.states[0, :] = self.states[0, :] * self.units['DistUnit']

    def _set_earth_units(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Set coefficients for transformations to/from Earth units.
        """
        ku = kiam.units('Earth')
        self.units['mu'] = 1.0
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']  # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']  # m/s^2

    def _set_moon_units(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Set coefficients for transformations to/from Moon units.
        """
        ku = kiam.units('Moon')
        self.units['mu'] = 1.0
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']  # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']  # m/s^2

    def _set_earth_moon_units(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Set coefficients for transformations to/from Earth--Moon units.
        """
        ku = kiam.units('Earth', 'Moon')
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']  # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']  # m/s^2

    def _set_sun_earth_units(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Set coefficients for transformations to/from Sun--Earth units.
        """
        ku = kiam.units('Sun', 'Earth')
        self.units['DistUnit'] = ku['DistUnit']  # km
        self.units['VelUnit'] = ku['VelUnit']  # km/s
        self.units['TimeUnit'] = ku['TimeUnit']  # days
        self.units['AccUnit'] = ku['AccUnit']  # m/s^2

    def _set_dim_units(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Set coefficients for transformations to/from dimensional units.
        """
        ku = kiam.units('Earth')
        self.units['mu'] = ku['GM']  # km^3/s^2
        self.units['DistUnit'] = 1.0  # km
        self.units['VelUnit'] = 1.0  # km/s
        self.units['TimeUnit'] = 1.0  # days
        self.units['AccUnit'] = 1.0  # m/s^2

    # Auxilary model routines.
    def _set_model_units(self, units_name: str) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Sets the units used in model of the trajectory.

        Parameters:
        -----------
        `units_name` : str

        Name of the units.

        Options:
        'earth', 'moon', 'earth_moon', 'sun_earth'
        """

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

    def _set_model_sources(self) -> None:
        """
        FOR THE TOOLBOX DEVELOPERS ONLY.
        Sets the sources taken into account in the model.
        """

        self.model['sources'] = kiam.prepare_sources_dict()

        for source in self.model['sources_list']:
            self.model['sources'][source.lower()] = True


def traj2dict(tr: Trajectory) -> dict:
    """
    Converts trajectory object to dictionary.

    Parameters:
    -----------
    `tr` : Trajectory

    The Trajectory object needed to convert into dictionary.

    Returns:
    --------
    `d` : dict

    The dictionary containing the attributes of the given Trajectory object.
    """

    d = {'vars': tr.vars, 'states': tr.states, 'times': tr.times,
         'system': tr.system, 'units_name': tr.units_name, 'jds': tr.jds,
         'initialDate': tr.initialDate, 'finalDate': tr.finalDate,
         'units': tr.units, 'parts': tr.parts, 'model': tr.model}
    return d


def dict2traj(d: dict) -> Trajectory:
    """
    Converts dictionary to trajectory object.

    Parameters:
    -----------
    `d` : dict

    The dictionary containing the attributes of some Trajectory object.

    Returns:
    --------
    `tr` : Trajectory

    The Trajectory object converted from the dictionary.
    """
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
