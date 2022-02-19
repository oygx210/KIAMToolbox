import kiam


class Model:

    def __init__(self, variables, model_type, primary, sources_cell):

        self.data = {}
        self.info = {}
        self.units = {}
        self.sources = {}
        self.data = {}

        primary = primary.lower()

        self.vars = variables
        self.type = model_type
        self.primary = primary
        self.data['sources_cell'] = sources_cell

        if model_type == 'r2bp':

            if (primary != 'earth') and (primary != 'moon'):
                raise Exception('Earth or Moon as primary are required for r2bp type of model.')

            self.data['grav_parameter'] = 1.0
            self.set_units(primary)
            if primary == 'earth':
                self.system = 'gcrs'
            elif primary == 'moon':
                self.system = 'scrs'
            else:
                raise Exception('Unknown primary')

            self.eqs = lambda t, s: kiam.r2bp(t, s)

            self.info['ode_name'] = 'r2bp'
            self.info['arguments'] = [self.data['grav_parameter'], False]

        elif model_type == 'cr3bp_fb':

            if primary == 'earthmoon':
                ku = kiam.units('Earth', 'Moon')
                self.data['mass_parameter'] = ku['mu']
                self.set_units('earth_moon')
                self.system = 'gsrf_em'
            elif primary == 'sunearth':
                ku = kiam.units('Sun', 'Earth')
                self.data['mass_parameter'] = ku['mu']
                self.set_units('sun_earth')
                self.system = 'hsrf_se'
            else:
                raise Exception('Unknown primary.')

            self.eqs = lambda t, s: kiam.cr3bp_fb(t, s, ku.mu, False)
            self.info['ode_name'] = 'kcr3bp_fb'
            self.info['arguments'] = [self.data['mass_parameter'], False]

        elif model_type == 'cr3bp_sb':

            if primary == 'earthmoon':
                ku = kiam.units('Earth', 'Moon')
                self.data['mass_parameter'] = ku['mu']
                self.set_units('earth_moon')
                self.system = 'ssrf_em'
            elif primary == 'sunearth':
                ku = kiam.units('Sun', 'Earth')
                self.data['mass_parameter'] = ku['mu']  # Sun - (Earth + Moon)
                self.set_units('sun_earth')
                self.system = 'gsrf_se'
            else:
                raise Exception('Unknown primary.')

            self.eqs = lambda t, s: kiam.cr3bp_sb(t, s, ku.mu, False)

            self.info['ode_name'] = 'cr3bp_sb'
            self.info['arguments'] = [self.data['mass_parameter'], False]

        elif model_type == 'nbp':

            self.set_sources()
            if variables == 'rv' and primary == 'earth':
                self.set_units('earth')
                self.system = 'gcrs'
                self.eqs = lambda t, s: self.nbp_rv_earth(t, s, False)
                self.info['ode_name'] = 'nbp_rv_earth'
                self.info['arguments'] = [self.sources, self.data, self.units, True]
            elif variables == 'rv_stm' and primary == 'earth':
                self.set_units('earth')
                self.system = 'gcrs'
                self.eqs = lambda t, s: self.nbp_rv_earth(t, s, True)
                self.info['ode_name'] = 'nbp_rv_earth'
                self.info['arguments'] = [self.sources, self.data, self.units, True]
            elif variables == 'rv' and primary == 'moon':
                self.set_units('moon')
                self.system = 'scrs'
                self.eqs = lambda t, s: self.nbp_rv_moon(t, s, False)
                self.info['ode_name'] = 'nbp_rv_moon'
                self.info['arguments'] = [self.sources, self.data, self.units, False]
            elif variables == 'rv_stm' and primary == 'moon':
                self.set_units('moon')
                self.system = 'scrs'
                self.eqs = lambda t, s: self.nbp_rv_moon(t, s, True)
                self.info['ode_name'] = 'nbp_rv_moon'
                self.info['arguments'] = [self.sources, self.data, self.units, True]
            elif variables == 'rvm' and primary == 'earth':
                self.set_units('earth')
                self.system = 'gcrs'
                self.eqs = lambda t, s: self.nbp_rvm_earth(t, s, False)
            elif variables == 'rvm' and primary == 'moon':
                self.set_units('moon')
                self.system = 'scrs'
                self.eqs = lambda t, s: self.nbp_rvm_moon(t, s, False)
            elif variables == 'ee' and primary == 'earth':
                self.set_units('earth')
                self.system = 'gcrs'
                self.eqs = lambda t, s: self.nbp_ee_earth(t, s, False)
            elif variables == 'ee' and primary == 'moon':
                self.set_units('moon')
                self.system = 'scrs'
                self.eqs = lambda t, s: self.nbp_ee_moon(t, s, False)
            else:
                raise Exception('Unknown model.')

    def set_units(self, units_name):

        if units_name == 'earth':

            self.units['name'] = 'earth'
            units = kiam.units('Earth')
            self.units['mu'] = units['GM']

        elif units_name == 'moon':

            self.units['name'] = 'moon'
            units = kiam.units('Moon')
            self.units['mu'] = units['GM']

        elif units_name == 'earth_moon':

            self.units['name'] = 'earth_moon'
            units = kiam.units('Earth', 'Moon')

        elif units_name == 'sun_earth':

            self.units['name'] = 'sun_earth'
            units = kiam.units('Sun', 'Earth')

        else:

            raise Exception('Unknown units_name.')

        self.units['DistUnit'] = units['DistUnit']  # km
        self.units['VelUnit'] = units['VelUnit']  # km/s
        self.units['TimeUnit'] = units['TimeUnit']  # days
        self.units['AccUnit'] = units['AccUnit']  # m/s^2

        _, star, planet, moon, _ = kiam.astro_const()

        self.units['RSun'] = star['Sun']['MeanRadius'] / units['DistUnit']  # 695700 km
        self.units['REarth'] = planet['Earth']['EquatorRadius'] / units['DistUnit']  # 6378.1366 km
        self.units['RMoon'] = moon['Moon']['MeanRadius'] / units['DistUnit']  # 1737.4 km

    def set_sources(self):

        self.sources = kiam.prepare_sources_dict()

        for source in self.data['sources_cell']:
            self.sources[source.lower()] = True

    def nbp_rv_earth(self, t, s, stm_req):
        return kiam.nbp_rv_earth(t, s, stm_req, self.sources, self.data, self.units)

    def nbp_rv_moon(self, t, s, stm_req):
        return kiam.nbp_rv_moon(t, s, stm_req, self.sources, self.data, self.units)

    def nbp_rvm_earth(self, t, s, stm_req):
        return kiam.nbp_rvm_earth(t, s, stm_req, self.sources, self.data, self.units)

    def nbp_rvm_moon(self, t, s, stm_req):
        return kiam.nbp_rvm_moon(t, s, stm_req, self.sources, self.data, self.units)

    def nbp_ee_earth(self, t, s, stm_req):
        return kiam.nbp_ee_earth(t, s, stm_req, self.sources, self.data, self.units)

    def nbp_ee_moon(self, t, s, stm_req):
        return kiam.nbp_ee_moon(t, s, stm_req, self.sources, self.data, self.units)
