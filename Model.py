import kiam


class Model:

    def __init__(self, variables, model_type, primary, sources_cell, jd_zero):

        self.data = {}
        self.info = {}
        self.units = {}

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

