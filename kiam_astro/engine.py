"""

Sources: https://fakel-russia.com/en/productions
"""

class Engine:

    def __init__(self):

        self.name = None
        self.engine_class = None  # 'attitude', 'orbital'
        self.engine_type = None  # 'chemical', 'electric'

        self.force = None
        self.specific_impulse = None
        self.input_power = None
        self.efficiency = None

    def get_force_dependence(self):
        return self.force

    def get_specific_impulse_dependence(self):
        return self.specific_impulse

    def get_input_power_dependence(self):
        return self.input_power

    def get_efficiency_dependence(self):
        return self.efficiency

class SPT50(Engine):

    def __init__(self):
        super(SPT50, self).__init__()

        self.name = 'spt50'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'
        self.discharge_voltage = 180  # V
        self.discharge_current = 1.25  # A
        self.discharge_power = 225  # W

        self.force = 14.0e-03  # N
        self.specific_impulse = 860  # s
        self.power_to_thrust_ratio = 16.1e+03  # W/N
        self.min_lifetime_h = 1217/24  # days
        self.min_lifetime_cycles = 3011  # cycles
        self.mass = 1.23  # kg
        self.dimensions = '160 x 120 x 91'  # mm

class SPT50M(Engine):

    def __init__(self, option=None):
        super(SPT50M, self).__init__()

        self.name = 'spt50m'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        if option == 'Xe_low':

            self.propellant = 'Xe'

            self.discharge_voltage = 180  # V
            self.discharge_current = 1.25  # A
            self.discharge_power = 225  # W

            self.force = 14.8e-03  # N
            self.specific_impulse = 930  # s
            self.power_to_thrust_ratio = 15.2e+03  # W/N
            self.min_lifetime_h = 5000/24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 1.32  # kg
            self.dimensions = '169 x 120 x 88'  # mm

        elif option == 'Xe_high':

            self.propellant = 'Xe'

            self.discharge_voltage = 300  # V
            self.discharge_current = 1.0  # A
            self.discharge_power = 300  # W

            self.force = 18.0e-03  # N
            self.specific_impulse = 1250  # s
            self.power_to_thrust_ratio = 16.7e+03  # W/N
            self.min_lifetime_h = 5000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 1.32  # kg
            self.dimensions = '169 x 120 x 88'  # mm

class SPT70(Engine):

    def __init__(self):
        super(SPT70, self).__init__()

        self.name = 'spt70'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300  # V
        self.discharge_current = 2.23  # A
        self.discharge_power = 670  # W

        self.force = 39.0e-03  # N
        self.specific_impulse = 1470  # s
        self.power_to_thrust_ratio = 16.1e+03  # W/N
        self.min_lifetime_h = 3100 / 24  # days
        self.min_lifetime_cycles = 3000  # cycles
        self.mass = 1.5  # kg
        self.dimensions = '198 x 146 x 98'  # mm

class SPT70M(Engine):

    def __init__(self, option=None):
        super(SPT70M, self).__init__()

        self.name = 'spt70m'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        if option == 'Xe_low':

            self.propellant = 'Xe'

            self.discharge_voltage = 300  # V
            self.discharge_current = 2.00  # A
            self.discharge_power = 600  # W

            self.force = 36.0e-03  # N
            self.specific_impulse = 1430  # s
            self.power_to_thrust_ratio = 15.2e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Xe_med':

            self.propellant = 'Xe'

            self.discharge_voltage = 300  # V
            self.discharge_current = 2.67  # A
            self.discharge_power = 800  # W

            self.force = 48.0e-03  # N
            self.specific_impulse = 1530  # s
            self.power_to_thrust_ratio = 16.7e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Xe_high':

            self.propellant = 'Xe'

            self.discharge_voltage = 300  # V
            self.discharge_current = 3.33  # A
            self.discharge_power = 1000  # W

            self.force = 59.0e-03  # N
            self.specific_impulse = 1600  # s
            self.power_to_thrust_ratio = 16.9e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_low':

            self.propellant = 'Kr'

            self.discharge_voltage = 300  # V
            self.discharge_current = 2.00  # A
            self.discharge_power = 600  # W

            self.force = 28.0e-03  # N
            self.specific_impulse = 1380  # s
            self.power_to_thrust_ratio = 21.4e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_med':

            self.propellant = 'Kr'

            self.discharge_voltage = 300  # V
            self.discharge_current = 2.67  # A
            self.discharge_power = 800  # W

            self.force = 37.0e-03  # N
            self.specific_impulse = 1490  # s
            self.power_to_thrust_ratio = 21.6e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_high':

            self.propellant = 'Kr'

            self.discharge_voltage = 300  # V
            self.discharge_current = 3.33  # A
            self.discharge_power = 1000  # W

            self.force = 47.0e-03  # N
            self.specific_impulse = 1560  # s
            self.power_to_thrust_ratio = 21.3e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        else:

            raise Exception('Unknown option.')

class SPT100B(Engine):

    def __init__(self):
        super(SPT100B, self).__init__()

        self.name = 'spt100B'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300  # V
        self.discharge_current = 4.5  # A
        self.discharge_power = 1350  # W

        self.force = 83.0e-03  # N
        self.specific_impulse = 1540  # s
        self.efficiency = 0.45
        self.power_to_thrust_ratio = 16.3e+03  # W/N
        self.min_lifetime_h = 9000 / 24  # days
        self.min_lifetime_cycles = 8800  # cycles
        self.mass = 3.5  # kg
        self.dimensions = '225 x 150 x 125'  # mm

class SPT100BM(Engine):

    def __init__(self):
        super(SPT100BM, self).__init__()

        self.name = 'spt100BM'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300  # V
        self.discharge_current = 4.5  # A
        self.discharge_power = 1350  # W

        self.force = 90.0e-03  # N
        self.specific_impulse = 1600  # s
        self.efficiency = 0.52
        self.power_to_thrust_ratio = 15.0e+03  # W/N
        self.min_lifetime_h = 9000 / 24  # days
        self.min_lifetime_cycles = 9000  # cycles
        self.mass = 4.2  # kg
        self.dimensions = '200 x 142 x 110'  # mm

class SPT140D(Engine):

    def __init__(self):
        super(SPT140D, self).__init__()

        self.name = 'spt140D'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300  # V
        self.discharge_current = 15  # A
        self.discharge_power = 4500  # W

        self.force = 290.0e-03  # N
        self.specific_impulse = 1750  # s
        self.efficiency = 0.53
        self.power_to_thrust_ratio = 15.5e+03  # W/N
        self.min_lifetime_h = 15000 / 24  # days
        self.min_lifetime_cycles = 7000  # cycles
        self.mass = 8.5  # kg
        self.dimensions = '305 x 249 x 109'  # mm
