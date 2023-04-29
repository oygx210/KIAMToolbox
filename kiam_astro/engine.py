"""
This Python module is a part of the KIAM Astrodynamics Toolbox developed in
Keldysh Institute of Applied Mathematics (KIAM), Moscow, Russia.

The module provides the Engine base class, which can be used to implement its own
classes abstracting the concept of an engine in astrodynamics. Based on the Engine
class, the module implements several classes for abstracting existing engines.
Users can create their own classes by inheriting from the Engine base class or
from the classes presented in this module.

The toolbox is licensed under the MIT License.

The GitHub page of the project:
https://github.com/shmaxg/KIAMToolbox.

Sources: https://fakel-russia.com/en/productions
"""
from typing import Union

class Engine:

    def __init__(self):
        """
        Initializes the instance of the Engine base class.
        The object properties:

        `name` : str

        The name of the engine.

        `engine_class` : str

        The class of the eingine. Options: 'attitude', 'orbital'.

        `engine_type` : str

        The type of the eingine. Options: 'chemical', 'electric'.

        `force` : float, None

        The force in N. Can be varied by implementing the get_force method.

        `specific_impulse` : float, None

        The specific impulse in s. Can be varied by implementing the get_specific_impulse method.

        `input_power` : float, None

        The input power in W. Can be varied by implementing the get_input_power method.

        `efficiency` : float, None

        The efficiency in [0, 1]. Can be varied by implementing the get_efficiency method.

        """

        self.name: str = ''
        self.engine_class: str = ''  # 'attitude', 'orbital'
        self.engine_type: str = ''  # 'chemical', 'electric'

        self.force: Union[float, None] = None  # N
        self.specific_impulse: Union[float, None] = None  # s
        self.input_power: Union[float, None] = None  # W
        self.efficiency: Union[float, None] = None  # in [0, 1]

    def get_force(self):
        """
        Returns the force of the engine in N. Can be overrided.

        Returns:
        --------

        `force` : float

        The force of the engine in N.
        """
        return self.force

    def get_specific_impulse(self):
        """
        Returns the specific impulse of the engine in s. Can be overrided.

        Returns:
        --------

        `specific_impulse` : float

        The specific impulse of the engine in s.
        """
        return self.specific_impulse

    def get_input_power(self):
        """
        Returns the input power of the engine in W. Can be overrided.

        Returns:
        --------

        `input_power` : float

        The input power of the engine in W.
        """
        return self.input_power

    def get_efficiency(self):
        """
        Returns the efficiency of the engine in [0, 1]. Can be overrided.

        Returns:
        --------

        `efficiency` : float

        The efficiency of the engine in [0, 1].
        """
        return self.efficiency

    def describe(self):
        dims = {
            'force': 'N',
            'specific_impulse': 's',
            'discharge_voltage': 'V',
            'discharge_current': 'A',
            'discharge_power': 'W',
            'power_to_thrust_ratio': 'W/N',
            'min_lifetime_h': 'days',
            'mass': 'kg',
            'dimensions': 'mm'
        }
        for a in dir(self):
            if not a.startswith('__') and not callable(getattr(self, a)):
                dimension = dims.get(a)
                dimension = '' if dimension is None else dimension
                print(f'{a} = {eval(f"self.{a}")} {dimension}')


class SPT50(Engine):
    """
    The class abstracting the SPT-50 engine: https://fakel-russia.com/produkciya
    This is a regime of the SPT-50M engine.
    """

    def __init__(self):
        super(SPT50, self).__init__()

        self.name = 'spt50'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'
        self.discharge_voltage = 180.0  # V
        self.discharge_current = 1.25  # A
        self.discharge_power = 225.0  # W

        self.force = 14.0e-03  # N
        self.specific_impulse = 860.0  # s
        self.power_to_thrust_ratio = 16.1e+03  # W/N
        self.min_lifetime_h = 1217 / 24  # days
        self.min_lifetime_cycles = 3011  # cycles
        self.mass = 1.23  # kg
        self.dimensions = '160 x 120 x 91'  # mm


class SPT50M(Engine):
    """
    The class abstracting the SPT-50M engine: https://fakel-russia.com/en/productions
    """

    def __init__(self, option: str = ''):
        """
        Parameters:
        -----------

        `option` : str

        The engine's regime.

        Options: 'Xe_low' (low force, Xenon), 'Xe_high' (high force, Xenon)

        """
        super(SPT50M, self).__init__()

        self.name = 'spt50m'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        if option == 'Xe_low':

            self.propellant = 'Xe'

            self.discharge_voltage = 180.0  # V
            self.discharge_current = 1.25  # A
            self.discharge_power = 225.0  # W

            self.force = 14.8e-03  # N
            self.specific_impulse = 930.0  # s
            self.power_to_thrust_ratio = 15.2e+03  # W/N
            self.min_lifetime_h = 5000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 1.32  # kg
            self.dimensions = '169 x 120 x 88'  # mm

        elif option == 'Xe_high':

            self.propellant = 'Xe'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 1.0  # A
            self.discharge_power = 300.0  # W

            self.force = 18.0e-03  # N
            self.specific_impulse = 1250.0  # s
            self.power_to_thrust_ratio = 16.7e+03  # W/N
            self.min_lifetime_h = 5000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 1.32  # kg
            self.dimensions = '169 x 120 x 88'  # mm


class SPT70(Engine):
    """
    The class abstracting the SPT-70 engine: https://fakel-russia.com/en/productions
    """

    def __init__(self):
        super(SPT70, self).__init__()

        self.name = 'spt70'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300.0  # V
        self.discharge_current = 2.23  # A
        self.discharge_power = 670.0  # W

        self.force = 39.0e-03  # N
        self.specific_impulse = 1470.0  # s
        self.power_to_thrust_ratio = 16.1e+03  # W/N
        self.min_lifetime_h = 3100 / 24  # days
        self.min_lifetime_cycles = 3000  # cycles
        self.mass = 1.5  # kg
        self.dimensions = '198 x 146 x 98'  # mm


class SPT70M(Engine):
    """
    The class abstracting the SPT-70M engine: https://fakel-russia.com/en/productions
    """

    def __init__(self, option: str = ''):
        """
        Parameters:
        -----------

        `option` : str

        The engine's regime.

        Options: 'Xe_low' (low force, Xenon), 'Xe_med' (medium force, Xenon),
        'Xe_high' (high force, Xenon), 'Kr_low' (low force, Krypton),
        'Kr_med' (medium force, Krypton), 'Kr_high' (high force, Krypton)

        """
        super(SPT70M, self).__init__()

        self.name = 'spt70m'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        if option == 'Xe_low':

            self.propellant = 'Xe'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 2.00  # A
            self.discharge_power = 600.0  # W

            self.force = 36.0e-03  # N
            self.specific_impulse = 1430.0  # s
            self.power_to_thrust_ratio = 15.2e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Xe_med':

            self.propellant = 'Xe'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 2.67  # A
            self.discharge_power = 800.0  # W

            self.force = 48.0e-03  # N
            self.specific_impulse = 1530.0  # s
            self.power_to_thrust_ratio = 16.7e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Xe_high':

            self.propellant = 'Xe'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 3.33  # A
            self.discharge_power = 1000.0  # W

            self.force = 59.0e-03  # N
            self.specific_impulse = 1600.0  # s
            self.power_to_thrust_ratio = 16.9e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_low':

            self.propellant = 'Kr'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 2.00  # A
            self.discharge_power = 600.0  # W

            self.force = 28.0e-03  # N
            self.specific_impulse = 1380.0  # s
            self.power_to_thrust_ratio = 21.4e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_med':

            self.propellant = 'Kr'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 2.67  # A
            self.discharge_power = 800.0  # W

            self.force = 37.0e-03  # N
            self.specific_impulse = 1490.0  # s
            self.power_to_thrust_ratio = 21.6e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        elif option == 'Kr_high':

            self.propellant = 'Kr'

            self.discharge_voltage = 300.0  # V
            self.discharge_current = 3.33  # A
            self.discharge_power = 1000.0  # W

            self.force = 47.0e-03  # N
            self.specific_impulse = 1560.0  # s
            self.power_to_thrust_ratio = 21.3e+03  # W/N
            self.min_lifetime_h = 7000 / 24  # days
            self.min_lifetime_cycles = 11000  # cycles
            self.mass = 2.6  # kg
            self.dimensions = '200 x 128 x 94'  # mm

        else:

            raise Exception('Unknown option.')


class SPT100B(Engine):
    """
    The class abstracting the SPT-100B engine: https://fakel-russia.com/en/productions
    """

    def __init__(self):
        super(SPT100B, self).__init__()

        self.name = 'spt100B'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300.0  # V
        self.discharge_current = 4.5  # A
        self.discharge_power = 1350.0  # W

        self.force = 83.0e-03  # N
        self.specific_impulse = 1540.0  # s
        self.efficiency = 0.45
        self.power_to_thrust_ratio = 16.3e+03  # W/N
        self.min_lifetime_h = 9000 / 24  # days
        self.min_lifetime_cycles = 8800  # cycles
        self.mass = 3.5  # kg
        self.dimensions = '225 x 150 x 125'  # mm


class SPT100BM(Engine):
    """
    The class abstracting the SPT-100BM engine: https://fakel-russia.com/en/productions
    """

    def __init__(self):
        super(SPT100BM, self).__init__()

        self.name = 'spt100BM'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300.0  # V
        self.discharge_current = 4.5  # A
        self.discharge_power = 1350.0  # W

        self.force = 90.0e-03  # N
        self.specific_impulse = 1600.0  # s
        self.efficiency = 0.52
        self.power_to_thrust_ratio = 15.0e+03  # W/N
        self.min_lifetime_h = 9000 / 24  # days
        self.min_lifetime_cycles = 9000  # cycles
        self.mass = 4.2  # kg
        self.dimensions = '200 x 142 x 110'  # mm


class SPT140D(Engine):
    """
    The class abstracting the SPT-140D engine: https://fakel-russia.com/en/productions
    """

    def __init__(self):
        super(SPT140D, self).__init__()

        self.name = 'spt140D'
        self.engine_class = 'orbital'
        self.engine_type = 'electric'

        self.propellant = 'Xe'

        self.discharge_voltage = 300.0  # V
        self.discharge_current = 15.0  # A
        self.discharge_power = 4500.0  # W

        self.force = 290.0e-03  # N
        self.specific_impulse = 1750.0  # s
        self.efficiency = 0.53
        self.power_to_thrust_ratio = 15.5e+03  # W/N
        self.min_lifetime_h = 15000 / 24  # days
        self.min_lifetime_cycles = 7000  # cycles
        self.mass = 8.5  # kg
        self.dimensions = '305 x 249 x 109'  # mm
