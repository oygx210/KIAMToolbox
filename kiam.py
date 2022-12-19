"""
This Python module is a part of the KIAM Astrodynamics Toolbox developed in
Keldysh Institute of Applied Mathematics (KIAM), Moscow, Russia.

The module serves as a safe and convenient interface to Fortran-compiled
astrodynamical routines and provides instruments for performing translations
between variables, coordinate systems, and time descriptions, propagating the
trajectories in various models, and getting fast answers on typical
questions about the two and n-body problems. It also contains some plotting
routins and useful matrix linear algebra operations.

The toolbox is licensed under the MIT License.

The GitHub page of the project:
https://github.com/shmaxg/KIAMToolbox.
"""

import FKIAMToolbox as fkt
import jdcal
import datetime
import math
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
from typing import Union, Any

mpl.rcParams['figure.dpi'] = 150
# mpl.use('Qt5Agg')  # troubleshoots on some systems

# General mathematics and tools (documented without examples)
def invadotb(a: numpy.ndarray, b: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate `a^(-1)*b` for matrices `a` and `b`.

    Parameters:
    -----------
    `a` : numpy.ndarray, shape (n, n)

    A square matrix.

    `b` : numpe.ndarray, shape (n, n)

    A square matrix.

    Returns:
    --------
    `c` : numpy.ndarray, shape (n, n)

    The matrix that equals `a^(-1)*b`.
    """
    return numpy.linalg.solve(a, b)
def dotainvb(a: numpy.ndarray, b: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate `a*b^(-1)` for matrices `a` and `b`.

    Parameters:
    -----------
    `a` : numpy.ndarray, shape (n, n)

    A square matrix.

    `b` : numpe.ndarray, shape (n, n)

    A square matrix.

    Returns:
    --------
    `c` : numpy.ndarray, shape (n, n)

    The matrix that equals `a*b^(-1)`
    """
    c = numpy.linalg.solve(b.T, a.T)
    return c.T
def eye2vec(n: int) -> numpy.ndarray:
    """
    Vector form of an identity matrix.

    Parameters:
    -----------
    `n` : int

    The number of rows and columns in the identity matrix.

    Returns:
    --------
    `a` : numpy.ndarray, shape (n**2,)

    Vector form of the identity matrix.
    """
    return numpy.reshape(numpy.eye(n), (n**2,))
def mat2vec(a: numpy.ndarray) -> numpy.ndarray:
    """
    Square matrix to vector form translation.

    Parameters:
    -----------
    a : numpy.ndarray, shape (n, n)

    A square matrix.

    Returns:
    --------
    v : numpy.ndarray, shape (n**2,)

    Vector form of the matrix.

    Vector structure (Fortran/MATLAB order): [a11, a21, a31, ... ]

    Examples:
    ---------
    ```
    >> a = numpy.array([[1, 2], [3, 4]])

    >> v = kiam.mat2vec(a)

    [1 3 2 4]
    ```
    """
    v = numpy.reshape(a, (a.size,), order='F')
    return v
def vec2mat(v: numpy.ndarray) -> numpy.ndarray:
    """
    Vector to square matrix translation.

    Parameters:
    -----------
    v : numpy.ndarray, shape (n**2,)

    Vector.

    Returns:
    --------

    a : numpy.ndarray, shape (n, n)

    A square matrix.

    Matrix structure (Fortran/MATLAB order): [[v1, v2, ..., vn], [v_(n+1), ...]].T

    Examples:
    ---------
    ```
    >> v = numpy.array([1, 2, 3, 4])

    >> m = kiam.vec2mat(v)

    [[1 3]

    [2 4]]
    ```
    """
    a = numpy.reshape(v, (int(round(numpy.sqrt(v.size))), -1), order='F')
    return a
def to_float(*args: Any) -> tuple:
    """
    Convert all arguments to the float64 type.

    Parameters:
    -----------
    `*args`

    Arguments separated by comma to convert to float64.

    Returns:
    --------
    `float_args` : tuple

    Converted to float64 arguments.
    """
    args_float = tuple(numpy.array(arg, dtype='float64') for arg in args)
    return args_float
def sind(x: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Sine of a degree argument.

    Parameters:
    -----------
    `x` : float, numpy.ndarray

    Angle or an array of angles in degrees.

    Returns:
    --------
    `s` : float, numpy.ndarray

    A sine or array of sines of angles in degrees.
    """
    return numpy.sin(x/180*numpy.pi)
def cosd(x: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Cosine of a degree argument.

    Parameters:
    -----------
    `x` : float, numpy.ndarray

    Angle or an array of angles in degrees.

    Returns:
    --------
    `s` : float, numpy.ndarray

    A cosine or array of cosines of angles in degrees.
    """
    return numpy.cos(x/180*numpy.pi)
def tand(x: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Tangent of a degree argument.

    Parameters:
    -----------
    `x` : float, numpy.ndarray

    Angle or an array of angles in degrees.

    Returns:
    --------
    `s` : float, numpy.ndarray

    A tangent or array of tangents of angles in degrees.
    """
    return numpy.tan(x/180*numpy.pi)
def cotand(x: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Coangent of a degree argument.

    Parameters:
    -----------
    `x` : float, numpy.ndarray

    Angle or an array of angles in degrees.

    Returns:
    --------
    `s` : float, numpy.ndarray

    A cotangent or array of cotangents of angles in degrees.
    """
    return 1/numpy.tan(x/180*numpy.pi)

# Plotting functions (not documented, will be rewritten)
def plot(x, y, ax=None, style='-', xlabel='', ylabel='', linewidth=1.0, show=False, saveto=False):
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, style, linewidth=linewidth)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.5, linestyle=':')
    # ax.axis('equal')
    if saveto:
        plt.savefig(saveto, dpi=300)
    if show:
        plt.show()
    return ax
def plot3(x, y, z, ax=None, style='-', xlabel='', ylabel='', zlabel='', linewidth=1.0, show=False, saveto=False):
    if ax is None:
        ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, style, linewidth=linewidth)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.xaxis._axinfo["grid"]['linestyle'] = ":"
    ax.yaxis._axinfo["grid"]['linestyle'] = ":"
    ax.zaxis._axinfo["grid"]['linestyle'] = ":"
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    if saveto:
        plt.savefig(saveto, dpi=300)
    if show:
        plt.show()
    return ax
def polar_plot(phi, r, ax=None, rmax=None, style='-', linewidth=1.0, show=False, saveto=False):
    if ax is None:
        _, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(phi, r, style, linewidth=linewidth)
    if rmax is not None:
        ax.set_rmax(rmax)
    ax.grid(True, alpha=0.5, linestyle=':')
    if saveto:
        plt.savefig(saveto, dpi=300)
    if show:
        plt.show()
    return ax
def histogram(x, num_bins=None, density=False, ax=None, xlabel='', ylabel='', show=False, saveto=False):
    if ax is None:
        _, ax = plt.subplots()
    n, bins, patches = ax.hist(x, bins=num_bins, density=density)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.5, linestyle=':')
    if saveto:
        plt.savefig(saveto, dpi=300)
    if show:
        plt.show()
    return ax, {'n': n, 'bins': bins, 'patches': patches}
def boxplot(x, ax=None, xlabel='', ylabel='', show=False, saveto=False):
    if ax is None:
        _, ax = plt.subplots()
    data = ax.boxplot(x, vert=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.5, linestyle=':')
    if saveto:
        plt.savefig(saveto, dpi=300)
    if show:
        plt.show()
    return ax, data
def save_plot(saveto):
    plt.savefig(saveto, dpi=300)

# Ephemeris information (documented with examples)
def jd2time(jd: float) -> datetime.datetime:
    """
    Julian date to usual date and time.

    Parameters:
    -----------
    `jd` : float

    Julian date

    Returns:
    --------
    `time` : datetime.datetime

    Date and time object of type datetime.datetime

    Examples:
    ---------
    ```
    >> jd2time(2459905.5)

    2022-11-22 00:00:00
    ```
    """
    gcal = jdcal.jd2gcal(2400000.5, jd - 2400000.5)
    frac_hours = gcal[3] * 24
    hours = math.floor(frac_hours)
    frac_minutes = (frac_hours - hours) * 60
    minutes = math.floor(frac_minutes)
    frac_seconds = (frac_minutes - minutes) * 60
    seconds = int(round(frac_seconds, 0))
    if seconds != 60:
        return datetime.datetime(gcal[0], gcal[1], gcal[2], hours, minutes, seconds)
    else:
        return datetime.datetime(gcal[0], gcal[1], gcal[2], hours, minutes) + datetime.timedelta(minutes=1)
def time2jd(time: datetime.datetime) -> float:
    """
    Usual date and time to Julian date.

    Parameters:
    -----------
    `time` : datetime.datetime

    Date and time object of type datetime.datetime

    Returns:
    --------
    `jd` : float

    Julian date

    Examples:
    ---------
    ```
    >> time2jd(datetime.datetime(2022, 11, 22, 0, 0, 0, 0))

    2459905.5
    ```
    """
    return sum(jdcal.gcal2jd(time.year, time.month, time.day)) + time.hour / 24 + \
        time.minute / 1440 + time.second / 86400 + (time.microsecond / 1000000) / 86400
def juliandate(year: int, month: int, day: int, hour: int, minute: int, second: int) -> float:
    """
    Usual date to Julian date.

    Parameters:
    -----------
    `year` : int

    Year

    `month` : int

    Month

    `day` : int

    Day

    `hour` : int

    Hour

    `minute` : int

    Minute

    `second` : int

    Second

    Returns:
    --------
    `jd` : float

    Julian date

    Examples:
    ---------
    ```
    >> juliandate(2022, 11, 22, 0, 0, 0)

    2459905.5
    ```
    """
    # return fkt.ephemeris.juliandate(year, month, day, hour, minute, second)
    return time2jd(datetime.datetime(year, month, day, hour, minute, second, 0))
def planet_state(jd: float, center: str, target: str) -> numpy.ndarray:
    """
    Gives position and velocity of the planet at specified julian date
    wrt to the specified center (planet).

    Parameters:
    -----------
    `jd` : float

    Julian date

    `center` : str

    Name of the center planet

    `target` : str

    Name of the target planet

    Returns:
    --------
    `state` : numpy.ndarray, shape(6,)

    State of the target planet wrt the center planet.

    Position in km, velocity in km/s.

    Examples:
    ---------
    ```

    >> planet_state(kiam.juliandate(2022, 12, 3, 0, 0, 0), 'Earth', 'Moon')

    [ 3.76623766e+05  7.07472988e+04  1.01213236e+04
        -1.36269070e-01 8.97864551e-01  4.72492325e-01]
    ```
    """

    return fkt.ephemeris.planetstate(jd, target, center)

# Translations (documented with examples)
def deg2rad(deg: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Degrees to radians conversion.

    Parameters:
    -----------
    `deg` : float, numpy.ndarray

    Angle or array of angles in degrees.

    Returns:
    --------
    `rad` : float, numpy.ndarray

    Angle or array of angles in radians.

    Examples:
    ---------
    ```
    >> deg2rad(180)

    3.141592653589793
    ```
    """
    return deg/180*numpy.pi
def rad2deg(rad: Union[float, numpy.ndarray]) -> Union[float, numpy.ndarray]:
    """
    Radians to degrees conversion.

    Parameters:
    -----------
    `rad` : float, numpy.ndarray

    Angle or array of angles in radians.

    Returns:
    --------
    `deg` : float, numpy.ndarray

    Angle or array of angles in degrees.

    Examples:
    ---------
    ```
    >> rad2deg(3.141592)

    179.99996255206332
    ```
    """
    return rad/numpy.pi*180
def rv2oe(rv: numpy.ndarray, mu: float, grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Position and velocity to classical orbital elements.

    Parameters:
    -----------
    `rv` : numpy.ndarray, shape (6,), (6,n)

    6D phase vector or array of column 6D phase vectors containing position and velocity.

    Vector structure: [x, y, z, vx, dy, dz].

    `mu` : float

    Gravitational parameter.

    `grad_req` : bool

    Flag to calculate the derivatives of elements wrt position and velocity.

    Returns:
    --------
    `oe` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column vectors of classical orbital elements:

    a (semi-major axis),

    e (eccentricity),

    i (inclination),

    Omega (right ascension of the ascending node),

    omega (argument of pericenter),

    theta (true anomaly)

    `doe` : numpy.ndarray, shape (6,6), (6,6,n)

    6x6 matrix or 6x6xn array of partial derivatives of oe wrt rv (doe/drv).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    rv = numpy.array([1, 0, 0, 0.1, 1, 0.1])

    oe = kiam.rv2oe(rv, 1.0, False)

    oe, doe = kiam.rv2oe(rv, 1.0, True)
    ```
    """
    if rv.shape == (6,):
        out = fkt.translations.krv2oe(rv, mu, grad_req)
    elif len(rv.shape) == 2 and rv.shape[0] == 6:
        fkt.translations.rv_mat = rv
        fkt.translations.krv2oe_mat(mu, grad_req)
        out = (fkt.translations.oe_mat, fkt.translations.doe_mat)
    else:
        raise 'rv should be a 6D vector or a 6xn array of 6D column vectors.'
    return _return_if_grad_req(out, grad_req)
def oe2rv(oe: numpy.ndarray, mu: float, grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Classical orbital elements to position and velocity.

    Parameters:
    -----------
    `oe` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column vectors of classical orbital elements:

    a (semi-major axis),

    e (eccentricity),

    i (inclination),

    Omega (right ascension of the ascending node),

    omega (argument of pericenter),

    theta (true anomaly).

    `mu` : float

    Gravitational parameter.

    `grad_req` : bool

    Flag to calculate the derivatives of position and velocity wrt elements.

    Returns:
    --------
    `rv` : numpy.ndarray, shape (6,), (6,n)

    6D phase vector or array of 6D column vectors containing position and velocity.

    Vector structure: [x, y, z, vx, dy, dz].

    `drv` : numpy.ndarray, shape (6,6), (6,6,n)

    6x6 matrix or 6x6xn array of partial derivatives of rv wrt oe (drv/doe).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    oe = numpy.array([1, 0.1, 1.0, 0.0, 0.0, 0.0])

    rv = kiam.oe2rv(oe, 1.0, False)

    rv, drv = kiam.oe2rv(oe, 1.0, True)
    ```
    """
    if oe.shape == (6,):
        out = fkt.translations.koe2rv(oe, mu, grad_req)
    elif len(oe.shape) == 2 and oe.shape[0] == 6:
        fkt.translations.oe_mat = oe
        fkt.translations.koe2rv_mat(mu, grad_req)
        out = (fkt.translations.rv_mat, fkt.translations.drv_mat)
    else:
        raise 'oe should be a 6D vector or a 6xn array of 6D column vectors.'
    return _return_if_grad_req(out, grad_req)
def rv2ee(rv: numpy.ndarray, mu: float, grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Position and velocity to equinoctial orbital elements.

    Parameters:
    -----------
    `rv` : numpy.ndarray, shape (6,), (6,)

    6D phase vector or array of 6D column vectors containing position and velocity.

    Vector structure: [x, y, z, vx, dy, dz]

    `mu` : float

    Gravitational parameter

    `grad_req` : bool

    Flag to calculate the derivatives of elements wrt position and velocity

    Returns:
    --------
    `ee` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column vectors of equinoctial orbital elements:

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination.

    `dee` : numpy.ndarray, shape (6,6), (6,6,n)

    6x6 matrix or 6x6xn array of partial derivatives of ee wrt rv (dee/drv).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    rv = numpy.array([1, 0, 0, 0, 1, 0])

    ee = kiam.rv2ee(rv, 1.0, False)

    ee, dee = kiam.rv2ee(rv, 1.0, True)
    ```
    """
    if rv.shape == (6,):
        out = fkt.translations.krv2ee(rv, mu, grad_req)
    elif len(rv.shape) == 2 and rv.shape[0] == 6:
        fkt.translations.rv_mat = rv
        fkt.translations.krv2ee_mat(mu, grad_req)
        out = (fkt.translations.ee_mat, fkt.translations.dee_mat)
    else:
        raise 'rv should be a 6D vector or a 6xn array of 6D column vectors.'
    return _return_if_grad_req(out, grad_req)
def ee2rv(ee: numpy.ndarray, mu: float, grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Equinoctial orbital elements to position and velocity.

    Parameters:
    -----------
    `ee` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column vectors of equinoctial orbital elements:

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination.

    `mu` : float

    Gravitational parameter.

    `grad_req` : bool

    Flag to calculate the derivatives of position and velocity wrt elemets.

    Returns:
    --------
    `rv` : numpy.ndarray, shape (6,), (6,n)

    6D phase vector or array of 6D column vectors containing position and velocity.

    Vector structure: [x, y, z, vx, dy, dz].

    `drv` : numpy.ndarray, shape (6,6), (6,6,n)

    6x6 matrix or 6x6xn array of partial derivatives of rv wrt ee (drv/dee).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    ee = numpy.array([1, 0, 0, 0, 0, 0])

    rv = kiam.ee2rv(ee, 1.0, False)

    rv, drv = kiam.ee2rv(ee, 1.0, True)
    ```
    """
    if ee.shape == (6,):
        out = fkt.translations.kee2rv(ee, mu, grad_req)
    elif len(ee.shape) == 2 and ee.shape[0] == 6:
        fkt.translations.ee_mat = ee
        fkt.translations.kee2rv_mat(mu, grad_req)
        out = (fkt.translations.rv_mat, fkt.translations.drv_mat)
    else:
        raise 'ee should be a 6D vector or a 6xn array of 6D column vectors.'
    return _return_if_grad_req(out, grad_req)
def cart2sphere(cart: numpy.ndarray) -> numpy.ndarray:
    """
    Cartesian coordinates to spherical coordinates.

    Parameters:
    -----------
    `cart` : numpy.ndarray, shape (3,), (3, n)

    3D vector or 3xn array of column 3D vectors of Cartesian coordinates

    Vector structure: [x, y, z]

    Returns:
    --------
    `sphere` : numpy.ndarray, shape (3,), (3, n)

    3D vector or 3xn array of column 3D vectors of spherical coordinates

    Vector structure: [r, phi, theta], where

    phi in [-pi, pi],

    theta in [0, pi],

    x = r*cos(theta)*cos(phi),

    y = r*cos(theta)*sin(phi),

    z = r*sin(theta).

    Examples:
    ---------
    ```
    >> cart = numpy.array([1, 0, 0])

    >> sphere = kiam.cart2sphere(cart)

    [1.         0.         1.57079633]
    ```
    """
    if cart.shape == (3,):
        out = fkt.translations.kcart2sphere(numpy.reshape(cart, (3, 1)))
        return out[:, 0]
    elif len(cart.shape) == 2 and cart.shape[0] == 3:
        return fkt.translations.kcart2sphere(cart)
    else:
        raise 'cart should be a 3D vector or a 3xn array of vectors.'
def sphere2cart(sphere: numpy.ndarray) -> numpy.ndarray:
    """
    Spherical coordinates to Cartesian coordinates.

    Parameters:
    -----------
    `sphere` : numpy.ndarray, shape (3,), (3, n)

    3D vector or 3xn array of column 3D vectors of spherical coordinates

    Vector structure: [r, phi, theta], where

    phi in [-pi, pi],

    theta in [0, pi],

    x = r*cos(theta)*cos(phi),

    y = r*cos(theta)*sin(phi),

    z = r*sin(theta)

    Returns:
    --------
    `cart` : numpy.ndarray, shape (3,), (3, n)

    3D vector or 3xn array of column 3D vectors of Cartesian coordinates.

    Vector structure: [x, y, z].

    Examples:
    ---------
    ```
    >> sphere = numpy.array([1, 0, 0])

    >> cart = kiam.sphere2cart(sphere)

    [0. 0. 1.]
    ```
    """
    if sphere.shape == (3,):
        out = fkt.translations.ksphere2cart(numpy.reshape(sphere, (3, 1)))
        return out[:, 0]
    elif len(sphere.shape) == 2 and sphere.shape[0] == 3:
        return fkt.translations.ksphere2cart(sphere)
    else:
        raise 'sphere should be a 3D vector or a 3xn array of vectors.'
def cart2latlon(cart: numpy.ndarray) -> numpy.ndarray:
    """
    Cartesian coordinates to latitude and longitude.

    Parameters:
    -----------
    `cart` : numpy.ndarray, shape (3,), (3, n)

    3D vector or array of column 3D vectors of Cartesian coordinates.

    Vector structure: [x, y, z]

    Returns:
    --------
    `latlon` : numpy.ndarray, shape (2,), (2, n)

    2D Vector or array of column 2D vectors of latitude and longitude pairs.

    Vector structure: [lat, lon], where

    lat in [-pi/2, pi/2],

    lon in [-pi, pi].

    Examples:
    ---------
    ```
    >> cart = numpy.array([1, 0, 0])

    >> latlon = kiam.cart2latlon(cart)

    [0. 0.]
    ```
    """
    if cart.shape[0] != 3 or len(cart.shape) not in [1, 2]:
        raise 'cart should be a 3D vector or array of 3D column vectors.'
    if len(cart.shape) == 2:
        return fkt.translations.kcart2latlon(cart)
    elif len(cart.shape) == 1:
        return fkt.translations.kcart2latlon(numpy.reshape(cart, (3, 1)))[:, 0]
def latlon2cart(latlon: numpy.ndarray) -> numpy.ndarray:
    """
    Latitude and longitude to Cartesian coordinates.

    Parameters:
    -----------
    `latlon` : numpy.ndarray, shape (2,), (2, n)

    2D Vector or array of column 2D vectors of latitude and longitude pairs.

    Vector structure: [lat, lon], where

    lat in [-pi/2, pi/2],

    lon in [-pi, pi]

    Returns:
    --------
    `cart` : numpy.ndarray, shape (3,), (3, n)

    3D vector or array of column 3D vectors of Cartesian coordinates.

    Vector structure: [x, y, z].

    Examples:
    ---------
    ```
    >> latlon = numpy.array([0, 0])

    >> cart = kiam.latlon2cart(latlon)

    [1. 0. 0.]
    ```
    """
    if latlon.shape[0] != 2 or len(latlon.shape) not in [1, 2]:
        raise 'latlon should be a 2D vector or array of 2D column vectors.'
    if len(latlon.shape) == 2:
        return fkt.translations.klatlon2cart(latlon)
    elif len(latlon.shape) == 1:
        return fkt.translations.klatlon2cart(numpy.reshape(latlon, (2, 1)))[:, 0]
def itrs2gcrs(xitrs: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vector from ITRS c/s to GCRS c/s.

    Parameters:
    -----------
    `xitrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the ITRS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to column(s) of xitrs

    `grad_req` : bool

    Flag to calculate the derivatives of the GCRS vector wrt the ITRS vector

    Returns:
    --------
    `xgcrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the GCRS coordinate system

    `dxgcrs` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    3x3 or 6x6 matrix or 3x3xn or 6x6xn array of partial derivatives of xgcrs wrt xitrs (dxgcrs/dxitrs).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    jd = kiam.juliandate(2022, 11, 22, 0, 0, 0)

    xitrs = numpy.array([1, 0, 0])

    xgcrs = kiam.itrs2gcrs(xitrs, jd, False)

    xgcrs, dxgcrs = kiam.itrs2gcrs(xitrs, jd, True)
    ```
    """
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))

    if len(xitrs.shape) == 1 and jd.shape[0] != 1:
        raise 'If xitrs is a vector, then jd should be a single number.'
    elif len(xitrs.shape) == 2 and xitrs.shape[1] != jd.shape[0]:
        raise 'Number of columns in xitrs should equal number of elements in jd.'

    if xitrs.shape == (3,):
        out = fkt.translations.kitrs2gcrs(xitrs, jd)
        return _return_if_grad_req(out, grad_req)

    if xitrs.shape == (6,):
        xitrs = numpy.reshape(xitrs, (6, 1))
        fkt.translations.xitrs_mat = xitrs
        fkt.translations.jd_mat = jd
        fkt.translations.kitrs2gcrs_mat()
        out = (fkt.translations.xgcrs_mat[:, 0], fkt.translations.dxgcrs_mat[:, :, 0])
        return _return_if_grad_req(out, grad_req)

    if len(xitrs.shape) == 2 and xitrs.shape[0] in [3, 6]:
        fkt.translations.xitrs_mat = xitrs
        fkt.translations.jd_mat = jd
        fkt.translations.kitrs2gcrs_mat()
        out = (fkt.translations.xgcrs_mat, fkt.translations.dxgcrs_mat)
    else:
        raise 'xitrs should be a 3D or 6D vector or an array of column 3D or 6D vectors.'
    return _return_if_grad_req(out, grad_req)
def gcrs2itrs(xgcrs: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vector from GCRS c/s to ITRS c/s.

    Parameters:
    -----------
    `xgcrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the GCRS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to column(s) in xgcrs

    `grad_req` : bool

    Flag to calculate the derivatives of the ITRS vector wrt the GCRS vector

    Returns:
    --------
    `xitrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the ITRS coordinate system

    `dxitrs` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    3x3 or 6x6 matrix or 3x3xn or 6x6xn array of partial derivatives of xitrs wrt xgcrs (dxitrs/dxgcrs).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    jd = kiam.juliandate(2022, 11, 22, 0, 0, 0)

    xgcrs = numpy.array([1, 0, 0])

    xitrs = kiam.gcrs2itrs(xgcrs, jd, False)

    xitrs, dxitrs = kiam.gcrs2itrs(xgcrs, jd, True)
    ```
    """
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))

    if len(xgcrs.shape) == 1 and jd.shape[0] != 1:
        raise 'If xgcrs is a vector, then jd should be a single number.'
    elif len(xgcrs.shape) == 2 and xgcrs.shape[1] != jd.shape[0]:
        raise 'Number of columns in xgcrs should equal number of elements in jd.'

    if xgcrs.shape == (3,):
        out = fkt.translations.kgcrs2itrs(xgcrs, jd)
        return _return_if_grad_req(out, grad_req)

    if xgcrs.shape == (6,):
        xgcrs = numpy.reshape(xgcrs, (6, 1))
        fkt.translations.xgcrs_mat = xgcrs
        fkt.translations.jd_mat = jd
        fkt.translations.kgcrs2itrs_mat()
        out = (fkt.translations.xitrs_mat[:, 0], fkt.translations.dxitrs_mat[:, :, 0])
        return _return_if_grad_req(out, grad_req)

    if len(xgcrs.shape) == 2 and xgcrs.shape[0] in [3, 6]:
        fkt.translations.xgcrs_mat = xgcrs
        fkt.translations.jd_mat = jd
        fkt.translations.kgcrs2itrs_mat()
        out = (fkt.translations.xitrs_mat, fkt.translations.dxitrs_mat)
    else:
        raise 'xgcrs should be a 3D or 6D vector or an array of column 3D or 6D vectors.'
    return _return_if_grad_req(out, grad_req)
def scrs2pa(xscrs: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vector from SCRS c/s to PA c/s.

    Parameters:
    -----------
    `xscrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SCRS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to column(s) in xscrs

    `grad_req` : bool

    Flag to calculate the derivatives of the PA vector wrt the SCRS vector

    Returns:
    --------
    `xpa` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the PA coordinate system

    `dxpa` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    3x3 or 6x6 matrix or 3x3xn or 6x6xn array of partial derivatives of xpa wrt xscrs (dxpa/dxscrs).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    jd = kiam.juliandate(2022, 11, 22, 0, 0, 0)

    xscrs = numpy.array([1, 0, 0])

    xpa = kiam.scrs2pa(xscrs, jd, False)

    xpa, dxpa = kiam.scrs2pa(xscrs, jd, True)
    ```
    """
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))

    if len(xscrs.shape) == 1 and jd.shape[0] != 1:
        raise 'If xscrs is a vector, then jd should be a single number.'
    elif len(xscrs.shape) == 2 and xscrs.shape[1] != jd.shape[0]:
        raise 'Number of columns in xscrs should equal number of elements in jd.'

    if xscrs.shape == (3,):
        out = fkt.translations.kscrs2pa(xscrs, jd)
        return _return_if_grad_req(out, grad_req)

    if xscrs.shape == (6,):
        xscrs = numpy.reshape(xscrs, (6, 1))
        fkt.translations.xscrs_mat = xscrs
        fkt.translations.jd_mat = jd
        fkt.translations.kscrs2pa_mat()
        out = (fkt.translations.xpa_mat[:, 0], fkt.translations.dxpa_mat[:, :, 0])
        return _return_if_grad_req(out, grad_req)

    if len(xscrs.shape) == 2 and xscrs.shape[0] in [3, 6]:
        fkt.translations.xscrs_mat = xscrs
        fkt.translations.jd_mat = jd
        fkt.translations.kscrs2pa_mat()
        out = (fkt.translations.xpa_mat, fkt.translations.dxpa_mat)
    else:
        raise 'xscrs should be a 3D or 6D vector or an array of column 3D or 6D vectors.'
    return _return_if_grad_req(out, grad_req)
def scrs2mer(xscrs: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vectors from SCRS c/s to MER c/s.

    Parameters:
    -----------
    `xscrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SCRS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to vector or columns in xscrs

    `grad_req` : bool

    Flag to calculate the derivatives of the MER vector wrt the SCRS vector

    Returns:
    --------
    `xmer` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the MER coordinate system.

    `dxmer` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    Matrix or array of matrices of partial derivatives of xmer wrt xscrs (dxmer/dxscrs).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    xscrs = numpy.array([1, 0, 0])

    jd = kiam.juliandate(2022, 11, 22, 0, 0, 0)

    xmer = kiam.scrs2mer(xscrs, jd, False)

    dxmer = kiam.scrs2mer(xscrs, jd, True)
    ```
    """
    initial_xscrs_shape = xscrs.shape
    if len(initial_xscrs_shape) == 1:
        xscrs = numpy.reshape(xscrs, (xscrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    chunk = 10000
    dim = xscrs.shape[0]
    ncols = xscrs.shape[1]
    if dim != 3 and dim != 6:
        raise 'xscrs should be a 3D or 6D vector or 3xn or 6xn array of vectors.'
    if xscrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xscrs should equal number of elements in jd.'
    xmer = numpy.empty((dim, ncols))
    if grad_req:
        dxmer = numpy.empty((dim, dim, ncols))
        for i in range(int(numpy.ceil(ncols / chunk))):
            cols = numpy.arange(i * chunk, min((i + 1) * chunk, ncols))
            xmer[:, cols], dxmer[:, :, cols] = fkt.translations.kscrs2mer(xscrs[:, cols], jd[cols])
        if len(initial_xscrs_shape) == 1:
            return xmer[:, 0], dxmer[:, :, 0]
        else:
            return xmer, dxmer
    else:
        for i in range(int(numpy.ceil(ncols / chunk))):
            cols = numpy.arange(i * chunk, min((i + 1) * chunk, ncols))
            xmer[:, cols], _ = fkt.translations.kscrs2mer(xscrs[:, cols], jd[cols])
        if len(initial_xscrs_shape) == 1:
            return xmer[:, 0]
        else:
            return xmer
def mer2scrs(xmer: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vectors from MER c/s to SCRS c/s.

    Parameters:
    -----------
    `xmer` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the MER coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to vector or columns in xmer

    `grad_req` : bool

    Flag to calculate the derivatives of the SCRS vector wrt the MER vector

    Returns:
    --------
    `xscrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SCRS coordinate system

    `dxscrs` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    Array of matrices of partial derivatives of xscrs wrt xmer (dxscrs/dxmer).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    xmer = numpy.array([1, 0, 0])

    jd = kiam.juliandate(2022, 11, 22, 0, 0, 0)

    xscrs = kiam.mer2scrs(xmer, jd, False)

    dxscrs = kiam.mer2scrs(xmer, jd, True)
    ```
    """
    initial_xmer_shape = xmer.shape
    if len(initial_xmer_shape) == 1:
        xmer = numpy.reshape(xmer, (xmer.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    chunk = 10000
    dim = xmer.shape[0]
    ncols = xmer.shape[1]
    if dim != 3 and dim != 6:
        raise 'xmer should be a 3D or 6D vector or 3xn or 6xn array of vectors.'
    if xmer.shape[1] != jd.shape[0]:
        raise 'number of columns in xmer should equal number of elements in jd.'
    xscrs = numpy.empty((dim, ncols))
    if grad_req:
        dxscrs = numpy.empty((dim, dim, ncols))
        for i in range(int(numpy.ceil(ncols / chunk))):
            cols = numpy.arange(i * chunk, min((i + 1) * chunk, ncols))
            xscrs[:, cols], dxscrs[:, :, cols] = fkt.translations.kmer2scrs(xmer[:, cols], jd[cols])
        if len(initial_xmer_shape) == 1:
            return xscrs[:, 0], dxscrs[:, :, 0]
        else:
            return xscrs, dxscrs
    else:
        for i in range(int(numpy.ceil(ncols / chunk))):
            cols = numpy.arange(i * chunk, min((i + 1) * chunk, ncols))
            xscrs[:, cols], _ = fkt.translations.kmer2scrs(xmer[:, cols], jd[cols])
        if len(initial_xmer_shape) == 1:
            return xscrs[:, 0]
        else:
            return xscrs
def scrs2gcrs(xscrs: numpy.ndarray, jd: Union[float, numpy.ndarray], dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from SCRS c/s to GCRS c/s.

    Parameters:
    -----------
    `xscrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the SCRS coordinate system

    Vector structure: [x, y, z, vx, vy, vz]

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xscrs

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xgcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the GCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    Example 1 (6D -> 6D):

    ku = kiam.units('earth', 'moon')

    xscrs = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xgcrs = kiam.scrs2gcrs(xscrs, jd, ku['DistUnit'], ku['VelUnit'])

    Example 2 (6x1 -> 6x1):

    ku = kiam.units('earth', 'moon')

    xscrs = numpy.array([[1, 0, 0, 0, 1, 0]]).T

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xgcrs = kiam.scrs2gcrs(xscrs, jd, ku['DistUnit'], ku['VelUnit'])
    ```
    """
    initial_xscrs_shape = xscrs.shape
    if len(initial_xscrs_shape) == 1:
        xscrs = numpy.reshape(xscrs, (xscrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xscrs.shape[0]
    if dim != 6:
        raise 'xscrs should be a 6D vector or 6xn array of vectors.'
    if xscrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xscrs should equal number of elements in jd.'
    xgcrs = fkt.translations.kscrs2gcrs(xscrs, jd, dist_unit, vel_unit)
    if len(initial_xscrs_shape) == 1:
        return xgcrs[:, 0]
    else:
        return xgcrs
def gcrs2scrs(xgcrs: numpy.ndarray, jd: Union[float, numpy.ndarray], dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from GCRS c/s to SCRS c/s.

    Parameters:
    -----------
    `xgcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the GCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz]

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xgcrs

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xscrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the SCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    Example 1 (6D -> 6D):

    ku = kiam.units('earth', 'moon')

    xgcrs = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xscrs = kiam.gcrs2scrs(xgcrs, jd, ku['DistUnit'], ku['VelUnit'])

    Example 2 (6x1 -> 6x1)

    ku = kiam.units('earth', 'moon')

    xgcrs = numpy.array([[1, 0, 0, 0, 1, 0]]).T

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xscrs = kiam.gcrs2scrs(xgcrs, jd, ku['DistUnit'], ku['VelUnit'])
    ```
    """
    initial_xgcrs_shape = xgcrs.shape
    if len(initial_xgcrs_shape) == 1:
        xgcrs = numpy.reshape(xgcrs, (xgcrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xgcrs.shape[0]
    if dim != 6:
        raise 'xgcrs should be a 6D vector or 6xn array of vectors.'
    if xgcrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xgcrs should equal number of elements in jd.'
    xscrs = fkt.translations.kgcrs2scrs(xgcrs, jd, dist_unit, vel_unit)
    if len(initial_xgcrs_shape) == 1:
        return xscrs[:, 0]
    else:
        return xscrs
def hcrs2gcrs(xhcrs: numpy.ndarray, jd: Union[float, numpy.ndarray], dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from HCRS c/s to GCRS c/s.

    Parameters:
    -----------
    `xhcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the HCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xhcrs

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xgcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the GCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    Example 1 (6D -> 6D)

    ku = kiam.units('sun', 'earth')

    xhcrs = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xgcrs = kiam.hcrs2gcrs(xhcrs, jd, ku['DistUnit'], ku['VelUnit'])

    Example 2 (6x1 -> 6x1)

    ku = kiam.units('sun', 'earth')

    xhcrs = numpy.array([[1, 0, 0, 0, 1, 0]]).T

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xgcrs = kiam.hcrs2gcrs(xhcrs, jd, ku['DistUnit'], ku['VelUnit'])
    ```
    """
    initial_xhcrs_shape = xhcrs.shape
    if len(initial_xhcrs_shape) == 1:
        xhcrs = numpy.reshape(xhcrs, (xhcrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xhcrs.shape[0]
    if dim != 6:
        raise 'xhcrs should be a 6D vector or 6xn array of vectors.'
    if xhcrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xhcrs should equal number of elements in jd.'
    xgcrs = fkt.translations.khcrs2gcrs(xhcrs, jd, dist_unit, vel_unit)
    if len(initial_xhcrs_shape) == 1:
        return xgcrs[:, 0]
    else:
        return xgcrs
def gcrs2hcrs(xgcrs: numpy.ndarray, jd: Union[float, numpy.ndarray], dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from GCRS c/s to HCRS c/s.

    Parameters:
    -----------
    `xgcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the GCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xgcrs

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xhcrs` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the HCRS coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    Example 1 (6D -> 6D):

    ku = kiam.units('sun', 'earth')

    xgcrs = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xhcrs = kiam.gcrs2hcrs(xgcrs, jd, ku['DistUnit'], ku['VelUnit'])

    Example 2 (6x1 -> 6x1):

    ku = kiam.units('sun', 'earth')

    xgcrs = numpy.array([[1, 0, 0, 0, 1, 0]]).T

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xhcrs = kiam.gcrs2hcrs(xgcrs, jd, ku['DistUnit'], ku['VelUnit'])
    ```
    """
    initial_xgcrs_shape = xgcrs.shape
    if len(initial_xgcrs_shape) == 1:
        xgcrs = numpy.reshape(xgcrs, (xgcrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xgcrs.shape[0]
    if dim != 6:
        raise 'xgcrs should be a 6D vector or 6xn array of vectors.'
    if xgcrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xgcrs should equal number of elements in jd.'
    xhcrs = fkt.translations.kgcrs2hcrs(xgcrs, jd, dist_unit, vel_unit)
    if len(initial_xgcrs_shape) == 1:
        return xhcrs[:, 0]
    else:
        return xhcrs
def scrs2sors(xscrs: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vectors from SCRS c/s to SORS c/s.

    Parameters:
    -----------
    `xscrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SCRS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xscrs

    `grad_req` : bool

    Flag to calculate the derivatives of the SORS vector wrt the SCRS vector

    Returns:
    --------
    `xsors` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SORS coordinate system

    `dxsors` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    Matrix or array of matrices of partial derivatives of xsors wrt xscrs (dxsors/dxscrs).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    xscrs = numpy.array([1, 0, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xsors = kiam.scrs2sors(xscrs, jd, False)

    xsors, dxsors = kiam.scrs2sors(xscrs, jd, True)
    ```
    """
    initial_xscrs_shape = xscrs.shape
    if len(initial_xscrs_shape) == 1:
        xscrs = numpy.reshape(xscrs, (xscrs.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xscrs.shape[0]
    if dim != 3 and dim != 6:
        raise 'xscrs should be a 3D or 6D vector or 3xn or 6xn array of vectors.'
    if xscrs.shape[1] != jd.shape[0]:
        raise 'number of columns in xscrs should equal number of elements in jd.'
    xsors, dxsors = fkt.translations.kscrs2sors(xscrs, jd)
    if grad_req:
        if len(initial_xscrs_shape) == 1:
            return xsors[:, 0], dxsors[:, :, 0]
        else:
            return xsors, dxsors
    else:
        if len(initial_xscrs_shape) == 1:
            return xsors[:, 0]
        else:
            return xsors
def sors2scrs(xsors: numpy.ndarray, jd: Union[float, numpy.ndarray], grad_req: bool = False) -> Union[numpy.ndarray, tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Translate vectors from SORS c/s to SCRS c/s.

    Parameters:
    -----------
    `xsors` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SORS coordinate system

    `jd` : float, numpy.ndarray, shape (n,)

    Julian dates corresponding to columns in xsors

    `grad_req` : bool

    Flag to calculate the derivatives of the SCRS vector wrt the SORS vector

    Returns:
    --------
    `xscrs` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the SCRS coordinate system

    `dxscrs` : numpy.ndarray, shape (3,3), (6,6), (3,3,n), (6,6,n)

    Matrix or array of matrices of partial derivatives of xscrs wrt xsors (dxscrs/dxsors).

    Returns only if `grad_req = True`.

    Examples:
    ---------
    ```
    xsors = numpy.array([1, 0, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    xscrs = kiam.sors2scrs(xsors, jd, False)

    xscrs, dxscrs = kiam.sors2scrs(xsors, jd, True)
    ```
    """
    initial_xsors_shape = xsors.shape
    if len(initial_xsors_shape) == 1:
        xsors = numpy.reshape(xsors, (xsors.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xsors.shape[0]
    if dim != 3 and dim != 6:
        raise 'xsors should be a 3D or 6D vector or 3xn or 6xn array of vectors.'
    if xsors.shape[1] != jd.shape[0]:
        raise 'number of columns in xsors should equal number of elements in jd.'
    xscrs, dxscrs = fkt.translations.ksors2scrs(xsors, jd)
    if grad_req:
        if len(initial_xsors_shape) == 1:
            return xscrs[:, 0], dxscrs[:, :, 0]
        else:
            return xscrs, dxscrs
    else:
        if len(initial_xsors_shape) == 1:
            return xscrs[:, 0]
        else:
            return xscrs
def ine2rot(xine: numpy.ndarray, t: Union[float, numpy.ndarray], t0: Union[float, numpy.ndarray]) -> numpy.ndarray:
    """
    Translate phase vectors from INE c/s to ROT c/s.

    Parameters:
    -----------
    `xine` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the INE coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `t` : float, numpy.ndarray, shape (1,), (n,)

    Time(s) corresponding to column(s) of xine

    `t0` : float, numpy.ndarray, shape (1,), (n,)

    Time(s) of INE and ROT c/s coincidence for each column of xine.

    Returns:
    --------
    `xrot` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the ROT coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    xine = numpy.array([1, 0, 0, 0, 1, 0])

    t = 1.0

    t0 = 0.0

    xrot = kiam.ine2rot(xine, t, t0)
    ```
    """
    initial_xine_shape = xine.shape
    if len(initial_xine_shape) == 1:
        xine = numpy.reshape(xine, (xine.shape[0], 1))
    if type(t) == float or t.shape == ():
        t = numpy.reshape(t, (1,))
    if type(t0) == float or t0.shape == ():
        t0 = numpy.reshape(t0, (1,))
    dim = xine.shape[0]
    if dim != 6:
        raise 'xine should be a 6D vector or 6xn array of vectors.'
    if xine.shape[1] != t.shape[0]:
        raise 'number of columns in xine should equal number of elements in t.'
    if xine.shape[1] != t0.shape[0] and 1 != t0.shape[0]:
        raise 'number of elements in t0 should equal 1 or number of elements in xine.'
    xrot = fkt.translations.kine2rot(xine, t, t0)
    if len(initial_xine_shape) == 1:
        return xrot[:, 0]
    else:
        return xrot
def rot2ine(xrot: numpy.ndarray, t: Union[float, numpy.ndarray], t0: Union[float, numpy.ndarray]) -> numpy.ndarray:
    """
    Translate phase vectors from ROT c/s to INE c/s.

    Parameters:
    -----------
    `xrot` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the ROT coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `t` : float, numpy.ndarray, shape (n,)

    Time(s) corresponding to column(s) of xrot

    `t0` : float, numpy.ndarray, shape (1,), (n,)

    Time(s) of of INE and ROT c/s coincidence for each column of xrot.

    Returns:
    --------
    `xine` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the INE coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    xrot = numpy.array([1, 0, 0, 0, 0, 0])

    t = 1.0

    t0 = 0.0

    xine = kiam.rot2ine(xrot, t, t0)
    ```
    """
    initial_xrot_shape = xrot.shape
    if len(initial_xrot_shape) == 1:
        xrot = numpy.reshape(xrot, (xrot.shape[0], 1))
    if type(t) == float or t.shape == ():
        t = numpy.reshape(t, (1,))
    if type(t0) == float or t0.shape == ():
        t0 = numpy.reshape(t0, (1,))
    dim = xrot.shape[0]
    if dim != 6:
        raise 'xrot should be a 6D vector or 6xn array of vectors.'
    if xrot.shape[1] != t.shape[0]:
        raise 'number of columns in xrot should equal number of elements in t.'
    if xrot.shape[1] != t0.shape[0] and 1 != t0.shape[0]:
        raise 'number of elements in t0 should equal 1 or number of elements in xrot.'
    xine = fkt.translations.krot2ine(xrot, t, t0)
    if len(initial_xrot_shape) == 1:
        return xine[:, 0]
    else:
        return xine
def ine2rot_eph(xine: numpy.ndarray, jd: Union[float, numpy.ndarray], first_body: str, secondary_body: str, dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from INEEPH c/s to ROTEPH c/s.

    Parameters:
    -----------
    `xine` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the INEEPH coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to column(s) in xine

    `first_body` : str

    Name of the first primary body

    Options: 'sun', 'mercury', 'venus', 'earth', 'moon', 'mars',
    'jupiter', 'saturn', 'uranus', 'neptune'

    `secondary_body` : str

    Name of the secondary primary body

    Options: 'sun', 'mercury', 'venus', 'earth', 'moon', 'mars',
    'jupiter', 'saturn', 'uranus', 'neptune'

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xrot` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the ROTEPH coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    xine = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    ku = kiam.units('earth', 'moon')

    xrot = kiam.ine2rot_eph(xine, jd, 'earth', 'moon', ku['DistUnit'], ku['VelUnit'])
    ```
    """
    first_body = first_body.capitalize()
    secondary_body = secondary_body.capitalize()
    if first_body == secondary_body:
        raise 'Bodies should be different.'
    initial_xine_shape = xine.shape
    if len(initial_xine_shape) == 1:
        xine = numpy.reshape(xine, (xine.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xine.shape[0]
    if dim != 6:
        raise 'xine should be a 6D vector or 6xn array of vectors.'
    if xine.shape[1] != jd.shape[0]:
        raise 'number of columns in xine should equal number of elements in jd.'
    xrot = fkt.translations.kine2roteph(xine, jd, first_body, secondary_body, dist_unit, vel_unit)
    if len(initial_xine_shape) == 1:
        return xrot[:, 0]
    else:
        return xrot
def rot2ine_eph(xrot: numpy.ndarray, jd: Union[float, numpy.ndarray], first_body: str, secondary_body: str, dist_unit: float, vel_unit: float) -> numpy.ndarray:
    """
    Translate phase vectors from ROTEPH c/s to INEEPH c/s.

    Parameters:
    -----------
    `xrot` : numpy.ndarray, shape (6,), (6,n)

    6D vector array or 6D column phase vectors in the ROTEPH coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    `jd` : float, numpy.ndarray, shape (n,)

    Julian date(s) corresponding to column(s) in xrot

    `first_body` : str

    Name of the first primary body

    Options: 'sun', 'mercury', 'venus', 'earth', 'moon', 'mars',
    'jupiter', 'saturn', 'uranus', 'neptune'

    `secondary_body` : str

    Name of the secondary primary body

    Options: 'sun', 'mercury', 'venus', 'earth', 'moon', 'mars',
    'jupiter', 'saturn', 'uranus', 'neptune'

    `dist_unit` : float

    The unit of distance in km

    `vel_unit` : float

    The unit of velocity in km/s

    Returns:
    --------
    `xine` : numpy.ndarray, shape (6,), (6,n)

    6D vector or array of 6D column phase vectors in the INEEPH coordinate system.

    Vector structure: [x, y, z, vx, vy, vz].

    Examples:
    ---------
    ```
    xrot = numpy.array([1, 0, 0, 0, 1, 0])

    jd = kiam.juliandate(2022, 12, 6, 0, 0, 0)

    ku = kiam.units('earth', 'moon')

    xine = kiam.ine2rot_eph(xrot, jd, 'earth', 'moon', ku['DistUnit'], ku['VelUnit'])
    ```
    """
    initial_xrot_shape = xrot.shape
    if len(initial_xrot_shape) == 1:
        xrot = numpy.reshape(xrot, (xrot.shape[0], 1))
    if type(jd) == float or jd.shape == ():
        jd = numpy.reshape(jd, (1,))
    dim = xrot.shape[0]
    if dim != 6:
        raise 'xrot should be a 6D vector or 6xn array of vectors.'
    if xrot.shape[1] != jd.shape[0]:
        raise 'number of columns in xrot should equal number of elements in jd.'
    xine = fkt.translations.krot2ineeph(xrot, jd, first_body, secondary_body, dist_unit, vel_unit)
    if len(initial_xrot_shape) == 1:
        return xine[:, 0]
    else:
        return xine
def mer2lvlh(xmer: numpy.ndarray, lat: float, lon: float) -> numpy.ndarray:
    """
    Translate phase vector(s) from MER c/s to LVLH c/s.

    Parameters:
    -----------
    `xmer` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the MER coordinate system

    `lat` : float

    Latitude that specifies the LVLH c/s in radians

    `lon` : float

    Longitude that specifies the LVLH c/s in radians

    Returns:
    --------
    `xlvlh` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the LVLH coordinate system

    Examples:
    ---------
    ```
    xmer = numpy.array([1, 2, 3])

    lat = 0.0

    lon = 1.0

    xlvlh = kiam.mer2lvlh(xmer, lat, lon)
    ```
    """

    if xmer.shape[0] not in [3, 6]:
        raise 'xmer should be a 3D or 6D vector or array of column 3D or 6D vectors'

    if xmer.shape == (3,):
        return fkt.translations.kmer2lvlh(xmer, lat, lon)

    if xmer.shape == (6,):
        xmer = numpy.reshape(xmer, (6, 1))
        fkt.translations.xmer_mat = xmer
        fkt.translations.kmer2lvlh_mat(lat, lon)
        return fkt.translations.xlvlh_mat[:, 0]

    if len(xmer.shape) == 2:
        fkt.translations.xmer_mat = xmer
        fkt.translations.kmer2lvlh_mat(lat, lon)
        return fkt.translations.xlvlh_mat
def lvlh2mer(xlvlh: numpy.ndarray, lat: float, lon: float) -> numpy.ndarray:
    """
    Translate phase vector(s) from LVLH c/s to MER c/s.

    Parameters:
    -----------
    `xlvlh` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the LVLH coordinate system

    `lat` : float

    Latitude that specifies the LVLH c/s in radians

    `lon` : float

    Longitude that specifies the LVLH c/s in radians

    Returns:
    --------
    `xmer` : numpy.ndarray, shape (3,), (6,), (3,n), (6,n)

    3D vector, 6D vector or array of 3D or 6D column vectors in the MER coordinate system

    Examples:
    ---------
    ```
    xlvlh = numpy.array([1, 2, 3])

    lat = 0.0

    lon = 1.0

    xmer = kiam.lvlh2mer(xlvlh, lat, lon)
    ```
    """
    if xlvlh.shape[0] not in [3, 6]:
        raise 'xlvlh should be a 3D or 6D vector or array of column 3D or 6D vectors'

    if xlvlh.shape == (3,):
        return fkt.translations.klvlh2mer(xlvlh, lat, lon)

    if xlvlh.shape == (6,):
        xlvlh = numpy.reshape(xlvlh, (6, 1))
        fkt.translations.xlvlh_mat = xlvlh
        fkt.translations.klvlh2mer_mat(lat, lon)
        return fkt.translations.xmer_mat[:, 0]

    if len(xlvlh.shape) == 2:
        fkt.translations.xlvlh_mat = xlvlh
        fkt.translations.klvlh2mer_mat(lat, lon)
        return fkt.translations.xmer_mat

# Units and constants (documented with examples)
def units(*args: str) -> dict:
    """
    Get units of distance, velocity, and time.

    Parameters:
    -----------
    `*args`

    Name or a pair of names of a celestial bodies

    Options for a single argument: 'earth', 'moon', 'sun', 'mercury', 'venus',
    'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto'

    Options for two arguments: ('earth', 'moon'), ('sun', 'earth')

    Returns:
    --------
    `units_dict` : dict

    A dictionary containing the units of distance, velocity, and time.

    Examples:
    ---------
    ```
    un = kiam.units('earth')

    DU = un['DistUnit']  # Unit of distance for the earth system of units

    un = kiam.units('earth', 'moon')

    VU = un['VelUnit']  # Unit of velocity for the Earth-Moon system of units
    ```
    """
    units_info = {}
    if len(args) == 1:
        output = fkt.constantsandunits.kunits_onebody(args[0].lower())
        units_info['GM'] = output[0]
        units_info['DistUnit'] = output[1]
        units_info['VelUnit'] = output[2]
        units_info['TimeUnit'] = output[3]
        units_info['AccUnit'] = output[4]
    elif len(args) == 2:
        output = fkt.constantsandunits.kunits_twobody(args[0].lower(), args[1].lower())
        units_info['GM'] = output[0]
        units_info['mu'] = output[1]
        units_info['DistUnit'] = output[2]
        units_info['VelUnit'] = output[3]
        units_info['TimeUnit'] = output[4]
        units_info['AccUnit'] = output[5]
    else:
        raise Exception('Wrong number of arguments in units.')
    return units_info
def astro_const() -> tuple[dict, dict, dict, dict, dict]:
    """
    Get astronomical constants.

    Returns:
    --------
    `uni_const` : dict

    Universal constants containing the speed of light (SoL) in km/s,
    astronomical unit (AU) in km, constant of gravitation (G) in km^3/kg/s^2,
    standard acceleration due to gravity (g0) in m/s^2, degrees in 1 radian (RAD).

    `star` : dict

    Contains a dictionary that constants of the Sun: the gravitational parameter (GM) in km^3/s^2,
    the mean radius (MeanRadius) in km.

    `planet` : dict

    Contains constants of the planets (Mercury, Venus, Earth, Mars, Jupiter,
    Saturn, Uranus, Neptune). The keys of the dictionary are the names of the planets.
    Eack planet[planet_name] is also a dictionary that contains the
    gravitational parameter of the planet (GM) in km^3/s^2,
    the mean raidus (MeanRadius) in km, the equator radius (EquatorRadius) in km,
    the semi-major axis of the orbit around the Sun (SemimajorAxis) in km. For the Earth
    there are additionaly the obliquity of the ecliptic (Obliquity) in degrees
    and its time derivative (dObliquitydt) in arcsec/cy (cy = century years).

    `moon` : dict

    Contains constants of the moons (currently only of the Moon). The dictionary
    has a single key named Moon and moon['Moon'] is also a dictionary.
    That dictionary contains the gravitational parameter of the Moon (GM) in km^3/s^2,
    the mean raidus (MeanRadius) in km,
    the semi-major axis of the orbit around the Sun (SemimajorAxis) in km.

    `small_body` : dict

    Contains constants of the small celestial bodies (currently only of the Pluto).
    The dictionary has a single key named Pluto and small_body['Pluto'] is also a
    dictionary. That dictionary contains the
    gravitational parameter of the Pluto (GM) in km^3/s^2,
    the mean raidus (MeanRadius) in km, the equator radius (EquatorRadius) in km,
    the semi-major axis of the orbit around the Sun (SemimajorAxis) in km.

    Examples:
    ---------
    ```
    uni_const, star, planet, moon, small_body = astro_const()  # If you need all the dicts

    _, star, planet, _, _ = astro_const()  # If you need only star and planet dicts

    star['Sun']['MeanRadius']  # Mean radius of the Sun

    planet['Earth']['GM']  # Gravitational parameter of the Earth

    planet['Mars']['SemimajorAxis']  # Semi-major axis of the Mars's orbit.
    ```
    """

    uni_const = {}
    star = {'Sun': {}}
    planet = {'Mercury': {}, 'Venus': {}, 'Earth': {}, 'Mars': {},
              'Jupiter': {}, 'Saturn': {}, 'Uranus': {}, 'Neptune': {}}
    moon = {'Moon': {}}
    small_body = {'Pluto': {}}

    uni_const['SoL'] = fkt.constantsandunits.uniconst_sol
    uni_const['AU'] = fkt.constantsandunits.uniconst_au
    uni_const['G'] = fkt.constantsandunits.uniconst_g
    uni_const['g0'] = fkt.constantsandunits.uniconst_g0
    uni_const['RAD'] = fkt.constantsandunits.uniconst_rad

    star['Sun']['GM'] = fkt.constantsandunits.sun_gm
    star['Sun']['MeanRadius'] = fkt.constantsandunits.sun_meanradius

    planet['Earth']['OrbitsAround'] = fkt.constantsandunits.earth_orbitsaround
    planet['Earth']['GM'] = fkt.constantsandunits.earth_gm
    planet['Earth']['MeanRadius'] = fkt.constantsandunits.earth_meanradius
    planet['Earth']['EquatorRadius'] = fkt.constantsandunits.earth_equatorradius
    planet['Earth']['SemimajorAxis'] = fkt.constantsandunits.earth_semimajoraxis
    planet['Earth']['Obliquity'] = fkt.constantsandunits.earth_obliquity
    planet['Earth']['dObliquitydt'] = fkt.constantsandunits.earth_dobliquitydt

    planet['Mercury']['OrbitsAround'] = fkt.constantsandunits.mercury_orbitsaround
    planet['Mercury']['GM'] = fkt.constantsandunits.mercury_gm
    planet['Mercury']['MeanRadius'] = fkt.constantsandunits.mercury_meanradius
    planet['Mercury']['EquatorRadius'] = fkt.constantsandunits.mercury_equatorradius
    planet['Mercury']['SemimajorAxis'] = fkt.constantsandunits.mercury_semimajoraxis

    planet['Venus']['OrbitsAround'] = fkt.constantsandunits.venus_orbitsaround
    planet['Venus']['GM'] = fkt.constantsandunits.venus_gm
    planet['Venus']['MeanRadius'] = fkt.constantsandunits.venus_meanradius
    planet['Venus']['EquatorRadius'] = fkt.constantsandunits.venus_equatorradius
    planet['Venus']['SemimajorAxis'] = fkt.constantsandunits.venus_semimajoraxis

    planet['Mars']['OrbitsAround'] = fkt.constantsandunits.mars_orbitsaround
    planet['Mars']['GM'] = fkt.constantsandunits.mars_gm
    planet['Mars']['MeanRadius'] = fkt.constantsandunits.mars_meanradius
    planet['Mars']['EquatorRadius'] = fkt.constantsandunits.mars_equatorradius
    planet['Mars']['SemimajorAxis'] = fkt.constantsandunits.mars_semimajoraxis

    planet['Jupiter']['OrbitsAround'] = fkt.constantsandunits.jupiter_orbitsaround
    planet['Jupiter']['GM'] = fkt.constantsandunits.jupiter_gm
    planet['Jupiter']['MeanRadius'] = fkt.constantsandunits.jupiter_meanradius
    planet['Jupiter']['EquatorRadius'] = fkt.constantsandunits.jupiter_equatorradius
    planet['Jupiter']['SemimajorAxis'] = fkt.constantsandunits.jupiter_semimajoraxis

    planet['Saturn']['OrbitsAround'] = fkt.constantsandunits.saturn_orbitsaround
    planet['Saturn']['GM'] = fkt.constantsandunits.saturn_gm
    planet['Saturn']['MeanRadius'] = fkt.constantsandunits.saturn_meanradius
    planet['Saturn']['EquatorRadius'] = fkt.constantsandunits.saturn_equatorradius
    planet['Saturn']['SemimajorAxis'] = fkt.constantsandunits.saturn_semimajoraxis

    planet['Uranus']['OrbitsAround'] = fkt.constantsandunits.uranus_orbitsaround
    planet['Uranus']['GM'] = fkt.constantsandunits.uranus_gm
    planet['Uranus']['MeanRadius'] = fkt.constantsandunits.uranus_meanradius
    planet['Uranus']['EquatorRadius'] = fkt.constantsandunits.uranus_equatorradius
    planet['Uranus']['SemimajorAxis'] = fkt.constantsandunits.uranus_semimajoraxis

    planet['Neptune']['OrbitsAround'] = fkt.constantsandunits.neptune_orbitsaround
    planet['Neptune']['GM'] = fkt.constantsandunits.neptune_gm
    planet['Neptune']['MeanRadius'] = fkt.constantsandunits.neptune_meanradius
    planet['Neptune']['EquatorRadius'] = fkt.constantsandunits.neptune_equatorradius
    planet['Neptune']['SemimajorAxis'] = fkt.constantsandunits.neptune_semimajoraxis

    moon['Moon']['OrbitsAround'] = fkt.constantsandunits.moon_orbitsaround
    moon['Moon']['GM'] = fkt.constantsandunits.moon_gm
    moon['Moon']['MeanRadius'] = fkt.constantsandunits.moon_meanradius
    moon['Moon']['SemimajorAxis'] = fkt.constantsandunits.moon_semimajoraxis

    small_body['Pluto']['OrbitsAround'] = fkt.constantsandunits.pluto_orbitsaround
    small_body['Pluto']['GM'] = fkt.constantsandunits.pluto_gm
    small_body['Pluto']['MeanRadius'] = fkt.constantsandunits.pluto_meanradius
    small_body['Pluto']['EquatorRadius'] = fkt.constantsandunits.pluto_equatorradius
    small_body['Pluto']['SemimajorAxis'] = fkt.constantsandunits.pluto_semimajoraxis

    return uni_const, star, planet, moon, small_body

# Equations of motion (documented without examples)
def r2bp(t: float, s: numpy.ndarray) -> numpy.ndarray:
    """
    Right-hand side of the restricted two-body problem equations of motion.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,)

    Phase state vector containing position and velocity.

    Vector structure [x, y, z, vx, vy, vz].

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,)

    Gravity acceleration according to the two-body model.

    Vector structure [fx, fy, fz, fvx, fvy, fvz].
    """
    return fkt.equationsmodule.kr2bp(t, s)
def cr3bp_fb(t: float, s: numpy.ndarray, mu: float, stm_req: bool) -> numpy.ndarray:
    """
    Right-hand side of the circular restricted three-body problem equations of motion
    wrt the first primary body.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing position and velocity and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm_req = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm_req = True

    `mu` : float

    Mass parameter of the three-body system

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Gravitational acceleration according to the circular restricted
    three-body model of motion wrt the first primary body extended
    (if stm_req = True) by the derivative of the state-transition matrix.

    Vector structure:

    [fx, fy, fz, fvx, fvy, fvz] if stm_req = False

    [fx, fy, fz, fvx, fvy, fvz, fm11, fm21, fm31, ... ] if stm_req = True
    """
    fkt.equationsmodule.massparameter = mu
    fkt.equationsmodule.stm_required = stm_req
    return fkt.equationsmodule.cr3bp_fb(t, s)
def cr3bp_sb(t: float, s: numpy.ndarray, mu: float, stm_req: bool) -> numpy.ndarray:
    """
    Right-hand side of the circular restricted three-body problem equations of motion
    wrt the secondary primary body.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing position and velocity and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm_req = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm_req = True

    `mu` : float

    Mass parameter of the three-body system

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Gravitational acceleration according to the circular restricted
    three-body model of motion wrt the secondary primary body extended
    (if stm_req = True) by the derivative of the state-transition matrix.

    Vector structure:

    [fx, fy, fz, fvx, fvy, fvz] if stm_req = False

    [fx, fy, fz, fvx, fvy, fvz, fm11, fm21, fm31, ... ] if stm_req = True
    """
    fkt.equationsmodule.massparameter = mu
    fkt.equationsmodule.stm_required = stm_req
    return fkt.equationsmodule.cr3bp_sb(t, s)
def nbp_rv_earth(t: float, s: numpy.ndarray, stm_req: bool, sources: dict, data: dict, units_data: dict) -> numpy.ndarray:
    """
    Right-hand side of the n-body problem equations of motion wrt the Earth in terms of
    the position and velocity variables.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing position and velocity and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm_req = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm_req = True

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    `sources` : dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `data` : dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `units_data` : dict

    A dictionary that contains the units.

    The dictionary keys:

    'DistUnit' (The unit of distance in km)

    'VelUnit'  (The unit of velocity in km/s)

    'TimeUnit' (The unit of time in days)

    'AccUnit'  (The unit of acceleration in m/s^2)

    'RSun'     (The radius of the Sun in the units of distance)

    'REarth'   (The radius of the Earth in the units of distance)

    'RMoon'    (The radius of the Moon in the units of distance)

    The units dictionary can be created by the kiam.prepare_units_dict() function.

    The gravitational parameter in the specified units should be 1.0.

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Gravitational acceleration according to the specified n-body problem equations
    of motion extended (if stm_req = True) by the derivative of the
    state-transition matrix.

    Vector structure:

    [fx, fy, fz, fvx, fvy, fvz] if stm_req = False

    [fx, fy, fz, fvx, fvy, fvz, fm11, fm21, fm31, ... ] if stm_req = True
    """
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_rv_earth(t, s)
def nbp_rv_moon(t: float, s: numpy.ndarray, stm_req: bool, sources: dict, data: dict, units_data: dict) -> numpy.ndarray:
    """
    Right-hand side of the n-body problem equations of motion wrt the Moon in terms of
    the position and velocity variables.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing position and velocity and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm_req = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm_req = True

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    `sources` : dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `data` : dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `units_data` : dict

    A dictionary that contains the units.

    The dictionary keys:

    'DistUnit' (The unit of distance in km)

    'VelUnit'  (The unit of velocity in km/s)

    'TimeUnit' (The unit of time in days)

    'AccUnit'  (The unit of acceleration in m/s^2)

    'RSun'     (The radius of the Sun in the units of distance)

    'REarth'   (The radius of the Earth in the units of distance)

    'RMoon'    (The radius of the Moon in the units of distance)

    The units dictionary can be created by the kiam.prepare_units_dict() function.

    The gravitational parameter in the specified units should be 1.0.

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Gravitational acceleration according to the specified n-body problem equations
    of motion extended (if stm_req = True) by the derivative of the
    state-transition matrix.

    Vector structure:

    [fx, fy, fz, fvx, fvy, fvz] if stm_req = False

    [fx, fy, fz, fvx, fvy, fvz, fm11, fm21, fm31, ... ] if stm_req = True
    """
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_rv_moon(t, s)
def nbp_ee_earth(t: float, s: numpy.ndarray, stm_req: bool, sources: dict, data: dict, units_data: dict) -> numpy.ndarray:
    """
    Right-hand side of the n-body problem equations of motion wrt the Earth in terms of
    the equinoctial orbital elements.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing equinoctial elements and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [h, ex, ey, ix, iy, L] if stm_req = False,

    [h, ex, ey, ix, iy, L, m11, m21, m31, ...] if stm_req = True,

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    `sources` : dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `data` : dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `units_data` : dict

    A dictionary that contains the units.

    The dictionary keys:

    'DistUnit' (The unit of distance in km)

    'VelUnit'  (The unit of velocity in km/s)

    'TimeUnit' (The unit of time in days)

    'AccUnit'  (The unit of acceleration in m/s^2)

    'RSun'     (The radius of the Sun in the units of distance)

    'REarth'   (The radius of the Earth in the units of distance)

    'RMoon'    (The radius of the Moon in the units of distance)

    The units dictionary can be created by the kiam.prepare_units_dict() function.

    The gravitational parameter in the specified units should be 1.0.

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Phase state time dericatives according to the specified n-body problem equations
    of motion extended (if stm_req = True) by the derivative of the
    state-transition matrix.

    Vector structure:

    [fh, fex, fey, fix, fiy, fL] if stm_req = False

    [fh, fex, fey, fix, fiy, fL, fm11, fm21, fm31, ... ] if stm_req = True
    """
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_ee_earth(t, s)
def nbp_ee_moon(t: float, s: numpy.ndarray, stm_req: bool, sources: dict, data: dict, units_data: dict) -> numpy.ndarray:
    """
    Right-hand side of the n-body problem equations of motion wrt the Moon in terms of
    the equinoctial orbital elements.

    Parameters:
    -----------
    `t` : float

    Time

    `s` : numpy.ndarray, shape (6,), (42,)

    Phase state vector containing equinoctial elements and (if stm_req = True)
    vectorized state-transition matrix.

    Vector structure:

    [h, ex, ey, ix, iy, L] if stm_req = False,

    [h, ex, ey, ix, iy, L, m11, m21, m31, ...] if stm_req = True,

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination

    `stm_req` : bool

    Flag to calculate the derivative of the state-transition matrix

    `sources` : dict

    Dictionary that contains the perturbations that should be accounted.
    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `data` : dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `units_data` : dict

    A dictionary that contains the units.

    The dictionary keys:

    'DistUnit' (The unit of distance in km)

    'VelUnit'  (The unit of velocity in km/s)

    'TimeUnit' (The unit of time in days)

    'AccUnit'  (The unit of acceleration in m/s^2)

    'RSun'     (The radius of the Sun in the units of distance)

    'REarth'   (The radius of the Earth in the units of distance)

    'RMoon'    (The radius of the Moon in the units of distance)

    The units dictionary can be created by the kiam.prepare_units_dict() function.

    The gravitational parameter in the specified units should be 1.0.

    Returns:
    --------
    `f` : numpy.ndarray, shape (6,), (42,)

    Phase state time dericatives according to the specified n-body problem equations
    of motion extended (if stm_req = True) by the derivative of the
    state-transition matrix.

    Vector structure:

    [fh, fex, fey, fix, fiy, fL] if stm_req = False

    [fh, fex, fey, fix, fiy, fL, fm11, fm21, fm31, ... ] if stm_req = True
    """
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_ee_moon(t, s)

# Propagation routines (documented without examples)
def propagate_nbp(central_body: str, tspan: numpy.ndarray, x0: numpy.ndarray, sources_dict: dict, dat_dict: dict, stm: bool, variables: str) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Propagate trajectory in the n-body model of motion.

    Parameters:
    -----------
    `central_body`: str

    Name of the central body

    `tspan`: numpy.ndarray, shape (n,)

    Time nodes at which the solution is required

    `x0`: numpy.ndarray, shape (6,), (42,)

    Initial state containing:

    position and velocity (if variables = 'rv', stm = False),

    position and velocoty extended by vectorized state-transition matrix (if variables = 'rv_stm', stm = True),

    equinoctial orbital elements (if variables = 'ee', stm = False),

    equinoctial orbital elements extended by vectorized state-transition matrix (if variables = 'ee_stm', stm = True),

    Vector structure:

    [x, y, z, vx, vy, vz] if variables = 'rv' and stm = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if variables = 'rv_stm' and stm = True

    [h, ex, ey, ix, iy, L] if variables = 'ee' and stm = False

    [h, ex, ey, ix, iy, L, m11, m21, m31, ...] if variables = 'ee_stm' and stm = True

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination

    `sources_dict`: dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `dat_dict`: dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `stm`: bool

    Flag to calculate the derivative of the state-transition matrix

    `variables`: str

    Type of variables used to propagate the trajectory.

    If stm = False, then variables should be 'rv' or 'ee'.

    If stm = True, then variables should be 'rv_stm' or 'ee_stm'.


    Returns:
    --------
    `t` : numpy.ndarray, shape(n,)

    Times (nodes) in tspan at which the solution is obtained

    `y` : numpy.ndarray, shape(6, n), shape(42, n)

    Array of column trajectory phase states extended (if stm = True) by
    vectorized state-transition matrices.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm = False and variables = 'rv'

    [x, y, z, vx, vy, vz, m11, m21, m31, ... ] if stm = True and variables = 'rv_stm'

    [h, ex, ey, ix, iy, L] if stm = False and variables = 'ee'

    [h, ex, ey, ix, iy, L, m11, m21, m31, ... ] if stm = True and variables = 'ee_stm'

    h = sqrt(p/mu),

    ex = e*cos(Omega+omega),

    ey = e*sin(Omega+omega),

    ix = tan(i/2)*cos(Omega),

    iy = tan(i/2)*sin(Omega),

    L = theta + omega + Omega,

    where

    mu - gravitational parameter,

    p - semi-latus rectum (focal parameter),

    e - eccentricity,

    Omega - right ascension of the ascending node,

    omega - argument of pericenter,

    i - inclination
    """
    tspan, x0 = to_float(tspan, x0)
    neq = 42 if stm else 6
    if variables == 'rv_stm':
        variables = 'rv'
    elif variables == 'ee_stm':
        variables = 'ee'
    sources_vec = _sources_dict_to_vec(sources_dict)
    dat_vec = _dat_dict_to_vec(dat_dict)
    t, y = fkt.propagationmodule.propagate_nbp(central_body.lower(), tspan, x0, sources_vec, dat_vec,
                                               stm, variables, neq)
    return t, y
def propagate_r2bp(tspan: numpy.ndarray, x0: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Propagate trajectory in the two-body model of motion.

    Parameters:
    -----------
    `tspan`: numpy.ndarray, shape (n,)

    Time nodes at which the solution is required

    `x0`: numpy.ndarray, shape (6,)

    Initial state containing position and velocity.

    Vector structure: [x, y, z, vx, vy, vz]

    Returns:
    --------
    `t` : numpy.ndarray, shape(n,)

    Times (nodes) in tspan at which the solution is obtained

    `y` : numpy.ndarray, shape(6, n)

    Array of column trajectory phase states.

    Vector structure: [x, y, z, vx, vy, vz].
    """
    tspan, x0 = to_float(tspan, x0)
    t, y = fkt.propagationmodule.propagate_r2bp(tspan, x0)
    return t, y
def propagate_cr3bp(central_body: str, tspan: numpy.ndarray, x0: numpy.ndarray, mu: float, stm: bool) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Propagate trajectory in the circular restricted three-body model of motion.

    Parameters:
    -----------
    `central_body` : str

    First or secondary primary body as the central body

    `tspan`: numpy.ndarray, shape (n,)

    Time nodes at which the solution is required

    `x0` : numpy.ndarray, shape (6,), (42,)

    Initial state containing:

    position and velocity (if stm = False),

    position and velocoty extended by vectorized state-transition matrix (if stm = True),

    Vectory structure:

    [x, y, z, vx, vy, vz] if stm = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm = True

    `mu` : float

    Mass parameter of the three-body system

    `stm` : bool

    Flag to calculate the derivative of the state-transition matrix

    Returns:
    --------
    `t` : numpy.ndarray, shape(n,)

    Times (nodes) in tspan at which the solution is obtained

    `y` : numpy.ndarray, shape(6, n)

    Array of column trajectory phase states.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm = True
    """
    tspan, x0, mu = to_float(tspan, x0, mu)
    neq = 42 if stm else 6
    t, y = fkt.propagationmodule.propagate_cr3bp(central_body.lower(), tspan, x0, mu, stm, neq)
    return t, y
def propagate_br4bp(central_body: str, tspan: numpy.ndarray, x0: numpy.ndarray, mu: float, gm4b, a4b: float, theta0: float, stm: bool) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Propagate trajectory in the bi-circular restricted four-body model of motion.

    Parameters:
    -----------
    `central_body` : str

    First or secondary primary body as the central body

    `tspan`: numpy.ndarray, shape (n,)

    Time nodes at which the solution is required

    `x0` : numpy.ndarray, shape (6,), (42,)

    Initial state containing:

    position and velocity (if stm = False),

    position and velocoty extended by vectorized state-transition matrix (if stm = True),

    Vectory structure:

    [x, y, z, vx, vy, vz] if stm = False

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm = True

    `mu` : float

    Mass parameter of the three-body system

    `gm4b` : float

    Scaled gravitational parameter of the fourth (perturbing) body

    `a4b` : float

    Distance from the center of mass of the primary bodies to the fourth body
    in units where the distance between the primaries equals 1.

    `theta0` : float

    Initial value of the synodic phase - the angle between the direction to
    the fourth body from the center of mass of the primaries and the line
    connecting the primaties.

    `stm` : bool

    Flag to calculate the derivative of the state-transition matrix.

    Returns:
    --------
    `t` : numpy.ndarray, shape(n,)

    Times (nodes) in tspan at which the solution is obtained.

    `y` : numpy.ndarray, shape(6, n)

    Array of column trajectory phase states.

    Vector structure:

    [x, y, z, vx, vy, vz] if stm = False.

    [x, y, z, vx, vy, vz, m11, m21, m31, ...] if stm = True.
    """
    tspan, x0, mu, gm4b, a4b, theta0 = to_float(tspan, x0, mu, gm4b, a4b, theta0)
    neq = 42 if stm else 6
    t, y = fkt.propagationmodule.propagate_br4bp(central_body.lower(), tspan, x0, mu, gm4b, a4b, theta0, stm, neq)
    return t, y
def prepare_sources_dict() -> dict:
    """
    Auxilary function that returns a dictionary of perturbations.

    Returns:
    --------
    `sources` : dict

    The dictionary keys:

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
    """
    sources = {'sun': False, 'mercury': False, 'venus': False, 'earth': False,
               'moon': False, 'mars': False, 'jupiter': False, 'saturn': False,
               'uranus': False, 'neptune': False, 'srp': False, 'cmplxmoon': False,
               'atm': False, 'j2': False}
    return sources
def prepare_units_dict(units_name: str) -> dict:
    """
    Auxilary function that returns a dictionary of units.

    Parameters:
    -----------
    units_name : str

    A name of the units that should be used.

    Options:

    'dim': a dictionary with unity values will be returned

    'earth': a dictionary of the earth units will be returned

    'moon': a dictionary of the moon units will be returned

    Returns:
    --------
    units_dict : dict

    A dictionary that containes the following keys:

    'DistUnit': the unit of distance, km

    'VelUnit': the unit of velocity, km/s

    'TimeUnit': the unit of time, days

    'AccUnit': the units of acceleration, m/s^2

    'RSun': the mean radius of the Sun in DistUnit units

    'REarth': the mean radius of the Earth in DistUnit units

    'RMoon': the mean radius of the Moon in DistUnit units

    """

    _, star, planet, moon, _ = astro_const()
    RSun_km = star['Sun']['MeanRadius']
    REarth_km = planet['Earth']['MeanRadius']
    RMoon_km = moon['Moon']['MeanRadius']

    units_data = {'DistUnit': 1.0,
                  'VelUnit': 1.0,
                  'TimeUnit': 1.0,
                  'AccUnit': 1.0,
                  'RSun': RSun_km/1.0,
                  'REarth': REarth_km/1.0,
                  'RMoon': RMoon_km/1.0}

    if units_name == 'dim':
        return units_data

    if units_name in ['earth', 'moon']:
        ku = units(units_name)
        units_data['DistUnit'] = ku['DistUnit']
        units_data['VelUnit'] = ku['VelUnit']
        units_data['TimeUnit'] = ku['TimeUnit']
        units_data['AccUnit'] = ku['AccUnit']
        units_data['RSun'] = RSun_km / ku['DistUnit']
        units_data['REarth'] = REarth_km / ku['DistUnit']
        units_data['RMoon'] = RMoon_km / ku['DistUnit']
        return units_data

    raise 'Unknown units_name.'

# Visibility routines (documented without examples)
def is_visible(r_sat: numpy.ndarray, lat_deg: Union[int, float, numpy.ndarray],
               long_deg: Union[int, float, numpy.ndarray], body_radius: float,
               threshold_deg: Union[int, float, numpy.ndarray])\
        -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """
    Get visibility statuses (0 or 1) of a vector from a point on a sphere surface.

    Parameters:
    -----------
    `r_sat` : numpy.ndarray, shape (3,), (3,n)

    Radius-vector(s) around a sphere surface

    `lat_deg` : int, float, numpy.ndarray, shape (m,)

    Latitude of a point(s) on a surface in degrees

    `long_deg` : int, float, numpy.ndarray, shape (m,)

    Longitude of a point(s) on a surface in degrees

    `body_radius` : float

    Body radius

    `threshold_deg` : int, float, numpy.ndarray, shape (n,)

    Minimum angle below which the vector is not visible.

    Returns:
    --------
    `status` : numpy.ndarray, shape (n,m)

    Visibility statuses of the r_sat vectors from lat_deg/long_deg points.

    n - number of vectors in r_sat

    m - number of points on the surface

    `elev_deg` : numpy.ndarray, shape (n,m)

    Elevation angles in degrees

    n - number of vectors in r_sat

    m - number of points on the surface

    `azim_deg` : numpy.ndarray, shape (n,m)

    Azimuth angles in degrees

    n - number of vectors in r_sat

    m - number of points on the surface
    """

    if len(r_sat.shape) == 1:
        if r_sat.shape[0] != 3:
            raise 'r_sat as a vector should have 3 components.'
        r_sat = numpy.reshape(r_sat, (3, 1), order='F')
    if r_sat.shape[0] != 3 or len(r_sat.shape) != 2:
        raise 'r_sat as a matrix should have 3 rows and N columns.'
    r_sat = r_sat.copy().T / body_radius  # transpose for Fortran module to get n x 3

    if isinstance(lat_deg, (float, int)):
        lat_deg = numpy.array([lat_deg])
    if isinstance(long_deg, (float, int)):
        long_deg = numpy.array([long_deg])
    if len(lat_deg.shape) != 1:
        raise 'lat_deg should be a scalar or a vector.'
    if len(long_deg.shape) != 1:
        raise 'long_deg should be a scalar or a vector.'
    if lat_deg.shape[0] != long_deg.shape[0]:
        raise 'lat_deg and long_deg should have the same size.'
    lat_long = numpy.reshape(numpy.concatenate((lat_deg/180*numpy.pi, long_deg/180*numpy.pi), axis=0), (2, -1))

    threshold = threshold_deg / 180 * numpy.pi
    if isinstance(threshold, (float, int)):
        threshold = numpy.full((r_sat.shape[0],), threshold)
    if len(threshold.shape) != 1:
        raise 'threshold_deg should be a scalar or a vector'
    if threshold.shape[0] != r_sat.shape[0]:
        raise 'threshold_deg should have r_sat.shape[1] number of elements'

    fkt.visibilitymodule.r_sat = r_sat
    fkt.visibilitymodule.lat_long = lat_long
    fkt.visibilitymodule.threshold = threshold

    fkt.visibilitymodule.isvisible(r_sat.shape[0], lat_long.shape[1])

    status = fkt.visibilitymodule.status.copy()
    elev = fkt.visibilitymodule.elev.copy()
    azim = fkt.visibilitymodule.azim.copy()

    status[status == -1] = 1
    elev_deg = elev/numpy.pi*180
    azim_deg = azim/numpy.pi*180

    return status, elev_deg, azim_deg

# Save and load routines (documented without examples)
def save(variable: str, filename: str) -> None:
    """
    Saves a variable into a specified file.

    Parameters:
    -----------
    `variable` : str

    Variable to be saved.

    For limitations on variables see the pickle package
    https://docs.python.org/3/library/pickle.html

    `filename` : str

    A path to the file.
    """
    with open(filename, 'wb') as f:
        pickle.dump(variable, f)
def load(filename: str) -> Any:
    """
    Loads a variable from a specified file.

    Parameters:
    -----------
    `filename` : str

    A path to the file.

    Returns:
    --------
    `var` : Any

    A variable contained in the file.
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)

# General astrodynamics (documented without examples)
def get_period_hours(altitude_km: float, body: str) -> float:
    """
    Calculate the circular orbit period with a given altitude
    around a specified celestial body.

    Parameters:
    -----------
    `altitude_km` : float

    The altitude above the surface of the body in km.

    `body` : str

    The name of the celesial body.

    Returns:
    --------
    `period` : float

    The circular orbit period in hours.
    """
    ku = units(body)
    return 2*numpy.pi*numpy.sqrt((ku['DistUnit']+altitude_km)**3/ku['GM'])/3600.0
def get_altitude_km(period_hours: float, body: str) -> float:
    """
    Calculate the altitude of a circular orbit with a given period
    around a specified celestial body.

    Parameters:
    -----------
    `period` : float

    The circular orbit period in hours.

    `body`: str

    The name of the celesial body.

    Returns:
    --------
    `altitude_km` : float

    The altitude above the surface of the body in km.
    """
    ku = units(body)
    period = period_hours*3600
    return (period**2/(4*numpy.pi**2)*ku['GM'])**(1/3) - ku['DistUnit']
def get_circular_velocity_km_s(altitude_km: float, body: str) -> float:
    """
    Calculate the circular velocity at a given altitude
    around a specified celestial body.

    Parameters:
    -----------
    `altitude_km` : float

    The altitude above the surface of the body in km.

    `body` : str

    The name of the celesial body.

    Returns:
    --------
    `velocity` : float

    The circular velocity at the given altitude.
    """
    ku = units(body)
    return numpy.sqrt(ku['GM']/(ku['DistUnit'] + altitude_km))
def get_dv_hohmann(r1_nondim: float, r2_nondim: float) -> float:
    """
    Calculate delta-v in a Hohmann transfer.

    Parameters:
    -----------
    `r1_nondim` : float

    Nondimensional distance to the center of mass of the central body at the start.

    `r2_nondim` : float

    Nondimensional distance to the center of mass of the central body at the end.

    Returns:
    --------
    `dv` : float

    Nondimensional delta-v in the Hohmann transfer connecting r1_nondim and r2_nondim.
    It is assumed that the gravitational parameter equals 1.0.
    """
    dv1 = numpy.sqrt(1.0 / r1_nondim) * (numpy.sqrt(2 * r2_nondim / (r1_nondim + r2_nondim)) - 1)
    dv2 = numpy.sqrt(1.0 / r2_nondim) * (numpy.sqrt(2 * r1_nondim / (r1_nondim + r2_nondim)) - 1)
    dv_nondim = dv1 + dv2
    return dv_nondim
def get_tof_hohmann(r1_nondim: float, r2_nondim: float) -> float:
    """
    Calculate the time of flight in a Hohmann transfer.

    Parameters:
    -----------
    `r1_nondim` : float

    Nondimensional distance to the center of mass of the central body at the start.

    `r2_nondim` : float

    Nondimensional distance to the center of mass of the central body at the end.

    Returns:
    --------
    `tof` : float

    Nondimensional time of flight in the Hohmann transfer connecting
    r1_nondim and r2_nondim. It is assumed that the gravitational parameter
    equals 1.0.
    """
    a = (r1_nondim + r2_nondim) / 2
    tof_nondim = numpy.pi * (a ** 1.5)
    return tof_nondim

# Trofimov-Shirobokov model (documented without examples)
def get_order(altitude_thousands_km: float, approx_level: str = 'soft') -> int:
    """
    The minimum order and degree of the complex lunar gravitational field
    at a given altitude according to the Trofimov--Shirobokov model.

    Parameters:
    -----------
    `altitude_thousands_km` : float

    The altitude above the lunar surface in km.

    `approx_level` : str

    The level of approximation, can be 'soft' or 'hard'.

    Returns:
    --------
    `order` : int

    The order and degree of the complex lunar gravitational field.
    """
    if approx_level == 'soft':
        return numpy.floor((25.0 / altitude_thousands_km)**0.8)+1
    elif approx_level == 'hard':
        return numpy.floor((40.0 / altitude_thousands_km)**0.8)+1
    else:
        raise Exception('Unknown approx_level.')

# Auxilary protected methods (documented without examples)
def _set_nbp_parameters(stm_req: bool, sources: dict, data: dict, units_data: dict) -> None:
    """
    FOR THE TOOLBOX DEVELOPERS ONLY.
    Safe setting of the parameters for getting the right-hand side
    of equations of motion in ephemeris model. Writes to fkt.

    Parameters:
    -----------
    `stm_req`: bool

    Flag to calculate the derivative of the state-transition matrix

    `sources`: dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    `data`: dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    `units_data` : dict

    A dictionary that contains the units.

    The dictionary keys:

    'DistUnit' (The unit of distance in km)

    'VelUnit'  (The unit of velocity in km/s)

    'TimeUnit' (The unit of time in days)

    'AccUnit'  (The unit of acceleration in m/s^2)

    'RSun'     (The radius of the Sun in the units of distance)

    'REarth'   (The radius of the Earth in the units of distance)

    'RMoon'    (The radius of the Moon in the units of distance)
    """

    fkt.equationsmodule.stm_required = stm_req

    fkt.equationsmodule.atm = sources['atm']
    fkt.equationsmodule.j2 = sources['j2']
    fkt.equationsmodule.srp = sources['srp']
    fkt.equationsmodule.sun = sources['sun']
    fkt.equationsmodule.mercury = sources['mercury']
    fkt.equationsmodule.venus = sources['venus']
    fkt.equationsmodule.earth = sources['earth']
    fkt.equationsmodule.mars = sources['mars']
    fkt.equationsmodule.jupiter = sources['jupiter']
    fkt.equationsmodule.saturn = sources['saturn']
    fkt.equationsmodule.uranus = sources['uranus']
    fkt.equationsmodule.neptune = sources['neptune']
    fkt.equationsmodule.moon = sources['moon']
    fkt.equationsmodule.cmplxmoon = sources['cmplxmoon']

    fkt.equationsmodule.jd_zero = data['jd_zero']
    fkt.equationsmodule.order = data['order']
    fkt.equationsmodule.area = data['area']
    fkt.equationsmodule.mass = data['mass']

    fkt.equationsmodule.distunit = units_data['DistUnit']
    fkt.equationsmodule.velunit = units_data['VelUnit']
    fkt.equationsmodule.timeunit = units_data['TimeUnit']
    fkt.equationsmodule.accunit = units_data['AccUnit']
    fkt.equationsmodule.rsun = units_data['RSun']
    fkt.equationsmodule.rearth = units_data['REarth']
    fkt.equationsmodule.rmoon = units_data['RMoon']
def _return_if_grad_req(out: tuple[numpy.ndarray, numpy.ndarray], grad_req: bool) -> Union[tuple[numpy.ndarray, numpy.ndarray], numpy.ndarray]:
    """
    FOR THE TOOLBOX DEVELOPERS ONLY.
    Returns a vector or a (vector, gradient) pair as specified.

    Parameters:
    -----------
    `out` : tuple[numpy.ndarray, numpy.ndarray]

    A tuple containing the vector and a matrix.

    The matrix can be a gradient of the vector or a dump matrix.

    It is assumed that when grad_req = True, the matrix is a gradient (not dump).

    `grad_req` : bool

    Flag to calculate the gradient

    Returns:
    --------
    `out` : tuple[numpy.ndarray, numpy.ndarray]

    The (vector, gradient) pair if grad_req = True.

    `vec` : numpy.ndarray

    Only the vector = out[0] if grad_req = True.

    In this case the matrix out[1] (possibly dump) is ignored.
    """
    if grad_req:
        return out
    else:
        return out[0]
def _sources_dict_to_vec(sources_dict: dict) -> numpy.ndarray:
    """
    FOR THE TOOLBOX DEVELOPERS ONLY.
    Converts the sources dictionary to a vector that is understandable by
    the Fortran module that propagates the trajectory.

    Parameters:
    -----------
    `sources_dict` : dict

    Dictionary that contains the perturbations that should be accounted.

    The dictionary keys:

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

    If sources[key] = True, the corresponding perturbation will be accounted.

    If sources[key] = False, the corresponding perturbation will not be accounted.

    The sources dictionary with all False values can be created by
    the kiam.prepare_sources_dict() function.

    Returns:
    --------
    `sources_vec` : numpy.ndarray, shape(14,)

    Vectorized version of sources_dict.

    """
    sources_vec = numpy.zeros((len(sources_dict),))
    sources_vec[0] = int(sources_dict['sun'])
    sources_vec[1] = int(sources_dict['mercury'])
    sources_vec[2] = int(sources_dict['venus'])
    sources_vec[3] = int(sources_dict['earth'])
    sources_vec[4] = int(sources_dict['moon'])
    sources_vec[5] = int(sources_dict['mars'])
    sources_vec[6] = int(sources_dict['jupiter'])
    sources_vec[7] = int(sources_dict['saturn'])
    sources_vec[8] = int(sources_dict['uranus'])
    sources_vec[9] = int(sources_dict['neptune'])
    sources_vec[10] = int(sources_dict['srp'])
    sources_vec[11] = int(sources_dict['cmplxmoon'])
    sources_vec[12] = int(sources_dict['atm'])
    sources_vec[13] = int(sources_dict['j2'])
    return sources_vec
def _dat_dict_to_vec(dat_dict: dict) -> numpy.ndarray:
    """
    FOR THE TOOLBOX DEVELOPERS ONLY.
    Converts the data dictionary to a vector that is understandable by
    the Fortran module that propagates the trajectory.

    Parameters:
    -----------
    `dat_dict` : dict

    A dictionary that contains auxilary data.

    The dictionary keys:

    'jd_zero' (Julian date that corresponds to t = 0)

    'order'   (Order of the lunar complex gravitational field)

    'area'    (Area of the spacecraft to account in atmospheric drag and SRP, m^2)

    'mass'    (Mass of the spacecraft to account in atmospheric drag and SRP, kg)

    The data should be submitted even if the corresponding perturbations
    are not accounted.

    Returns:
    --------
    `dat_vec` : numpy.ndarray, shape(4,)

    Vectorized version of dat_dict.
    """
    dat_vec = numpy.zeros((4,))
    dat_vec[0] = dat_dict['jd_zero']
    dat_vec[1] = dat_dict['area']
    dat_vec[2] = dat_dict['mass']
    dat_vec[3] = dat_dict['order']
    return dat_vec
