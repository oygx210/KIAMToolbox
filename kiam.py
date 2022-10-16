import FKIAMToolbox as fkt
import jdcal
import datetime
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pickle

mpl.rcParams['figure.dpi'] = 150
mpl.use('Qt5Agg')

# General mathematics and tools.
def dotinvAB(A, B):
    return np.linalg.solve(A, B)
def dotAinvB(A, B):
    C = np.linalg.solve(B.T, A.T)
    return C.T
def eyevec(n):
    return np.reshape(np.eye(n), (n**2,))
def to_float(*args):
    args_float = tuple(np.array(arg, dtype='float64') for arg in args)
    return args_float
def sind(x):
    return np.sin(x/180*np.pi)
def cosd(x):
    return np.cos(x/180*np.pi)
def tand(x):
    return np.tan(x/180*np.pi)
def cotand(x):
    return 1/np.tan(x/180*np.pi)

# Plotting functions
def plot(x, y, ax=None, style='-', xlabel='', ylabel='', linewidth=1.0, show=False, saveto=False):
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, style, linewidth=linewidth)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.5, linestyle=':')
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

# Translations
def deg2rad(deg):
    return deg/180*np.pi
def rad2deg(rad):
    return rad/np.pi*180
def jd2time(jd):
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
def time2jd(time):
    return sum(jdcal.gcal2jd(time.year, time.month, time.day)) + time.hour / 24 + \
           time.minute / 1440 + time.second / 86400 + (time.microsecond / 1000000) / 86400
def juliandate(year, month, day, hour, minute, second):
    return fkt.ephemeris.juliandate(year, month, day, hour, minute, second)
def rv2oe(rv, mu, grad_req=False):
    out = fkt.translations.krv2oe(rv, mu, grad_req)
    return _return_if_grad_req(out, grad_req)
def oe2rv(oe, mu, grad_req=False):
    out = fkt.translations.koe2rv(oe, mu, grad_req)
    return _return_if_grad_req(out, grad_req)
def rv2ee(rv, mu, grad_req=False):
    out = fkt.translations.krv2ee(rv, mu, grad_req)
    return _return_if_grad_req(out, grad_req)
def ee2rv(ee, mu, grad_req=False):
    out = fkt.translations.kee2rv(ee, mu, grad_req)
    return _return_if_grad_req(out, grad_req)
def cart2sphere(cart):
    return fkt.translations.kcart2sphere(cart)
def sphere2cart(sphere):
    return fkt.translations.ksphere2cart(sphere)
def cart2latlon(cart):
    if len(cart.shape) == 2:
        return fkt.translations.kcart2latlon(cart)
    elif len(cart.shape) == 1:
        return fkt.translations.kcart2latlon(np.reshape(cart, (3, 1)))[:, 0]
def latlon2cart(latlon):
    if len(latlon.shape) == 2:
        return fkt.translations.klatlon2cart(latlon)
    elif len(latlon.shape) == 1:
        return fkt.translations.klatlon2cart(np.reshape(latlon, (2, 1)))[:, 0]
def itrs2gcrs(xitrs, jd, grad_req=False):
    out = fkt.translations.kitrs2gcrs(xitrs, jd)
    return _return_if_grad_req(out, grad_req)
def gcrs2itrs(xgcrs, jd, grad_req=False):
    out = fkt.translations.kgcrs2itrs(xgcrs, jd)
    return _return_if_grad_req(out, grad_req)
def scrs2pa(xscrs, jd, grad_req=False):
    out = fkt.translations.kscrs2pa(xscrs, jd)
    return _return_if_grad_req(out, grad_req)
def scrs2mer(xscrs, jd, grad_req=False):
    chunk = 10000
    dim = xscrs.shape[0]
    ncols = xscrs.shape[1]
    xmer = np.empty((dim, ncols))
    dxmer = np.empty((dim, dim, ncols))
    for i in range(int(np.ceil(ncols/chunk))):
        cols = np.arange(i*chunk, min((i+1)*chunk, ncols))
        xmer[:, cols], dxmer[:, :, cols] = fkt.translations.kscrs2mer(xscrs[:, cols], jd[cols])
    if grad_req:
        return xmer, dxmer
    else:
        return xmer
def mer2scrs(xmer, jd, grad_req=False):
    chunk = 10000
    dim = xmer.shape[0]
    ncols = xmer.shape[1]
    xscrs = np.empty((dim, ncols))
    dxscrs = np.empty((dim, dim, ncols))
    for i in range(int(np.ceil(ncols / chunk))):
        cols = np.arange(i * chunk, min((i + 1) * chunk, ncols))
        xscrs[:, cols], dxscrs[:, :, cols] = fkt.translations.kmer2scrs(xmer[:, cols], jd[cols])
    if grad_req:
        return xscrs, dxscrs
    else:
        return xscrs
def scrs2gcrs(xscrs, jd, DistUnit, VelUnit):
    return fkt.translations.kscrs2gcrs(xscrs, jd, DistUnit, VelUnit)
def gcrs2scrs(xgcrs, jd, DistUnit, VelUnit):
    return fkt.translations.kgcrs2scrs(xgcrs, jd, DistUnit, VelUnit)
def hcrs2gcrs(xhcrs, jd, DistUnit, VelUnit):
    return fkt.translations.khcrs2gcrs(xhcrs, jd, DistUnit, VelUnit)
def gcrs2hcrs(xgcrs, jd, DistUnit, VelUnit):
    return fkt.translations.kgcrs2hcrs(xgcrs, jd, DistUnit, VelUnit)
def scrs2sors(xscrs, jd, grad_req=False):
    out = fkt.translations.kscrs2sors(xscrs, jd)
    return _return_if_grad_req(out, grad_req)
def sors2scrs(xsors, jd, grad_req=False):
    out = fkt.translations.ksors2scrs(xsors, jd)
    return _return_if_grad_req(out, grad_req)
def ine2rot(xine, t, t0):
    return fkt.translations.kine2rot(xine, t, t0)
def rot2ine(xrot, t, t0):
    return fkt.translations.krot2ine(xrot, t, t0)
def ine2rotEph(xine, jd, first_body, secondary_body, DistUnit, VelUnit):
    return fkt.translations.kine2roteph(xine, jd, first_body, secondary_body, DistUnit, VelUnit)
def rot2ineEph(xrot, jd, first_body, secondary_body, DistUnit, VelUnit):
    return fkt.translations.krot2ineeph(xrot, jd, first_body, secondary_body, DistUnit, VelUnit)
def mer2lvlh(xmer, lat, lon):
    if len(xmer.shape) == 1 and xmer.shape[0] == 3:
        return fkt.translations.kmer2lvlh(xmer, lat, lon)
    elif len(xmer.shape) == 2 and xmer.shape[0] == 3:
        XLVLH = np.zeros((3, xmer.shape[1]))
        for i in range(xmer.shape[1]):
            XLVLH[:, i] = fkt.translations.kmer2lvlh(xmer[:, i], lat, lon)
        return XLVLH
    else:
        raise 'xmer should be a 3d vector or a 3xN matrix.'
def lvlh2mer(xlvlh, lat, lon):
    if len(xlvlh.shape) == 1 and xlvlh.shape[0] == 3:
        return fkt.translations.kxlvlh2mer(xlvlh, lat, lon)
    elif len(xlvlh.shape) == 2 and xlvlh.shape[0] == 3:
        XMER = np.zeros((3, xlvlh.shape[1]))
        for i in range(xlvlh.shape[1]):
            XMER[:, i] = fkt.translations.kxlvlh2mer(xlvlh[:, i], lat, lon)
        return XMER
    else:
        raise 'xlvlh should be a 3d vector or a 3xN matrix.'

# Units and constants.
def units(*args):
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
def astro_const():
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

# Equations of motion.
def r2bp(t, s):
    return fkt.equationsmodule.kr2bp(t, s)
def cr3bp_fb(t, s, mu, stm_req):
    fkt.equationsmodule.massparameter = mu
    fkt.equationsmodule.stm_required = stm_req
    return fkt.equationsmodule.cr3bp_fb(t, s)
def cr3bp_sb(t, s, mu, stm_req):
    fkt.equationsmodule.massparameter = mu
    fkt.equationsmodule.stm_required = stm_req
    return fkt.equationsmodule.cr3bp_sb(t, s)
def nbp_rv_earth(t, s, stm_req, sources, data, units_data):
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_rv_earth(t, s)
def nbp_rv_moon(t, s, stm_req, sources, data, units_data):
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_rv_moon(t, s)
def nbp_rvm_earth(t, s, stm_req, sources, data, units_data):
    pass
    # _set_nbp_parameters(stm_req, sources, data, units_data)
    # return fkt.equationsmodule.knbp_rvm_earth(t, s)
def nbp_rvm_moon(t, s, stm_req, sources, data, units_data):
    pass
    # _set_nbp_parameters(stm_req, sources, data, units_data)
    # return fkt.equationsmodule.knbp_rvm_moon(t, s)
def nbp_ee_earth(t, s, stm_req, sources, data, units_data):
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_ee_earth(t, s)
def nbp_ee_moon(t, s, stm_req, sources, data, units_data):
    _set_nbp_parameters(stm_req, sources, data, units_data)
    return fkt.equationsmodule.knbp_ee_moon(t, s)

# Propagation routines.
def propagate_nbp(central_body, tspan, x0, sources_dict, dat_dict, stm, variables):
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
def propagate_r2bp(tspan, x0):
    tspan, x0 = to_float(tspan, x0)
    t, y = fkt.propagationmodule.propagate_r2bp(tspan, x0)
    return t, y
def propagate_cr3bp(central_body, tspan, x0, mu, stm):
    tspan, x0, mu = to_float(tspan, x0, mu)
    neq = 42 if stm else 6
    t, y = fkt.propagationmodule.propagate_cr3bp(central_body.lower(), tspan, x0, mu, stm, neq)
    return t, y
def propagate_br4bp(central_body, tspan, x0, mu, GM4b, a4b, theta0, stm):
    tspan, x0, mu, GM4b, a4b, theta0 = to_float(tspan, x0, mu, GM4b, a4b, theta0)
    neq = 42 if stm else 6
    t, y = fkt.propagationmodule.propagate_br4bp(central_body.lower(), tspan, x0, mu, GM4b, a4b, theta0, stm, neq)
    return t, y
def prepare_sources_dict():
    sources = {}
    sources['sun'] = False
    sources['mercury'] = False
    sources['venus'] = False
    sources['earth'] = False
    sources['moon'] = False
    sources['mars'] = False
    sources['jupiter'] = False
    sources['saturn'] = False
    sources['uranus'] = False
    sources['neptune'] = False
    sources['srp'] = False
    sources['cmplxmoon'] = False
    sources['atm'] = False
    sources['j2'] = False
    return sources

# Visibility routines.
def is_visible(r_sat, lat_deg, long_deg, body_radius, threshold_deg):

    # r_sat:         3d-vector or 3 x N matrix
    # lat_deg:       scalar or vector
    # long_deg:      scalar or vector
    # threshold_deg: scalar or vector

    if len(r_sat.shape) == 1:
        if r_sat.shape[0] != 3:
            raise 'r_sat as a vector should have 3 components.'
        r_sat = np.reshape(r_sat, (3, 1), order='F')
    if r_sat.shape[0] != 3 or len(r_sat.shape) != 2:
        raise 'r_sat as a matrix should have N rows and 3 columns.'
    r_sat = r_sat.copy().T / body_radius

    if isinstance(lat_deg, (float, int)):
        lat_deg = np.array([lat_deg])
    if isinstance(long_deg, (float, int)):
        long_deg = np.array([long_deg])
    if len(lat_deg.shape) != 1:
        raise 'lat_deg should be a scalar or a vector.'
    if len(long_deg.shape) != 1:
        raise 'long_deg should be a scalar or a vector.'
    if lat_deg.shape[0] != long_deg.shape[0]:
        raise 'lat_deg and long_deg should have the same size.'
    lat_long = np.reshape(np.concatenate((lat_deg/180*np.pi, long_deg/180*np.pi), axis=0), (2, -1))

    threshold = threshold_deg / 180 * np.pi
    if isinstance(threshold, (float, int)):
        threshold = np.full((r_sat.shape[0],), threshold)
    if len(threshold.shape) != 1:
        raise 'threshold_deg should be a scalar or a vector'
    if threshold.shape[0] != r_sat.shape[0]:
        raise 'threshold_deg should have r_sat.shape[0] number of elements'

    fkt.visibilitymodule.r_sat = r_sat
    fkt.visibilitymodule.lat_long = lat_long
    fkt.visibilitymodule.threshold = threshold

    fkt.visibilitymodule.isvisible(r_sat.shape[0], lat_long.shape[1])

    status = fkt.visibilitymodule.status.copy()
    elev = fkt.visibilitymodule.elev.copy()
    azim = fkt.visibilitymodule.azim.copy()

    status[status == -1] = 1
    elev_deg = elev/np.pi*180
    azim_deg = azim/np.pi*180

    return status, elev_deg, azim_deg

# Save and load routines.
def save(variable, filename):
    with open(filename, 'wb') as f:
        pickle.dump(variable, f)
def load(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

# General astrodynamics.
def period_hours(altitude_km, body):
    ku = units(body)
    return 2*np.pi*np.sqrt((ku['DistUnit']+altitude_km)**3/ku['GM'])/3600.0
def altitude_km(period_hours, body):
    ku = units(body)
    period = period_hours*3600
    return (period**2/(4*np.pi**2)*ku['GM'])**(1/3) - ku['DistUnit']
def circular_velocity_km_s(altitude_km, body):
    ku = units(body)
    return np.sqrt(ku['GM']/(ku['DistUnit'] + altitude_km))
def dv_hohmann(r1_nondim, r2_nondim):
    dv1 = np.sqrt(1.0 / r1_nondim) * (np.sqrt(2 * r2_nondim / (r1_nondim + r2_nondim)) - 1)
    dv2 = np.sqrt(1.0 / r2_nondim) * (np.sqrt(2 * r1_nondim / (r1_nondim + r2_nondim)) - 1)
    dv_nondim = dv1 + dv2
    return dv_nondim
def tof_hohmann(r1_nondim, r2_nondim):
    a = (r1_nondim + r2_nondim) / 2
    tof_nondim = np.pi * (a ** 1.5)
    return tof_nondim

# Trofimov-Shirobokov model.
def get_order(altitude_thousands_km, approx_level='soft'):
    if approx_level == 'soft':
        return np.floor((25.0 / altitude_thousands_km)**0.8)+1
    elif approx_level == 'hard':
        return np.floor((40.0 / altitude_thousands_km)**0.8)+1
    else:
        raise Exception('Unknown approx_level.')

# Auxilary protected methods.
def _set_nbp_parameters(stm_req, sources, data, units_data):
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
def _return_if_grad_req(out, grad_req):
    if grad_req:
        return out
    else:
        return out[0]
def _sources_dict_to_vec(sources_dict):
    sources_vec = np.zeros((len(sources_dict),))
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
def _dat_dict_to_vec(dat_dict):
    dat_vec = np.zeros((4,))
    dat_vec[0] = dat_dict['jd_zero']
    dat_vec[1] = dat_dict['area']
    dat_vec[2] = dat_dict['mass']
    dat_vec[3] = dat_dict['order']
    return dat_vec
