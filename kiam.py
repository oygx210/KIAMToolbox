import FKIAMToolbox as fkt
import jdcal
import datetime
import math

# Translations.
def jd2time(jd):
    gcal = jdcal.jd2gcal(2400000.5, jd - 2400000.5)
    frac_hours = gcal[3]*24
    hours = math.floor(frac_hours)
    frac_minutes = (frac_hours - hours)*60
    minutes = math.floor(frac_minutes)
    frac_seconds = (frac_minutes - minutes)*60
    seconds = int(round(frac_seconds, 0))
    if seconds != 60:
        return datetime.datetime(gcal[0], gcal[1], gcal[2], hours, minutes, seconds)
    else:
        return datetime.datetime(gcal[0], gcal[1], gcal[2], hours, minutes) + datetime.timedelta(minutes=1)
def time2jd(time):
    return sum(jdcal.gcal2jd(time.year, time.month, time.day)) + time.hour/24 + \
                             time.minute/1440 + time.second/86400 + (time.microsecond/1000000)/86400
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
    return fkt.translations.sphere2cart(sphere)
def cart2latlon(cart):
    return fkt.translations.kcart2latlon(cart)
def latlon2cart(latlon):
    return fkt.translations.latlon2cart(latlon)
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
    out = fkt.translations.kscrs2mer(xscrs, jd)
    return _return_if_grad_req(out, grad_req)
def mer2scrs(xmer, jd, grad_req=False):
    out = fkt.translations.kmer2scrs(xmer, jd)
    return _return_if_grad_req(out, grad_req)
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

# Units and constants.
def units(*args):
    units_info = {}
    if len(args) == 1:
        output = fkt.constantsandunits.kunits_onebody(args[0])
        units_info['GM'] = output[0]
        units_info['DistUnit'] = output[1]
        units_info['VelUnit'] = output[2]
        units_info['TimeUnit'] = output[3]
        units_info['AccUnit'] = output[4]
    elif len(args) == 2:
        output = fkt.constantsandunits.kunits_twobody(args[0], args[1])
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

    moon['Moon']['OrbitsAround'] = fkt.constantsandunits.neptune_orbitsaround
    moon['Moon']['GM'] = fkt.constantsandunits.neptune_gm
    moon['Moon']['MeanRadius'] = fkt.constantsandunits.neptune_meanradius
    moon['Moon']['SemimajorAxis'] = fkt.constantsandunits.neptune_semimajoraxis

    small_body['Pluto']['OrbitsAround'] = fkt.constantsandunits.neptune_orbitsaround
    small_body['Pluto']['GM'] = fkt.constantsandunits.neptune_gm
    small_body['Pluto']['MeanRadius'] = fkt.constantsandunits.neptune_meanradius
    small_body['Pluto']['EquatorRadius'] = fkt.constantsandunits.neptune_equatorradius
    small_body['Pluto']['SemimajorAxis'] = fkt.constantsandunits.neptune_semimajoraxis

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