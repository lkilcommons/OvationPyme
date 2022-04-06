"""

"""
import datetime
from collections import OrderedDict

import numpy as np
import functools

from geospacepy import special_datetime, sun
from nasaomnireader.omnireader import omni_interval
from logbook import Logger
log = Logger('OvationPyme.ovation_utilites')

#_ovation_prime_omni_cadence = 'hourly' #Ovation Prime was created using hourly SW

def cache_omni_interval(cadence):
    """Decorator which decorates functions with call signature
    func(dt,oi) which calculate something from a given omni interval
    Implements on-the-fly creation of an omni_interval, cacheing it
    as a function parameter, and then creating a new one if requested
    dateimte is out of range
    """
    cache = {}

    def cache_omni_interval_decorator(func):
        
        tol_hrs_before=4
        tol_hrs_after=1
        new_interval_days_before_dt = 1.5
        new_interval_days_after_dt = 1.5

        def _dt_within_range(dt,oi):
            """
            Check that dt is more than tol_hrs_before hours after the start of
            the omni_interval, and more that tol_hrs_after before the end of
            """
            st_hrs_before_dt = (dt-oi.startdt).total_seconds()/3600.
            ed_hrs_after_dt = (oi.enddt-dt).total_seconds()/3600.
            in_before_range = st_hrs_before_dt > tol_hrs_before
            in_after_range = ed_hrs_after_dt > tol_hrs_after
            return in_before_range and in_after_range

        def cache_omni_interval_wrapper(dt):

            #print("Cached OMNI called for {}".format(dt))

            if 'omni_interval_{}'.format(cadence) in cache:
                cached_oi = cache['omni_interval_{}'.format(cadence)]
                need_new_oi = not _dt_within_range(dt,cached_oi)
            else:
                need_new_oi = True

            if need_new_oi:
                startdt = dt-datetime.timedelta(days=new_interval_days_before_dt)
                enddt = dt+datetime.timedelta(days=new_interval_days_after_dt)

                oi = omni_interval(startdt,enddt,cadence,silent=True)

                #Save to cache
                cache['omni_interval_{}'.format(cadence)] = oi
                log.debug("Created new solar wind interval: {}-{}".format(oi.startdt,
                                                                        oi.enddt))
            else:
                #Load from cache
                oi = cache['omni_interval_{}'.format(cadence)]
                
                log.debug("Using cached solar wind interval: {}-{}".format(oi.startdt,
                                                                        oi.enddt))


            return func(dt,oi)

        return cache_omni_interval_wrapper

    return cache_omni_interval_decorator

def calc_coupling(Bx, By, Bz, V):
    """
    Empirical Formula for dF/dt
    (i.e. the Newell coupling)
    """
    Ec = np.zeros_like(Bx)
    Ec.fill(np.nan)
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    BT = np.sqrt(By**2 + Bz**2)
    bztemp = Bz
    bztemp[Bz == 0] = 0.001
    #Caculate clock angle (theta_c = t_c)
    tc = np.arctan2(By,bztemp)
    neg_tc = BT*np.cos(tc)*Bz < 0
    tc[neg_tc] = tc[neg_tc] + np.pi
    sintc = np.abs(np.sin(tc/2.))
    Ec = (V**1.33333)*(sintc**2.66667)*(BT**0.66667)
    return Ec

@cache_omni_interval('1min')
def read_solarwind(dt,oi):
    """Get the solar wind parameters involved in the Newell coupling
    function at an hourly cadence (regardless of the cadence of the
    omni_interval input oi)
    """
    if oi.cadence == 'hourly':
        velvar, densvar = 'V', 'N'
    else:
        velvar, densvar = 'flow_speed', 'proton_density'

    #Variable name/key for return dict and associated omni_interval key
    swvars = OrderedDict(Bx='BX_GSE',
                         By='BY_GSM',
                         Bz='BZ_GSM',
                         V=velvar,
                         Ni=densvar)

    sw = OrderedDict() 
    sw['jd']=special_datetime.datetimearr2jd(oi['Epoch']).flatten()
    for swkey,oikey in swvars.items():
        sw[swkey]=oi[oikey]

    #Newell coupling    
    sw['Ec']=calc_coupling(sw['Bx'], sw['By'], sw['Bz'], sw['V'])
    return sw

@cache_omni_interval('1min')
def hourly_solarwind_for_average(dt,oi):
    """
    Takes a solarwind (sw) OrderedDict (output of read_solarwind)
    made using omni data at a sub-hourly cadence (5min or 1min),
    and averages the data to an hourly cadence relative to time dt
    (dt must be within the range of the data in sw). Returns a
    new OrderedDict with 'n_hours_in_average' values for each
    solar wind parameter
    """
    n_hours_in_average = 4 #number of hourly datapoints (4 previous)
    
    target_jd = special_datetime.datetime2jd(dt)
    
    sw = read_solarwind(dt)

    #Julian date in days to time relative to target time in hours
    #with positive values indicating time before the target
    hours_before_target = -1*(sw['jd']-target_jd)*24.
    
    sw4avg = OrderedDict()
    for swkey,swdata in sw.items():
        hourly_swdata = []
        for hour in range(n_hours_in_average)[::-1]:
            hourmask = np.logical_and(hours_before_target>=hour,
                                      hours_before_target<(hour+1))    
            if swkey == 'jd':
                hourly_swdata.append(np.nanmax(swdata[hourmask]))
            else:
                hourly_swdata.append(np.nanmean(swdata[hourmask]))
        sw4avg[swkey]=np.array(hourly_swdata)
    return sw4avg

@cache_omni_interval('1min')
def calc_avg_solarwind(dt,oi):
    """
    Calculates a weighted average of several
    solar wind variables n_hours (4 by default) backward
    in time from the closest hourly OMNIWeb
    datum to datetime dt

    oi is an optional omnireader.omni_interval
    instance from which to read the data. If this
    is None (default), will create a new omni_interval
    """
    prev_hour_weight=0.65

    sw4avg = hourly_solarwind_for_average(dt)
    n = sw4avg['jd'].size #Number of hourly datapoints to be averaged
    weights = [prev_hour_weight**n_hours_back for n_hours_back in range(n)[::-1]] #reverse the range

    #Calculate weighted averages
    avgsw = OrderedDict()
    avgsw['Bx'] = np.nansum(sw4avg['Bx']*weights)/np.sum(weights)
    avgsw['By'] = np.nansum(sw4avg['By']*weights)/np.sum(weights)
    avgsw['Bz'] = np.nansum(sw4avg['Bz']*weights)/np.sum(weights)
    avgsw['V'] = np.nansum(sw4avg['V']*weights)/np.sum(weights)
    avgsw['Ec'] = np.nansum(sw4avg['Ec']*weights)/np.sum(weights)

    return avgsw

@cache_omni_interval('hourly')
def get_daily_f107(dt,oi):
    """
    Since OvationPyme uses hourly OMNI data
    I just do the mean for all of the 1 hour values for the day
    of dt (used for calculating solar conductance)
    """
    omjd = special_datetime.datetimearr2jd(oi['Epoch']).flatten()
    omf107 = oi['F10_INDEX']
    jd = special_datetime.datetime2jd(dt)
    imatch = np.nanargmin(np.abs(omjd-jd))
    return omf107[imatch]

def calc_dF(dt):
    """dF==newell coupling for Ovation Prime"""
    return calc_avg_solarwind(dt)['Ec']

def robinson_auroral_conductance(numflux, eavg):
    """Robinson empirical formula for auroral conductance from
    energy flux and average energy of precipitating electrons

    INPUTS
    ------
        numflux, np.ndarray
            Electron number flux in #/cm^2/s
        eavg, np.ndarray
            Electron average energy in keV

    RETURNS
    -------
        sigp, np.ndarray
            Pedersen conductance [Mho]
        sigh, np.ndarray
            Hall conductance [Mho]
    """
    #From E. Cousins IDL code
    #Implement the Robinson formula
    #Assume all of the particles come in at the average energy??
    #Under maxwellian assumption (eavg=eflux/nflux), so this is valid
    keV_to_ergs = 1.6022e-9
    energyflux = numflux*eavg*keV_to_ergs
    #energyflux_grid *= 1.6022e-9
    sigp = 40.*eavg/(16+eavg**2) * np.sqrt(energyflux)
    sigh = 0.45*eavg**0.85*sigp
    return sigp,sigh

def brekke_moen_solar_conductance(dt,glats,glons,f107):
    """
    Estimate the solar conductance using methods from:

    Cousins, E. D. P., T. Matsuo, and A. D. Richmond (2015), Mapping
    high-latitude ionospheric electrodynamics with SuperDARN and AMPERE

    --which cites--

    Asgeir Brekke, Joran Moen,
    Observations of high latitude ionospheric conductances

    INPUTS
    ------
        dt, datetime.datetime
            Universal date and time to calculate solar conductances for

        glats, np.ndarray
            GeoDETIC latitudes

        glons, np.ndaray
            Geographic longitudes

    OUTPUTS
    -------

        sigp, np.ndarray
            Pedersen conductance at locations (glats,glons)

        sigh, np.ndarray
            Hall conductance at locations (glats,glons)

        f107, float
            F10.7 index (daily) used to calcuate solar conductance
    """
    szas_rad = sun.solar_zenith_angle(special_datetime.datetime2jd(dt), glats, glons)
    szas = np.rad2deg(szas_rad)

    sigp,sigh = np.zeros_like(glats),np.zeros_like(glats)

    cos65 = np.cos(65/180.*np.pi)
    sigp65  = .5*(f107*cos65)**(2./3)
    sigh65  = 1.8*np.sqrt(f107)*cos65
    sigp100 = sigp65-0.22*(100.-65.)

    in_band = szas <= 65.
    #print "%d/%d Zenith Angles < 65" % (np.count_nonzero(in_band),len(in_band))
    sigp[in_band] = .5*(f107*np.cos(szas_rad[in_band]))**(2./3)
    sigh[in_band] = 1.8*np.sqrt(f107)*np.cos(szas_rad[in_band])

    in_band = np.logical_and(szas >= 65.,szas < 100.)
    #print "%d/%d Zenith Angles > 65 and < 100" % (np.count_nonzero(in_band),len(in_band))
    sigp[in_band] = sigp65-.22*(szas[in_band]-65.)
    sigh[in_band] = sigh65-.27*(szas[in_band]-65.)

    in_band = szas > 100.
    sigp[in_band] = sigp100-.13*(szas[in_band]-100.)
    sigh[in_band] = sigh65-.27*(szas[in_band]-65.)

    sigp[sigp<.4] = .4
    sigh[sigh<.8] = .8

    #correct for inverse relationship with magnetic field from AMIE code
    #(conductance_models.f90)
    theta = np.radians(90.-glats)
    bbp = np.sqrt(1. - 0.99524*np.sin(theta)**2)*(1. + 0.3*np.cos(theta)**2)
    bbh = np.sqrt(1. - 0.01504*(1.-np.cos(theta)) - 0.97986*np.sin(theta)**2)*(1.0+0.5*np.cos(theta)**2)
    sigp = sigp*1.134/bbp
    sigh = sigh*1.285/bbh

    return sigp,sigh
