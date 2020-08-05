"""

"""
import datetime
from collections import OrderedDict

import numpy as np
import functools

from geospacepy import special_datetime, sun
from nasaomnireader.omnireader import omni_interval

_ovation_prime_omni_cadence = 'hourly' #Ovation Prime was created using hourly SW

def cache_omni_interval(func):
    """Decorator which decorates functions with call signature
    func(dt,oi) which calculate something from a given omni interval
    Implements on-the-fly creation of an omni_interval, cacheing it
    as a function parameter, and then creating a new one if requested
    dateimte is out of range
    """
    cache = {}

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

        if 'omni_interval' in cache:
            cached_oi = cache['omni_interval']
            need_new_oi = not _dt_within_range(dt,cached_oi)
        else:
            need_new_oi = True

        if need_new_oi:
            startdt = dt-datetime.timedelta(days=new_interval_days_before_dt)
            enddt = dt+datetime.timedelta(days=new_interval_days_after_dt)

            oi = omni_interval(startdt,
                                            enddt,
                                            _ovation_prime_omni_cadence,
                                            silent=True)
            #Save to cache
            cache['omni_interval'] = oi
            # print("Created new solar wind interval: {}-{}".format(oi.startdt,
            #                                                         oi.enddt))
        else:
            #Load from cache
            oi = cache['omni_interval']
            
            # print("Using cached solar wind interval: {}-{}".format(oi.startdt,
            #                                                         oi.enddt))


        return func(dt,oi)

    return cache_omni_interval_wrapper

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

@cache_omni_interval
def read_solarwind(dt,oi):
    #Find closest time to dt in julian date instead of datetime
    #(comparisons with arrays of datetime are slow)
    if oi.cadence == 'hourly':
        velvar, densvar = 'V', 'N'
    else:
        velvar, densvar = 'flow_speed', 'proton_density'

    jd = special_datetime.datetimearr2jd(oi['Epoch']).flatten()
    Bx, By, Bz = oi['BX_GSE'], oi['BY_GSM'], oi['BZ_GSM']
    V,Ni = oi[velvar],oi[densvar]
    Ec = calc_coupling(Bx, By, Bz, V)

    sw = OrderedDict()
    sw['jd']=jd
    sw['Bx']=Bx
    sw['By']=Bx
    sw['Bz']=Bx
    sw['V']=Bx
    sw['Ni']=Ni
    sw['Ec']=Ec
    return sw

@cache_omni_interval
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
    n_hours = 4       # hours previous to integrate over
    prev_hour_weight = 0.65    # reduce weighting by factor of wh each hour back

    jd = special_datetime.datetime2jd(dt)

    sw = read_solarwind(dt)
    om_jd = sw['jd']

    i_first = np.nanargmin(np.abs(om_jd-jd))

    #Get indicies of all data in the average
    om_in_avg = list(range(i_first-n_hours, i_first+1))
    weights = [prev_hour_weight*n_hours_back for n_hours_back in reversed(list(range(n_hours+1)))] #reverse the list

    #Calculate weighted averages
    avgsw = OrderedDict()
    avgsw['Bx'] = np.nansum(sw['Bx'][om_in_avg]*weights)/len(om_in_avg)
    avgsw['By'] = np.nansum(sw['By'][om_in_avg]*weights)/len(om_in_avg)
    avgsw['Bz'] = np.nansum(sw['Bz'][om_in_avg]*weights)/len(om_in_avg)
    avgsw['V'] = np.nansum(sw['V'][om_in_avg]*weights)/len(om_in_avg)
    avgsw['Ec'] = np.nansum(sw['Ec'][om_in_avg]*weights)/len(om_in_avg)

    return avgsw

@cache_omni_interval
def get_daily_f107(dt,oi):
    """
    Since OvationPyme uses hourly OMNI data
    I just do the mean for all of the 1 hour values for the day
    of dt (used for calculating solar conductance)
    """
    omjd = special_datetime.datetimearr2jd(oi['Epoch']).flatten()
    omf107 = oi['F10_INDEX']
    jd = special_datetime.datetime2jd(dt)
    imatch = omjd==np.floor(jd)
    return np.nanmean(omf107[imatch])

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
