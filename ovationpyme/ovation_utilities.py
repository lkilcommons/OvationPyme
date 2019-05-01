"""

"""
import numpy as np
import datetime
from geospacepy import omnireader, special_datetime

def calc_avg_solarwind(dt,oi=None):
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

    if oi is None:
        oi = omnireader.omni_interval(dt-datetime.timedelta(hours=n_hours+2),dt+datetime.timedelta(hours=n_hours+1),'5min',silent=True)
    
    if oi.cadence == 'hourly':
        velvar,densvar = 'V','N'
    else:
        velvar,densvar = 'flow_speed','proton_density'

    #Find closest time to dt in julian date instead of datetime (comparisons with arrays of datetime are slow)
    jd = special_datetime.datetime2jd(dt)
    om_jd = special_datetime.datetimearr2jd(oi['Epoch'])
    Bx,By,Bz = oi['BX_GSE'],oi['BY_GSM'],oi['BZ_GSM']
    V,Ni = oi[velvar],oi[densvar]
    Ec = calc_coupling(Bx,By,Bz,V)

    i_first = np.nanargmin(np.abs(om_jd-jd))
    
    #Get indicies of all data in the average
    om_in_avg = range(i_first-n_hours,i_first+1)
    weights = [prev_hour_weight*n_hours_back for n_hours_back in reversed(range(n_hours+1))] #reverse the list
    
    #Calculate weighted averages
    avgsw = dict()
    avgsw['Bx'] = np.nansum(Bx[om_in_avg]*weights)/len(om_in_avg)
    avgsw['By'] = np.nansum(By[om_in_avg]*weights)/len(om_in_avg)
    avgsw['Bz'] = np.nansum(Bz[om_in_avg]*weights)/len(om_in_avg)
    avgsw['V'] = np.nansum(V[om_in_avg]*weights)/len(om_in_avg)
    avgsw['Ec'] = np.nansum(Ec[om_in_avg]*weights)/len(om_in_avg)

    return avgsw

def calc_coupling(Bx,By,Bz,V):
    """
    
    Empirical Formula for dF/dt
    (i.e. the Newell coupling)

    """
    Ec = np.zeros_like(Bx)
    Ec.fill(np.nan)
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    BT = np.sqrt(By**2 + Bz**2)
    bztemp = Bz
    bztemp[Bz == 0] = .001
    #Caculate clock angle (theta_c = t_c)
    tc = np.arctan2(By,bztemp)
    neg_tc = BT*np.cos(tc)*Bz < 0 
    tc[neg_tc] = tc[neg_tc] + np.pi
    sintc = np.abs(np.sin(tc/2.))
    Ec = (V**1.33333)*(sintc**2.66667)*(BT**0.66667)
    return Ec

