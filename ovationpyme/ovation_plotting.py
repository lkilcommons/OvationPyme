import numpy as np
import matplotlib.pyplot as plt

def latlt2polar(lat,lt,hemisphere):
    """
    Converts an array of latitude and lt points to polar for a top-down dialplot (latitude in degrees, LT in hours)
    i.e. makes latitude the radial quantity and MLT the azimuthal 

    get the radial displacement (referenced to down from northern pole if we want to do a top down on the north, 
    or up from south pole if visa-versa)
    """
    from numpy import pi
    if hemisphere=='N':
        r = 90.-lat
    elif hemisphere=='S':
        r = 90.-(-1*lat)
    else:
        raise ValueError('%s is not a valid hemisphere, N or S, please!' % (hemisphere))
    #convert lt to theta (azimuthal angle) in radians
    theta = lt/24. * 2*pi

    #the pi/2 rotates the coordinate system from
    #theta=0 at negative y-axis (local time) to
    #theta=0 at positive x axis (traditional polar coordinates)
    return r,theta

def polar2dial(ax):
    """
    Turns a matplotlib axes polar plot into a dial plot
    """
    #Rotate the plot so that noon is at the top and midnight
    #is at the bottom, and fix the labels so radial direction
    #is latitude and azimuthal direction is local time in hours
    ax.set_theta_zero_location('S')
    theta_label_values = np.array([0.,3.,6.,9.,12.,15.,18.,21.])*180./12
    theta_labels = ['%d:00' % (int(th/180.*12)) for th in theta_label_values.flatten().tolist()]
    ax.set_thetagrids(theta_label_values,labels=theta_labels)

    r_label_values = 90.-np.array([80.,70.,60.,50.])
    r_labels = [r'$%d^{o}$' % (int(90.-rv)) for rv in r_label_values.flatten().tolist()]
    ax.set_rgrids(r_label_values,labels=r_labels)
    ax.set_rlim([0.,40.])

def pcolor_flux(ax,mlatgrid,mltgrid,fluxgrid,hemisphere,**pcolor_kwargs):
    mlats,mlts = mlatgrid.flatten(),mltgrid.flatten()
    flux = fluxgrid.flatten()
    if 'vmin' not in pcolor_kwargs:
        pcolor_kwargs['vmin'] = np.nanpercentile(flux,5)
    if 'vmax' not in pcolor_kwargs:
        pcolor_kwargs['vmax'] = np.nanpercentile(flux,95)
    r,theta = latlt2polar(mlats,mlts,hemisphere)
    rgrid = r.reshape(mlatgrid.shape)
    thetagrid = theta.reshape(mltgrid.shape)
    mappable = ax.pcolormesh(thetagrid,rgrid,fluxgrid,**pcolor_kwargs)
    return mappable