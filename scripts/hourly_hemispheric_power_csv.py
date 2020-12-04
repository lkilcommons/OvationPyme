import datetime,os
from collections import OrderedDict
import numpy as np
from ovationpyme.ovation_prime import FluxEstimator
from geospacepy.spherical_geometry import grid_surface_integral

def datetime_to_iso8601_str(dt):
    return dt.strftime('%Y%m%dT%H:%M:%S')

if __name__ == '__main__':

    Re = 6371200
    mWtoW = 1/1000.
    WtoGW = 1/1.0e9

    energy_or_number = 'energy'
    atypes = ['diff','mono','wave','ions']
    hemis = ['N','S']
    estimators = OrderedDict() 
    for atype in atypes:
        estimators[atype]=FluxEstimator(atype,energy_or_number)

    startdt = datetime.datetime(2015,1,1)
    enddt = datetime.datetime(2016,1,1)-datetime.timedelta(days=1)

    fn = 'ovationpyme_hourly_hemispheric_power_{}.csv'.format(startdt.year)
    homedir = os.path.expanduser('~/')
    csvfn = os.path.join(homedir,fn)

    if os.path.exists(csvfn):
        raise ValueError('File {} already exists'.format(csvfn))

    csv_column_names = []
    csv_column_names.append('Time (ISO8601)')
    for atype in atypes:
        for hemi in hemis:
            csv_column_names.append(atype+'_'+hemi)

    with open(csvfn,'w') as f:
        f.write(','.join(csv_column_names)+'\n')
        dt = startdt
        while dt < enddt:
            dtstr = datetime_to_iso8601_str(dt)
            csv_row_data = []
            csv_row_data.append(dtstr)
            for atype in atypes:
                for hemi in ['N','S']:
                    outs = estimators[atype].get_flux_for_time(dt,hemi=hemi)
                    grid_mlats,grid_mlts,energy_flux = outs
                    #Integrate flux over bins
                    intflux = grid_surface_integral(grid_mlats,grid_mlts,energy_flux,
                                                    Re,'hour')
                    intflux = intflux*mWtoW*WtoGW # Convert to GW
                    csv_row_data.append('{:0.3f}'.format(intflux))
            csvrowline = ','.join(csv_row_data)+'\n' 
            f.write(csvrowline)
            print(csvrowline)

            dt+=datetime.timedelta(hours=1)

