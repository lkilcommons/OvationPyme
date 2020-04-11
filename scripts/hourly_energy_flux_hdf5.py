import datetime,os
from collections import OrderedDict

import numpy as np
import h5py

from ovationpyme.ovation_prime import FluxEstimator

def datetime_to_iso8601_str(dt):
    return dt.strftime('%Y%m%dT%H:%M:%S')

if __name__ == '__main__':
    energy_or_number = 'energy'
    estimators = OrderedDict() 
    for atype in ['diff','mono','wave','ions']:
        estimators[atype]=FluxEstimator(atype,energy_or_number)

    startdt = datetime.datetime(2003,1,1)
    enddt = datetime.datetime(2004,1,1)

    fn = 'ovationpyme_hourly_energy_flux_{}.h5'.format(startdt.year)
    homedir = os.path.expanduser('~/')
    h5fn = os.path.join(homedir,fn)

    if os.path.exists(h5fn):
        raise ValueError('File {} already exists'.format(h5fn))

    dt = startdt
    while dt < enddt:
        with h5py.File(h5fn) as h5f:
            dtstr = datetime_to_iso8601_str(dt)
            print(dtstr)
            for hemi in ['N','S']:
                
                ion_flux_outs = estimators['ions'].get_flux_for_time(dt,
                                                                    hemi=hemi)
                grid_mlats,grid_mlts,ion_energy_flux = ion_flux_outs
                
                if hemi == 'N' and 'grid_mlats' not in h5f:
                    h5f.create_dataset('grid_mlats',data=grid_mlats)
                    h5f.create_dataset('grid_mlts',data=grid_mlts)

                electron_energy_flux = np.zeros_like(ion_energy_flux)
                for atype in ['diff','mono','wave']:
                    get_flux_outs = estimators[atype].get_flux_for_time(dt,
                                                                        hemi=hemi)
                    auroral_type_electron_energy_flux = get_flux_outs[2]
                    electron_energy_flux += auroral_type_electron_energy_flux
                
                h5grp = h5f.create_group('/'+dtstr+'/'+hemi)
                h5grp.create_dataset('ion_energy_flux',data=ion_energy_flux)
                h5grp.create_dataset('electron_energy_flux',data=electron_energy_flux)
            
            dt+=datetime.timedelta(hours=1)

