import datetime

import numpy as np
import matplotlib.pyplot as plt

import ovation_utilities
from geospacepy import satplottools, special_datetime


if __name__=='__main__':
    startdt = datetime.datetime(2013,3,16,2,0,0)
    enddt = datetime.datetime(2013,3,16,4,0,0)
    dt = datetime.datetime(2013,3,16,3,0,0)

    startjd = special_datetime.datetime2jd(startdt)
    endjd = special_datetime.datetime2jd(enddt)
    jd = special_datetime.datetime2jd(dt)

    f = plt.figure(figsize=(8,11))
    sw = ovation_utilities.read_solarwind(dt)
    avgsw = ovation_utilities.calc_avg_solarwind(dt)

    inwindow = np.logical_and(sw['jd']>startjd,sw['jd']<endjd)

    n_plots = len(avgsw.keys())
    axs = [f.add_subplot(n_plots,1,i) for i in range(1,n_plots+1)]
    for i,swvar in enumerate(avgsw):
        axs[i].plot(sw['jd'],sw[swvar],'k-',label=swvar)
        axs[i].plot(jd,avgsw[swvar],'bo',label='Average')
        axs[i].set_ylabel(swvar)
        axs[i].legend()

    f.suptitle('{}-{}'.format(startdt.strftime('%Y%m%d %H:%M'),
                              enddt.strftime('%Y%m%d %H:%M')))

    f.savefig('visual_test_averaged_solarwind.png')


