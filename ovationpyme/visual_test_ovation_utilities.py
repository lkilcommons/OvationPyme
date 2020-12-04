import datetime,sys

import numpy as np
import matplotlib.pyplot as plt

import ovation_utilities
from geospacepy import satplottools, special_datetime

from logbook import StreamHandler

if __name__=='__main__':
    StreamHandler(sys.stdout).push_application()

    startdt = datetime.datetime(2013,3,16,2,0,0)
    enddt = datetime.datetime(2013,3,16,4,0,0)
    dt = datetime.datetime(2013,3,16,3,0,0)

    jd = special_datetime.datetime2jd(dt)

    f = plt.figure(figsize=(8,11))
    sw = ovation_utilities.read_solarwind(dt)
    sw4avg = ovation_utilities.hourly_solarwind_for_average(dt)
    avgsw = ovation_utilities.calc_avg_solarwind(dt)

    n_plots = len(avgsw.keys())
    axs = [f.add_subplot(n_plots,1,i) for i in range(1,n_plots+1)]
    for i,swvar in enumerate(avgsw):
        axs[i].plot(sw['jd'],sw[swvar],'k-',label=swvar)
        axs[i].plot(sw4avg['jd'],sw4avg[swvar],'b.',label='hour average')
        axs[i].plot(jd,avgsw[swvar],'go',label='5 hour weighted average')
        axs[i].set_ylabel(swvar)
        axs[i].legend()

    f.suptitle('{}-{}'.format(startdt.strftime('%Y%m%d %H:%M'),
                              enddt.strftime('%Y%m%d %H:%M')))

    f.savefig('visual_test_averaged_solarwind.png')


