# Ovation Pyme
## A pure-python implementation of the Ovation Prime 2010 auroral precipitation model

You can run the notebooks in this repo on the cloud by clicking below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lkilcommons/OvationPyme/HEAD?labpath=notebooks%2FInteractiveWithHemisphericPower.ipynb)

## Introduction to The Model
Ovation Prime 2010, (described in [Newell et al., 2009](https://doi.org/10.1029/2009JA014326) and [Newell et al. 2010](https://doi.org/10.1029/2009JA014805)) is a model which predicts the total energy flux, total number flux,
and average characteristic energy of precipitating electrons and ions in the polar regions.
A Maxwellian particle energy (and velocity) distribution is assumed for the calculation of
average energy. The model was based on data from the Defense Meteorology Satellite Program
spacecraft, using data beginning in 1985. These spacecraft carry a particle detector (SSJ)
which is sensitive to particles with characteristic energies between 30 eV and 30 keV.

## Provenance of this implementation
Ovation Pyme is a complete translation of the IDL (a proprietary programming language,
owned by ExelisVis) version released in an open source format on
[Sourceforge](https://sourceforge.net/projects/ovation-prime/) by Janet Machol, Rob Redmon and Nathan Case
of NOAA National Center for Environmental Information (NCEI).
__Contributions and comments are very welcome__.

## Problems and limitations of Ovation Prime 2010
These limitations are described in detail in [Newell et al. 2013](https://doi.org/10.1002/2014SW001056)
1. Model produces large spurious values in some isolated bins (salt-and-pepper noise)
2. Model is most valid for low to moderate geomagnetic activity (Kp<=6)
3. Model is paramaterized to use hourly solar wind data (only one unique result per hour)

## Differences in this implementation compared to the IDL version

1. Compute estimated ionospheric height-integrated conductivity (conductance) using 
the technique used by Ellen Cousins et. al. (2015), which uses the Robinson emperical
relationship mapping conductance to a function of total electron energy flux and average energy, using
Ovation's diffuse aurora.
2. Interpolation of data to arbitrary latitude and longitude grids

## Verison Restrictions
Use Python versions >= 2.7

## Installation Instructions
1. Clone or download the [nasaomnireader](https://github.com/lkilcommons/nasaomnireader) library
3. From the nasaomnireader directory: `python setup.py install`
4. Clone or download the OvationPyme repostiory
5. From the OvationPyme directory: `python setup.py install`

## Tests
Unit tests are written for the py.test framework. If you have this installed,
you can run the tests by issuing `py.test` from the command line in the 'ovationpyme'
directory.

Functional tests (which make several plots to show things are working properly) 
can be run by calling the package's test_plots.py as a script.
i.e. from the command line run:
`python ovationpyme/visual_test_ovation_prime.py`

The current test plots are:

1. A plot of the Northern and Southern polar electron energy flux for a fixed solar coupling value for one of the seasons (summer)

2. A plot of the combined (averaged) northern and southern electron energy flux for the hemispherically appropriate summer
(the data for boreal summer for the north, and for boreal winter (which is austral summer), for the south. This is traditionally how the model is run (combined hemispheres). 

3. A plot of the electron energy flux for 
combined hemispheres for a particular time. This
tests the ability of the code to automatically download solar wind data from the NASA Omniweb FTP server and calculate the 
Newell solar wind coupling function (see references). 

4. A plot of the ionospheric hall and pedersen conductances for the Northern Hemisphere.

## References

- Cousins, E. D. P., T. Matsuo, and A. D. Richmond (2015), Mapping high-latitude ionospheric electrodynamics with SuperDARN and AMPERE, J. Geophys. Res. Space Physics, 120, 5854–5870, doi:10.1002/2014JA020463.

- Machol, J. L., J. C. Green, R. J. Redmon, R. A. Viereck, and P. T. Newell (2012), Evaluation of OVATION
Prime as a Forecast Model for Visible Aurorae, Space Weather, 10, S03005,
doi:10.1029/2011SW000746.

- Newell, P.T., T. Sotirelis, K. Liou, C.‐I. Meng, and F.J. Rich (2007), A nearly universal solar wind‐
magnetosphere coupling function inferred from 10 magnetospheric state variables, J. Geophys. Res.,
112, A01206, doi:10.1029/2006JA012015.

- Newell, P. T., T. Sotirelis, and S. Wing (2009), Diffuse, monoenergetic, and broadband aurora: The global
precipitation budget, J. Geophys. Res., 114, A09207, doi:10.1029/2009JA014326.

- Newell, P. T., T. Sotirelis, and S. Wing (2010), Seasonal variations in diffuse, monoenergetic, and broadband aurora, J. Geophys. Res., 115, A03216, doi:10.1029/2009JA014805.

- Sotirelis, T. and P. T. Newell (2000), "Boundary-oriented electron precipitation model," J. Geophys. Res.,
105 (A8), 18,655-18,673.
