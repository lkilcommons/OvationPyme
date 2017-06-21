# Ovation Pyme
## A pure-python implementation of the Ovation Prime auroral precipitation model

## Introduction to The Model
Ovation Prime is a model which predicts the total energy flux, total number flux,
and average characteristic energy of precipitating electrons and ions in the polar regions.
A Maxwellian particle energy (and velocity) distribution is assumed for the calculation of
average energy. The model was based on data from the Defense Meteorology Satellite Program
spacecraft, using data beginning in 1985. These spacecraft carry a particle detector (SSJ)
which is sensitive to particles with characteristic energies between 30 eV and 30 keV.
The original model was created by Patrick Newell et. al. (John Hopkins Applied Physics Laboratory)

## Provenance of this implementation
Ovation Pyme is a complete translation of the IDL (a proprietary programming language,
owned by ExelisVis) version used by NOAA, which was released in an open source format 
[Sourceforge](https://sourceforge.net/projects/ovation-prime/) by Rob Redmon and Janet Machol
of NOAA National Center for Environmental Information (NCEI). 
Liam Kilcommons created this software to support the Pythonification of the Assimiliative Mapping
of Ionospheric Electrodynamics project (AMIEPy), and for the Python space science community. 
__Contributions and comments are very welcome__.

## Differences in this implementation compared to the IDL version
1. I correct what I belive to be an oversight in how flux is averaged between
the Northern and Southern hemispheres. In the original code, the hemispheres 
were combined (mean) after the fluxes for the various seasons were averged using a weight
appropriate to the day of the year the model was called for. This does not account for
the inverse hemispheric response to a particular season (northern activity maximizes in the
summer and southern activity in the winter). In OvationPyme, the seasonal weighting for day of year
DOY is done on both hemispheres seperately with the 
southern hemisphere weighting calculated for DOY_S = 365-DOY, and then the hemispheres are combined. 

2. Compute estimated ionospheric height-integrated conductivity (conductance) using 
the technique used by Ellen Cousins et. al. (2015), which uses the Robinson emperical
relationship mapping conductance to a function of total electron energy flux and average energy, using
Ovation's diffuse aurora.

## Future Features To Be Added (beyond IDL version):
2. Spatial numerical derivatives of conductance (derivatives in magnetic latitude and magnetic local time directions).

## Verison Restrictions
Use Python 2.7 (not 3) for the moment. Python 3 support is pending.

## Installation Instructions
1. Clone or download the (geospacepy-lite)[https://github.com/lkilcommons/geospacepy-lite] library
2. Edit the OMNI data download path in geospacepy-lite file geospacepy-config.py
3. From the geospacepy-lite directory: `python setup.py install`
4. Clone or download the OvationPyme repostiory
5. From the OvationPyme directory: `python setup.py install`

## Tests
Unit tests are written for the py.test framework. If you have this installed,
you can run the tests by issuing `py.test` from the command line in the 'ovationpyme'
directory. If you are not a programmer, this probably will not be very interesting.

Functional tests (which make several plots to show things are working properly) 
can be run by calling the package's test_plots.py as a script.
i.e. from the command line run:
`python ovationpyme/functionaltests.py`

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
precipitation budget, J. Geophys. Res., 114, A09207, doi:10.1029/2009JA014326.Newell, P. T., T. Sotirelis, and S. Wing (2010), Seasonal variations in diffuse, monoenergetic, and broadband aurora, J. Geophys. Res., 115, A03216, doi:10.1029/2009JA014805.

- Sotirelis, T. and P. T. Newell (2000), "Boundary-oriented electron precipitation model," J. Geophys. Res.,
105 (A8), 18,655-18,673.
