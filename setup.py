# Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Jul, 2016
# (C) 2016 Liam Kilcommons and AMIEPy Project Group (Tomoko Matsuo, PI)

import os
import glob

os.environ['DISTUTILS_DEBUG'] = "1"

from setuptools import setup, Extension
from setuptools.command import install as _install

setup(name='ovationpyme',
      version = "0.1.0",
      description = "Ovation Prime Auroral Model",
      #author = "VT SuperDARN Lab and friends",
      #author_email = "ajribeiro86@gmail.com",
      author = "Liam Kilcommons",
      author_email = 'liam.kilcommons@colorado.edu',
      url = "https://github.com/lkilcommons/ovationpyme",
      download_url = "https://github.com/lkilcommons/ovationpyme",
      long_description = "This is a Python port of the IDL implementation of Ovation Prime"+\
      " auroral particle flux model originally created by Newell et al (JHU)"+\
      " and packaged on Sourceforge by Redmon (NOAA NCEI), Machol, and Case "+\
      " for more information visit: https://sourceforge.net/projects/ovation-prime/",
      install_requires=['numpy','matplotlib','aacgmv2','geospacepy','logbook','scipy'],
      packages=['ovationpyme'],
      package_dir={'ovationpyme' : 'ovationpyme'},
      package_data={'ovationpyme': ['data/premodel/*.txt']}, #data names must be list
      license='LICENSE.txt',
      zip_safe = False,
      classifiers = [
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Programming Language :: Python"
            ]
      )
