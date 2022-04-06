import datetime
import pytest
import numpy as np
from numpy import testing as nptest
from ovationpyme.ovation_utilities import read_solarwind 
"""
Unit Tests for Ovation Prime Utilities
"""

#Basic does it run test
def test_read_solarwind():
    dt = datetime.datetime(2000,1,1,12,34,51)
    sw = read_solarwind(dt)
    assert 'Ec' in sw

