import pytest

import numpy as np
from numpy import testing as nptest

import ovationpyme
"""
Unit Tests for Ovation Prime Main Module

IDL script create_test_results.pro output:

Ran season_epoch_ec for atype=diff,jtype=e_E_flux,(year,doy,sod)=    2011     103    3600
je[      24,      60]:    0.0620033
dF_bin(      3134.17)=       5
prob_estimate: dF=      3134.17 i_mlt=      24 j_mlat=     140=      0.897351
b1p[      24,     140]:     0.826437
b2p[      24,     140]:  2.26261e-05
b1a[      24,     140]:    0.0634671
b2a[      24,     140]: -2.70994e-06
"""

@pytest.fixture()
def idl_call_results(request):
    results = {
        'atype': 'diff',
        'season': 'winter',
        'i_mlt': 24,
        'j_mlat': 140,
        'j_mlat_NS': 60,
        'je': 0.0549602,
        'b1p': 0.826437,
        'b2p': 2.26261e-05,
        'b1a': 0.0634671,
        'b2a': -2.70994e-06,
        'dF': 3134.17,
        'prob': 0.897351
        }
    return results

@pytest.fixture()
def seasonal_flux_estimator(request):
    season, atype, energy_or_number = 'winter', 'diff', 'energy'
    return ovationpyme.ovation_prime.SeasonalFluxEstimator(season, atype, energy_or_number)

@pytest.fixture()
def flux_estimator(request):
    atype, energy_or_number = 'diff', 'energy'
    return ovationpyme.ovation_prime.FluxEstimator(atype, energy_or_number)

def test_b1a_same_as_idl(seasonal_flux_estimator, idl_call_results):
    """
    Test the value of b1p are
    identical for a given i_mlt,j_mlat
    """
    j_mlat, i_mlt = idl_call_results['j_mlat'], idl_call_results['i_mlt']
    b1a, b2a = idl_call_results['b1a'], idl_call_results['b2a']
    b1py = seasonal_flux_estimator.b1a[i_mlt, j_mlat]
    assert abs(b1py-b1a)<0.000001

def test_b1_same_as_idl(seasonal_flux_estimator, idl_call_results):
    """
    Test the value of b1p are
    identical for a given i_mlt,j_mlat
    """
    j_mlat, i_mlt = idl_call_results['j_mlat'], idl_call_results['i_mlt']
    b1p, b2p = idl_call_results['b1p'], idl_call_results['b2p']
    b1py = seasonal_flux_estimator.b1p[i_mlt, j_mlat]
    assert abs(b1py-b1p)<0.000001

def test_prob_estimate_same_as_idl(seasonal_flux_estimator, idl_call_results):
    """
    Check that prob_estimate method of SeasonalFluxEstimator
    produces the same result for a given
    bin (i_mlat,j_mlt), coupling strength (dF), season,
    auroral type (atype) and fluxtype (jtype)
    """
    j_mlat, i_mlt = idl_call_results['j_mlat'], idl_call_results['i_mlt']
    idl_dF = idl_call_results['dF']
    idl_prob = idl_call_results['prob']
    py_prob = seasonal_flux_estimator.prob_estimate(idl_dF, i_mlt, j_mlat)
    assert abs(py_prob-idl_prob)<0.000001

def test_flux_same_as_idl(seasonal_flux_estimator, idl_call_results):
    """
    Check the seasonal flux is the same as produced
    by the IDL code for the same season and auroral type
    for a given i_mlt,j_mlat
    """
    j_mlat, i_mlt = idl_call_results['j_mlat'], idl_call_results['i_mlt']
    idl_dF = idl_call_results['dF']
    py_flux = seasonal_flux_estimator.estimate_auroral_flux(idl_dF, i_mlt, j_mlat)
    idl_flux = idl_call_results['je']
    assert py_flux == idl_flux
