import random
from importlib import resources

import numpy as np
import pandas as pd
import pytest

import binaryshift
from binaryshift import BinaryShift

try:
    from fitter import Model, Observations

    GCFIT = True
except ImportError:
    GCFIT = False


@pytest.fixture
def obs():
    return Observations("NGC0104")


@pytest.fixture
def model(obs):
    theta = [
        6.62,
        0.88,
        6.82,
        1.33,
        1.03,
        0.37,
        0.01,
        3.49,
        0.47,
        1.18,
        2.15,
        0.13,
        4.42,
    ]

    return Model(theta=theta, observations=obs)


@pytest.fixture
def binshift(model):
    return binaryshift.gcfit.from_gcfit(model)


# test that the isochrones are included and we can read them
def test_isochrones():
    with resources.files("binaryshift") / "resources" as path:
        isochrones = list(path.glob("*.dat"))
    test_isochrone = random.choice(isochrones)
    pd.read_csv(test_isochrone)


@pytest.fixture
def isochrone(model):
    return binaryshift.gcfit.get_isochrone(model)


def test_observed_mass(isochrone):
    # hardcoded values from testing
    observed = binaryshift.gcfit.get_observed_mass(isochrone=isochrone, mj=1, q=0.5)
    assert np.isclose(observed, 0.6735)


def test_rescale_densities(model):
    # rescale
    binaryshift.gcfit.rescale_densities(model)
    assert np.sum(model.rescaled_Sigmaj) >= 0
