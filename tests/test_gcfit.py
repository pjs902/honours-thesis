import random
import warnings
from importlib import resources

import numpy as np
import pandas as pd
import pytest

from binaryshift import BinaryShift
from binaryshift.gcfit import from_gcfit

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
    return from_gcfit(model)


@pytest.mark.skipif(GCFIT == False, reason="GCFit is not installed")
def test_from_gcfit(binshift):
    pass


# test that the isochrones are included and we can read them
def test_isochrones():
    with resources.files("binaryshift") / "resources" as path:
        isochrones = list(path.glob("*.dat"))
    test_isochrone = random.choice(isochrones)
    pd.read_csv(test_isochrone)
