import numpy as np
import pytest
import warnings

from binaryshift import BinaryShift
from binaryshift.gcfit import from_gcfit

from fitter import Observations
from fitter import Model


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


def test_from_gcfit(binshift):
    pass
