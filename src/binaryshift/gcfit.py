from importlib import resources

import numpy as np
import pandas as pd
import scipy as sp

"""
Module which contains all the functionality specifc to interacting with `GCFit`.
"""

from .binaryshift import BinaryShift


def from_gcfit(model):
    """
    Create a `BinaryShift` object from a `GCFit` `Model` object. Sets the correct flags to use units properly.
    """

    # make binshift instance
    binshift = BinaryShift(
        mj=model.mj.value, Mj=model.Mj.value, MF=model._mf, GCFit=True, model=model
    )

    return binshift


def update_masks(binshift, model):
    """
    TODO: Remake GCFit masks, specific BH/NS quantities
    """

    # TODO: So here we need to update the masks that GCFit uses internally, also need to update
    # anything else similar as well as add in our own binary masks

    pass


def get_isochrone(model):
    """
    Get the isochrone closest to the Fe/H of the model
    """

    # get Fe/H from model
    feh = model._mf.FeHe

    # get isochrone list
    with resources.files("binaryshift") / "resources" as path:
        isochrones = list(path.glob("*.dat"))

    fehs = [float(str(i).split("FEH=")[1].split(".dat")[0]) for i in isochrones]

    # get isochrone closest to Fe/H
    best = isochrones[np.abs(np.array(fehs) - feh).argmin()]
    return pd.read_csv(best, engine="pyarrow")


def rescale_densities(binshift, model):
    """
    TODO: Functionality to rescale density profiles for MF fitting
    """

    # first get isochrone
    isochrone = get_isochrone(model)

    # create copies of the original profiles
    rhoj_rescaled = model.rhoj.copy()

    # loop over each binary component in each binary bin, and rescale the matching density profile

    # add the rescaled profiles to the model
    model.rhoj_rescaled = rhoj_rescaled


def get_observed_mass(isochrone, mj, q):

    MS_isochrone = isochrone[0 : int(1463 / 7)]

    mass_to_lum = sp.interpolate.InterpolatedUnivariateSpline(
        x=MS_isochrone.star_mass,
        y=(10 ** (MS_isochrone.WFC3_UVIS_F814W / -2.5)),
        k=3,
        ext=2,
    )
    lum_to_mass = sp.interpolate.InterpolatedUnivariateSpline(
        x=(10 ** (MS_isochrone.WFC3_UVIS_F814W / -2.5)),
        y=MS_isochrone.star_mass,
        k=3,
        ext=2,
    )

    # first calculate the individual masses
    mb = mj / (1 + q)
    ma = mj - mb
    # then get observed mass
    lum1 = mass_to_lum(ma)
    lum2 = mass_to_lum(mb)

    observed_mass = lum_to_mass(lum1 + lum2)
    return float(observed_mass)
