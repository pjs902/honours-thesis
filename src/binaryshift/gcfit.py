from importlib import resources

import astropy.units as u
import numpy as np
import pandas as pd
import scipy as sp

from .binaryshift import BinaryShift

"""
Module which contains all the functionality specifc to interacting with `GCFit`.
"""


def from_gcfit(model):
    """
    Create a `BinaryShift` object from a `GCFit` `Model` object. Sets the correct flags to use units properly.
    """

    _mf = model._mf

    # Set bins that should be empty to empty
    cs = _mf.Ns[-1] > 10 * _mf.Nmin
    ms, Ms = _mf.ms[-1][cs], _mf.Ms[-1][cs]

    cr = _mf.Nr[-1] > 10 * _mf.Nmin
    mr, Mr = _mf.mr[-1][cr], _mf.Mr[-1][cr]

    # Collect mean mass and total mass bins
    mj = np.r_[ms, mr]
    Mj = np.r_[Ms, Mr]

    # make binshift instance
    binshift = BinaryShift(mj=mj, Mj=Mj, MF=_mf, GCFit=True, model=model)

    return binshift


def get_isochrone(model):
    """
    Get the isochrone closest to the Fe/H of the model
    """

    # get Fe/H from model
    feh = model._mf.FeHe

    # get isochrone list
    # this ugly stuff is python 3.8 to 3.9 compatability
    try:
        with resources.files("binaryshift") / "resources" as path:
            isochrones = list(path.glob("*.dat"))
    except AttributeError:
        with resources.path('binaryshift', 'resources') as path:
            isochrones = list(path.glob("*.dat"))

    fehs = [float(str(i).split("FEH=")[1].split(".dat")[0]) for i in isochrones]

    # get isochrone closest to Fe/H
    best = isochrones[np.abs(np.array(fehs) - feh).argmin()]
    return pd.read_csv(best)


def flatten(t):
    """
    Flatten a nested list
    https://stackoverflow.com/a/952952/13577036
    """
    return [item for sublist in t for item in sublist]


def rescale_densities(model):

    # make copy of density profile, Mj
    rescaled_Sigmaj = model.Sigmaj.copy()
    scaled_Mj = model.Mj.copy()

    # get binaries
    binaries = flatten(model._binshift.q_values)

    # get isochrone
    isochrone = get_isochrone(model)

    # loop through each binary population
    for binary in binaries:

        # find "observed" mass
        observed_mass = (
            get_observed_mass(isochrone=isochrone, mj=binary.mj, q=binary.q) << u.Msun
        )

        # find closest mass bin
        closest_idx = np.argmin(np.abs(model.mj[model._single_mask] - observed_mass))

        # scale the corresponding density bin

        # find scale factor
        scale_factor = (scaled_Mj[closest_idx] + (binary.Mj * u.Msun)) / scaled_Mj[
            closest_idx
        ]
        # apply scale
        scaled_Mj[closest_idx] *= scale_factor
        rescaled_Sigmaj[closest_idx] *= scale_factor

        # add rescaled density profiles to model
        model.rescaled_Sigmaj = rescaled_Sigmaj


def get_observed_mass(isochrone, mj, q):

    # TODO: this hardcoded stuff is not great, find a better way to do this
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
