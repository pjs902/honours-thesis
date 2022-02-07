"""
Module which contains all the functionality specifc to interacting with `GCFit`.
"""

from .binaryshift import BinaryShift
import numpy as np


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
    pass


def rescale_densities(binshift, model):
    """
    TODO: Functionality to rescale density profiles for MF fitting
    """
