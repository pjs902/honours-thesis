from ssptools import evolve_mf_3 as emf3
import numpy as np
import pytest

from binaryshift import BinaryShift


@pytest.fixture
# TODO: figure out how to ignore the SSPTools warnings
def f():
    # config for ssptools
    m123 = [0.1, 0.5, 1.0, 100]  # Slope breakpoints for initial mass function
    a12 = [-0.468, -1.178, -2.117]  # Slopes for initial mass function
    nbin12 = [5, 5, 20]

    # Output times for the evolution
    tout = np.array([11000])

    # Integration settings
    N0 = 5e5  # Normalization of stars
    Ndot = (
        -0.0001
    )  # Regulates how low mass objects are depleted default -20, 0 for 47 Tuc
    tcc = 0  # Core collapse time
    NS_ret = 0.1  # Initial neutron star retention
    BH_ret_int = 1  # Initial Black Hole retention
    BH_ret_dyn = 0.00235  # Dynamical Black Hole retention
    FeHe = -0.7  # Metallicity
    # ignore the warnings from SSPTools
    f = emf3.evolve_mf(
        m123=m123,
        a12=a12,
        nbin12=nbin12,
        tout=tout,
        N0=N0,
        Ndot=Ndot,
        tcc=tcc,
        NS_ret=NS_ret,
        BH_ret_int=BH_ret_int,
        BH_ret_dyn=BH_ret_dyn,
        FeHe=FeHe,
        natal_kicks=True,
        vesc=100,
    )
    return f


@pytest.fixture
def mj_Mj(f):
    cs = f.Ns[-1] > 10 * f.Nmin
    cr = f.Nr[-1] > 10 * f.Nmin
    mj = np.r_[f.ms[-1][cs], f.mr[-1][cr]]
    Mj = np.r_[f.Ms[-1][cs], f.Mr[-1][cr]]
    dms = f.mes[-1][1:] - f.mes[-1][0:-1]
    nms = len(f.ms[-1][cs])
    return mj, Mj


@pytest.fixture
def binshift(mj_Mj, f):
    mj, Mj = mj_Mj
    bs = BinaryShift(mj=mj, Mj=Mj, MF=f, verbose=True)

    # make sure the lengths match
    with pytest.raises(ValueError):
        BinaryShift(mj=[2], Mj=[3, 2], MF=f, verbose=True)
    return bs


# Test mass conservation


def test_shift_q(mj_Mj, binshift):

    # do the shifting
    mj, Mj = mj_Mj
    mj_new, Mj_new = binshift.shift_q(fb=[0.1, 0.1, 0.1], q=[0.3, 0.5, 0.8])

    # check mass conservation
    assert np.isclose(np.sum(Mj), np.sum(Mj_new))

    # check number conservation
    Ntotal_initial = np.sum(Mj[binshift.MS_mask]/mj[binshift.MS_mask])
    Nj = Mj_new / mj_new

    Ntotal_shifted = 2*np.sum(Nj[binshift.bin_mask]) + np.sum(Nj[binshift.MS_mask_new])

    assert np.isclose(Ntotal_initial, Ntotal_shifted)


    # these are vad values for q, fb
    with pytest.raises(ValueError):
        binshift.shift_q(fb=[0.3, 0, 4, 0.3], q=[3, 0.5, 0.8])

    with pytest.raises(ValueError):
        binshift.shift_q(fb=[2, 0.5], q=[0.5, 0.8])


def test_shift_equal(mj_Mj, binshift):

    # do the shifting
    mj, Mj = mj_Mj
    mj_new, Mj_new = binshift.shift_equal(fb=0.3)

    # check mass conservation
    assert np.isclose(np.sum(Mj), np.sum(Mj_new))

    # check number conservation
    Ntotal_initial = np.sum(Mj[binshift.MS_mask]/mj[binshift.MS_mask])
    Nj = Mj_new / mj_new

    Ntotal_shifted = 2*np.sum(Nj[binshift.bin_mask]) + np.sum(Nj[binshift.MS_mask_new])

    assert np.isclose(Ntotal_initial, Ntotal_shifted)


    # bad value for fb
    with pytest.raises(ValueError):
        binshift.shift_equal(fb=3)


def test_shift_solar(mj_Mj, binshift):
    binshift.shift_solar()


def test_shift_kroupa(mj_Mj, binshift):
    binshift.shift_kroupa()


def test_shift_dump(binshift):
    binshift.dump()
