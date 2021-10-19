import numpy as np
import warnings

__all__ = ["BinaryShift"]


class BinaryShift:
    def __init__(self, mj, Mj, MF, verbose=False):
        """
        Initialize an instance of `BinaryShift` with the provided mass bins and MF, (TODO)
        the IFMR will (eventually) be used to label the object type of each bin.
        """

        if len(mj) != len(Mj):
            raise ValueError("mj and Mj must have the same length.")

        self.mj = mj
        self.Mj = Mj
        self._mf = MF

        # IFMR stuff
        self._nms = np.sum(np.array(self._mf.Ns) > self._mf.Nmin) - 1
        self._mWD_max = self._mf.IFMR.predict(self._mf.IFMR.wd_m_max)
        self._mBH_min = self._mf.mBH_min

        # Use this to make the binary mask
        self._len_mj_init = len(mj)

        # Here are a bunch of masks for the bins
        self.MS_mask = np.ones_like(self.mj, dtype=bool)
        self.MS_mask[self._nms + 1 :] = False
        self.WD_mask = (self.mj <= self._mWD_max) & ~self.MS_mask
        self.NS_mask = (self.mj < self._mBH_min) & (self.mj > self._mWD_max)
        self.BH_mask = self.mj >= self._mBH_min

        self.verbose = verbose

    def dump(self):
        """
        Dump current config of `BinaryShift`.
        """
        # Mass bins
        print(f"{self.mj = }")

        print(f"{self.Mj = }")

        # These ones might be set
        try:
            print(f"{self.q = }")
        except:
            pass

        try:
            print(f"{self.fb = }")
        except:
            pass

        # IFMR Stuff
        print(f"{self._nms = }")

        print(f"{self._mWD_max = }")

        print(f"{self._mBH_min = }")

        # TODO: dump masks?

    def shift_equal(self, fb):
        """
        Shift mass to create binaries of equal mass, amount of mass shifted is determined by `fb`.
        """

        return self.shift_q([fb], [1.0])

    def shift_q_OLD(self, fb, q):
        """
        Shift mass in to binaries with mass fraction `q`, amount of mass shifted determined by `fb`.
        NOTE: This is the old version that (sort of) used mass fraction, this is broken.
        """
        warnings.warn("This is the old broken method", DeprecationWarning)
        self.fb = np.array(fb)
        self.q = np.array(q)

        # check for invalid values
        if np.any(self.q > 1.0) or np.any(self.q < 0):
            raise ValueError("q must be between 0 and 1.")

        if np.any(self.fb > 1.0) or np.any(self.fb < 0):
            raise ValueError("fb must be between 0 and 1.")

        # don't mess with original
        mj = self.mj.copy()
        Mj = self.Mj.copy()

        # loop through the binary mass ratios
        for fb, q in zip(self.fb, self.q):

            # loop through the MS mass bins
            for i in range(self._nms + 1):
                if self.verbose:
                    print()
                    print(f"current mj: {mj[i]:.3f}")

                # get mass of companion
                companion_mass = mj[i] * q
                if self.verbose:
                    print(f"{companion_mass = :.3f}")

                # if the companion is much smaller than the lightest MS bin, just skip it
                if companion_mass < np.min(mj) and (
                    np.abs(companion_mass - np.min(mj)) > 0.025
                ):
                    if self.verbose:
                        print(
                            f"companion mass {companion_mass:.3f} smaller than {np.min(mj):.3f}, skipping"
                        )
                    pass
                else:

                    # find closest bin to companion mass
                    companion_idx = np.argmin(np.abs(mj[: self._nms] - companion_mass))
                    if self.verbose:
                        print(f"closest {companion_idx = }")
                    # here change the mass of the companion to the mass of the closest bin
                    companion_mass = mj[companion_idx]
                    if self.verbose:
                        print(f"closest {companion_mass = :.3f}")

                    # mass of new bin
                    binary_mj = mj[i] + companion_mass
                    if self.verbose:
                        print(f"new mass: {binary_mj:.3f} ")

                    # get total mass of new binary bin, will be (fb * mj) + (fb * companion mass bin)
                    # here use self.Mj so we take from the unchanged Mj
                    primary_Mj = fb * self.Mj[i]
                    companion_Mj = fb * self.Mj[companion_idx]

                    binary_Mj = primary_Mj + companion_Mj

                    # add in new binary mean mass bin
                    mj = np.append(mj, binary_mj)

                    # add total mass to new binary bin
                    Mj = np.append(Mj, binary_Mj)

                    # remove mass from both primary, companion bins
                    Mj[i] -= primary_Mj
                    Mj[companion_idx] -= companion_Mj

        # set the binary mask
        self.bin_mask = np.array([False] * len(mj))
        self.bin_mask[self._len_mj_init :] = True
        # add the extra Falses onto the ends of the other masks
        nbins_added = len(self.bin_mask) - len(self.MS_mask)
        self.MS_mask_new = np.append(self.MS_mask, [False] * nbins_added)
        self.WD_mask_new = np.append(self.WD_mask, [False] * nbins_added)
        self.NS_mask_new = np.append(self.NS_mask, [False] * nbins_added)
        self.BH_mask_new = np.append(self.BH_mask, [False] * nbins_added)

        return mj, Mj

    def shift_q(self, fb, q):
        """
        Shift mass in to binaries with mass fraction `q`, amount of mass shifted determined by `fb`.
        NOTE: This is the new version that uses the number of stars/binaries instead of their masses.
        """

        self.fb = np.array(fb)
        self.q = np.array(q)

        # check for invalid values
        if np.any(self.q > 1.0) or np.any(self.q < 0):
            raise ValueError("q must be between 0 and 1.")

        if np.any(self.fb > 1.0) or np.any(self.fb < 0):
            raise ValueError("fb must be between 0 and 1.")

        # don't mess with original
        mj = self.mj.copy()
        Mj = self.Mj.copy()
        Nj = Mj / mj
        Nj_shifted = Nj.copy()

        # loop through the binary mass ratios
        for _fb, q in zip(self.fb, self.q):

            # NOTE: Here's a conversion from one way of counting to the other, works great for
            # larger fb, seems slightly off for smaller fb?
            fb = _fb / (1.0 + _fb)

            # loop through the MS mass bins
            for i in range(self._nms + 1):
                if self.verbose:
                    print()
                    print(f"current mj: {mj[i]:.3f}")

                # get mass of companion
                companion_mass = mj[i] * q
                if self.verbose:
                    print(f"{companion_mass = :.3f}")

                # if the companion is much smaller than the lightest MS bin, just skip it
                if companion_mass < np.min(mj) and (
                    np.abs(companion_mass - np.min(mj)) > 0.025
                ):
                    if self.verbose:
                        print(
                            f"companion mass {companion_mass:.3f} smaller than {np.min(mj):.3f}, skipping"
                        )
                    pass
                else:

                    # find closest bin to companion mass
                    companion_idx = np.argmin(np.abs(mj[: self._nms] - companion_mass))
                    if self.verbose:
                        print(f"closest {companion_idx = }")
                    # here change the mass of the companion to the mass of the closest bin
                    companion_mass = mj[companion_idx]
                    if self.verbose:
                        print(f"closest {companion_mass = :.3f}")

                    # mean mass of new bin
                    binary_mj = mj[i] + companion_mass
                    if self.verbose:
                        print(f"new mass: {binary_mj:.3f} ")

                    # get total N of new binary bin
                    binary_Nj = Nj[i] * fb  # / 2
                    if self.verbose:
                        print(f"current bin N: {Nj[i]:.3f} ")
                        print(f"binary N: {binary_Nj:.3f} ")

                    # add in new binary mean mass bin
                    mj = np.append(mj, binary_mj)

                    # add total N to new binary bin
                    Nj_shifted = np.append(Nj_shifted, binary_Nj)

                    # remove N from both primary, companion bins
                    Nj_shifted[i] -= binary_Nj

                    if Nj_shifted[i] < 0:
                        raise ValueError(
                            f"Value of fb is too high: bin {i} {mj[i] = } went negative"
                        )

                    Nj_shifted[companion_idx] -= binary_Nj

                    if Nj_shifted[companion_idx] < 0:
                        raise ValueError(
                            f"Value of fb is too high: bin {i} {mj[i] = } went negative"
                        )

        # set the binary mask
        self.bin_mask = np.array([False] * len(mj))
        self.bin_mask[self._len_mj_init :] = True
        # add the extra Falses onto the ends of the other masks
        nbins_added = len(self.bin_mask) - len(self.MS_mask)
        self.MS_mask_new = np.append(self.MS_mask, [False] * nbins_added)
        self.WD_mask_new = np.append(self.WD_mask, [False] * nbins_added)
        self.NS_mask_new = np.append(self.NS_mask, [False] * nbins_added)
        self.BH_mask_new = np.append(self.BH_mask, [False] * nbins_added)

        Mj = Nj_shifted * mj
        # TODO: here check if any bins are empty

        return mj, Mj

    def shift_kroupa(self):
        """
        (TODO) Shift mass according to `fb` and `q` determined by random draw from Kroupa IMF.
        """
        pass

    def shift_solar(self):
        """
        (TODO) Shift mass according to `fb` and `q` in the solar neighborhood.
        """
        pass
