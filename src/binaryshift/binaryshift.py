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

        # minimum possible q value based on the mass bins
        self._q_min = np.min(self.mj[self.MS_mask]) / np.max(self.mj[self.MS_mask])

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

    def shift_q(self, fb, q):
        """
        Shift mass in to binaries with mass fraction `q`, amount of mass shifted determined by `fb`.
        This is wrapped by `shift_solar` and `shift_flat` to make it easier to use for common cases.
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

                # TODO: Here what we actually want to do is to use these low mass bins as
                # companions, so when we reach a point where there are no availible low mass bin for
                # companions we should instead look for a high mass primary star that still
                # satisfies the q value.

                # TODO: So this is sort of working now, but it's not quite right. fb doesn't match
                # the requested value anymore expect for 30%, check to see if its still taking away
                # mass and stars properly

                # TODO: refactor this to be more readable, maybe catch the error somewhere?
                if companion_mass < np.min(mj[: self._nms + 1]) and (
                    np.abs(companion_mass - np.min(mj[: self._nms + 1])) > 0.025
                ):
                    if self.verbose:
                        print(
                            f"companion mass {companion_mass:.3f} smaller than {np.min(mj):.3f}, switching to companion star"
                        )
                    new_q = 1.0 / q
                    companion_mass = mj[i] * new_q
                    if self.verbose:
                        print(f"new (primary) {companion_mass = :.3f}")
                        print(f"{new_q = :.3f}")
                    if companion_mass > np.max(mj[: self._nms + 1]) and (
                        np.abs(companion_mass - np.max(mj[: self._nms + 1])) > 0.025
                    ):
                        if self.verbose:
                            print(
                                f"companion mass {companion_mass:.3f} larger than {np.max(mj[:self._nms+1]):.3f}, skipping companion star"
                            )
                        # go to next bin
                        continue

                # find closest bin to companion mass
                companion_idx = np.argmin(np.abs(mj[: self._nms + 1] - companion_mass))
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
                        f"Value of {fb=} is too high: bin {i} {mj[i] = } went negative"
                    )
                Nj_shifted[companion_idx] -= binary_Nj
                if Nj_shifted[companion_idx] < 0:
                    raise ValueError(
                        f"Value of {fb=} is too high: bin {i} {mj[i] = } went negative"
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
        # TODO: here check if any bins are empty (doesn't really seem to be needed)
        # cs = Nj_shifted > 10* self._mf.Nmin
        # mj = mj[cs]
        # Mj = Mj[cs]
        # print(f"{cs = }")
        # print(f"{self.BH_mask = }")

        # Here it might be nice to compute the "true fb" especially while we're troubleshooting
        self.fb_true = np.sum(Nj_shifted[self.bin_mask]) / (
            np.sum(Nj_shifted[self.bin_mask]) + np.sum(Nj_shifted[self.MS_mask_new])
        )

        return mj, Mj

    def shift_solar(self, fb):
        """
        Shift mass according to `fb` and `q` in the solar neighborhood.
        Values for solar binaries from Fisher et al. (2005) (10.1111/j.1365-2966.2005.09193.x)
        """

        # validate fb
        self.fb = float(fb)
        if self.fb > 1.0 or self.fb < 0:
            raise ValueError("fb must be between 0 and 1.")

        # frequencies from Table 3 of Fisher et al. (2005) (10.1111/j.1365-2966.2005.09193.x)
        freqs = np.array([29, 29, 30, 32, 31, 32, 36, 45, 27, 76])

        # get P(q) for each q value (we could hardcode this I guess)
        p_q = freqs / np.sum(freqs)

        # full list of q values
        q = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

        # Truncate q distribution based on smallest possible q value
        # get the total fb from bad q values
        q_mask = q > self._q_min
        pq_disallowed = np.sum(p_q[~q_mask])

        # remove the values from the q and p(q) distributions
        q = q[q_mask]
        p_q = p_q[q_mask]

        # add the removed probabilities to the p(q) distribution
        extra_pq = pq_disallowed / len(q)
        p_q = p_q + extra_pq
        assert np.isclose(np.sum(p_q), 1.0)

        # here find the individual fb for each q value by adjusting the total fb by P(q)
        fb = fb * p_q

        return self.shift_q(fb=fb, q=q)

    def shift_equal(self, fb):
        """
        Shift mass to create binaries of equal mass, amount of mass shifted is determined by `fb`.
        """

        return self.shift_q([fb], [1.0])

    def shift_flat(self, fb):
        """
        Shift mass according to a flat `q` distribution.
        """

        # validate fb
        self.fb = float(fb)
        if self.fb > 1.0 or self.fb < 0:
            raise ValueError("fb must be between 0 and 1.")

        # full list of q values
        q = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

        # Truncate q distribution based on smallest possible q value
        q = q[q > self._q_min]

        p_q = np.ones_like(q) / len(q)
        # here find the individual fb for each q value by adjusting the total fb by P(q)
        fb = self.fb * p_q
        assert np.isclose(np.sum(p_q), 1.0)

        return self.shift_q(fb=fb, q=q)
