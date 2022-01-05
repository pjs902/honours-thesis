import numpy as np
from collections import namedtuple

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
        self.MS_mask_original = np.ones_like(self.mj, dtype=bool)
        self.MS_mask_original[self._nms + 1 :] = False
        self.WD_mask_original = (self.mj <= self._mWD_max) & ~self.MS_mask_original
        self.NS_mask_original = (self.mj < self._mBH_min) & (self.mj > self._mWD_max)
        self.BH_mask_original = self.mj >= self._mBH_min

        # minimum possible q value based on the mass bins
        self._q_min = np.min(self.mj[self.MS_mask_original]) / np.max(self.mj[self.MS_mask_original])

        self.verbose = verbose

        # Keep track of the rebinning
        self.previous_rebin = None


    def _shift_q(self, fb, q):
        """
        Shift mass in to binaries with mass fraction `q`, amount of mass shifted determined by `fb`.
        This is wrapped by `shift_solar` and `shift_flat` to make it easier to use for common cases.

        NOTE: this is not a public method, it is used internally by the other recipes and assumes
        that `fb` has already been converted to the internal `fb_ratio` quantity which we use for
        moving the mass around

        TODO: make wrapper which lets you use this nicely without worrying about converting between
        `fb` and `fb_ratio`
        """

        fb_arr = np.array(fb)
        q_arr = np.array(q)

        # check for invalid values
        if np.any(q_arr > 1.0) or np.any(q_arr < 0):
            raise ValueError("q must be between 0 and 1.")

        if np.any(fb_arr > 1.0) or np.any(fb_arr < 0):
            raise ValueError("fb must be between 0 and 1.")

        # don't mess with original
        mj = self.mj.copy()
        Mj = self.Mj.copy()
        Nj = Mj / mj
        Nj_shifted = Nj.copy()

        self.binary_components = []

        #####################################
        # loop through the binary mass ratios
        #####################################
        for fb, q in zip(fb_arr, q_arr):

            # loop through the MS mass bins
            for i in range(self._nms + 1):
                if self.verbose:
                    print()
                    print(f"current mj: {mj[i]:.3f}")

                ###############################################
                # Choosing the companion bin to the current bin
                ###############################################

                # get mass of companion
                companion_mass = mj[i] * q
                if self.verbose:
                    print(f"{companion_mass = :.3f}")

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
                                f"companion mass {companion_mass:.3f} larger than {np.max(mj[:self._nms+1]):.3f}, skipping companion star, original {q=}"
                            )
                        # go to next bin
                        continue

                ##################################
                # Selecting the corresponding bins
                ##################################

                # find closest bin to companion mass
                companion_idx = np.argmin(np.abs(mj[: self._nms + 1] - companion_mass))
                if self.verbose:
                    print(f"closest {companion_idx = }")

                # here change the mass of the companion to the mass of the closest bin
                companion_mass = mj[companion_idx]
                if self.verbose:
                    print(f"closest {companion_mass = :.3f}")

                ########################
                # Moving the mass around
                ########################

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

                # update binary components
                self.binary_components.append((mj[i], mj[companion_idx]))

                # remove N from both primary, companion bins
                # NOTE: This is set up so that if a bin gets below zero, we set it to zero and
                # correct the binary_Nj by the corresponding amount, this lets us go to higher
                # fb without breaking but it does mean that above ~45% we see more than 5%
                # errors in fb, which is not ideal but for realistic values of fb it's fine.
                Nj_shifted[i] -= binary_Nj
                if Nj_shifted[i] < 0:
                    Nj_shifted[-1] += Nj_shifted[i]
                    Nj_shifted[i] = 0
                    continue

                Nj_shifted[companion_idx] -= binary_Nj
                if Nj_shifted[companion_idx] < 0:
                    Nj_shifted[-1] += Nj_shifted[companion_idx]
                    Nj_shifted[companion_idx] = 0
                    continue

        # set the binary mask
        self.bin_mask = np.array([False] * len(mj))
        self.bin_mask[self._len_mj_init :] = True
        # add the extra Falses onto the ends of the other masks
        nbins_added = len(self.bin_mask) - len(self.MS_mask_original)
        self.MS_mask = np.append(self.MS_mask_original, [False] * nbins_added)
        self.WD_mask = np.append(self.WD_mask_original, [False] * nbins_added)
        self.NS_mask = np.append(self.NS_mask_original, [False] * nbins_added)
        self.BH_mask = np.append(self.BH_mask_original, [False] * nbins_added)

        Mj = Nj_shifted * mj
        # TODO: here check if any bins are empty (might want to do this depending on
        # if we use the original bins from SSPTools without removing the empty bins
        # first in order to use all of the bin edges)
        # cs = Nj_shifted > 10* self._mf.Nmin
        # mj = mj[cs]
        # Mj = Mj[cs]
        # print(f"{cs = }")
        # print(f"{self.BH_mask = }")

        # Here it might be nice to compute the "true fb" especially while we're troubleshooting
        self.fb_true = np.sum(Nj_shifted[self.bin_mask]) / (
            np.sum(Nj_shifted[self.bin_mask]) + np.sum(Nj_shifted[self.MS_mask])
        )

        # keep these around for rebinning and such
        self.Nj_shifted = Nj_shifted
        self.Mj_shifted = Mj
        self.mj_shifted = mj

        # lets also compute the q value for each binary bin
        qs = []
        for masses in self.binary_components:
            q = np.min(masses) / np.max(masses)
            qs.append(q)
        self.q_values = qs

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

        # switch to the internal quantity we use for fb
        self.fb_ratio = self.fb / (1.0 + self.fb)

        # here find the individual fb for each q value by adjusting the total fb by P(q)
        fb_arr = self.fb_ratio * p_q

        return self._shift_q(fb=fb_arr, q=q)

    def shift_equal(self, fb):
        """
        Shift mass to create binaries of equal mass, amount of mass shifted is determined by `fb`.
        """

        # validate fb
        self.fb = float(fb)
        if self.fb > 1.0 or self.fb < 0:
            raise ValueError("fb must be between 0 and 1.")

        self.fb_ratio = fb / (1.0 + fb)
        return self._shift_q([self.fb_ratio], [1.0])

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

        # switch to the internal quantity we use for fb
        self.fb_ratio = self.fb / (1.0 + self.fb)

        fb_arr = self.fb_ratio * p_q
        assert np.isclose(np.sum(p_q), 1.0)

        return self._shift_q(fb=fb_arr, q=q)

    def rebin(self, bins=15):
        """
        Down-sample the number of binary bins to a reasonable number for running models. Testing
        seems to indicate 15 bins is still fast enough to run models while keeping lots of
        resolution, but this can be adjusted to fit the use-case.
        """



        # first we should make sure that we've already done the shifting and
        # mj_shifted and bin_mask exist
        if not hasattr(self, "mj_shifted") or not hasattr(self, "bin_mask"):
            raise ValueError("Must shift before you can rebin!")

        # now check that we haven't already rebinned
        if self.previous_rebin is not None:
            if self.previous_rebin != bins:
                raise ValueError("Cannot rebin to a different number of bins!")
            if self.previous_rebin == bins:
                return self.mj_shifted, self.Mj_shifted

        # now we can do the rebinning

        # first calculate the bins
        _, binned = np.histogram(self.mj_shifted[self.bin_mask], bins=bins)
        # switch from bin edges to bin centers
        binned = np.array(
            [(binned[i] + binned[i + 1]) / 2 for i in range(len(binned) - 1)]
        )

        # then rebin the data

        # empty arrays to hold the new data
        new_Mj_binned = np.zeros_like(binned)
        new_Nj_binned = np.zeros_like(binned)

        # rebinned indexes for each binary mj
        bin_idxs = [
            (np.abs(binned - mj)).argmin() for mj in self.mj_shifted[self.bin_mask]
        ]

        # hold the q values for each bin

        # NamedTuple to hold the info for binary population in a combined bin
        BinaryPop = namedtuple("BinaryPop", ["mj", "q", "Mj"])

        rebinned_q_vals = [[] for _ in range(bins)]

        # loop over each binary bin
        for i in range(len(self.mj_shifted[self.bin_mask])):
            new_Mj_binned[bin_idxs[i]] += self.Mj_shifted[self.bin_mask][i]
            new_Nj_binned[bin_idxs[i]] += self.Nj_shifted[self.bin_mask][i]
            # keep track of the q values for each bin
            rebinned_q_vals[bin_idxs[i]].append(
                BinaryPop(
                    self.mj_shifted[self.bin_mask][i],
                    self.q_values[i],
                    self.Mj_shifted[self.bin_mask][i],
                )
            )

        # get new mean masses of rebinned binaries
        new_mj_binned = new_Mj_binned / new_Nj_binned

        # then replace all the old values
        self.q_values = rebinned_q_vals

        # remove old binaries, append new ones
        self.mj_shifted = np.append(self.mj_shifted[~self.bin_mask], new_mj_binned)
        self.Mj_shifted = np.append(self.Mj_shifted[~self.bin_mask], new_Mj_binned)

        # then remake the masks
        # remove old binaries from masks, append new ones
        binned_binaires_added = len(new_mj_binned)
        self.BH_mask = np.append(
            self.BH_mask[~self.bin_mask],
            np.zeros(binned_binaires_added, dtype=bool),
        )
        self.MS_mask = np.append(
            self.MS_mask[~self.bin_mask],
            np.zeros(binned_binaires_added, dtype=bool),
        )
        self.WD_mask = np.append(
            self.WD_mask[~self.bin_mask],
            np.zeros(binned_binaires_added, dtype=bool),
        )
        self.NS_mask = np.append(
            self.NS_mask[~self.bin_mask],
            np.zeros(binned_binaires_added, dtype=bool),
        )
        self.bin_mask = np.append(
            self.bin_mask[~self.bin_mask], np.ones(binned_binaires_added, dtype=bool)
        )

        # now we can set the already_rebinned flag
        self.previous_rebin = bins

        # return mj, Mj here to be consistent with shift_q
        return self.mj_shifted, self.Mj_shifted
