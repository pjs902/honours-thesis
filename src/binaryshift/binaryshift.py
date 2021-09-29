import numpy as np


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
        self.mf = MF

        # IFMR stuff
        self.nms = np.sum(np.array(self.mf.Ns) > self.mf.Nmin) - 1
        self.mWD_max = self.mf.IFMR.predict(self.mf.IFMR.wd_m_max)
        self.mBH_min = self.mf.mBH_min

        # Here are a bunch of masks for the bins
        self.MS_mask = np.ones_like(self.mj, dtype=bool)
        self.MS_mask[self.nms + 1 :] = False
        self.WD_mask = (self.mj <= self.mWD_max) & ~self.MS_mask
        self.NS_mask = (self.mj < self.mBH_min) & (self.mj > self.mWD_max)
        self.BH_mask = self.mj >= self.mBH_min

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
        print(f"{self.nms = }")

        print(f"{self.mWD_max = }")

        print(f"{self.mBH_min = }")

    def shift_equal(self, fb):
        """
        Shift mass to create binaries of equal mass, amount of mass shifted is determined by `fb`.
        """

        if fb > 1.0 or fb < 0:
            raise ValueError("fb must be between 0 and 1.")

        self.q = 1
        self.fb = fb

        # don't mess with original
        mj = self.mj.copy()
        Mj = self.Mj.copy()

        # loop through the MS mass bins
        for i in range(self.nms):

            if self.verbose:
                print()
                print(f"current mj: {mj[i]:.3f}")

            # get mass of companion
            companion_mass = mj[i] * self.q
            if self.verbose:
                print(f"{companion_mass = :.3f}")

            # mass of new bin
            mj_bin = mj[i] + companion_mass
            if self.verbose:
                print(f"new mass: {mj_bin:.3f} ")

            # total mass in binaries for this new bin
            Mj_bin = Mj[i] * self.fb
            if self.verbose:
                print(f"total mass in binaries {Mj_bin:.3f}")

            # add in the new mean mass bin
            mj = np.append(mj, mj_bin)

            # add in the new total mass bin
            Mj = np.append(Mj, Mj_bin)

            # remove the mass from the old total mass bin
            Mj[i] -= Mj_bin
        return mj, Mj

    def shift_q(self, fb, q):
        """
        Shift mass in to binaries with mass function `q`, amount of mass shifted determined by `fb`.
        (TODO) Eventually this will support lists of `fb` and `q`.
        """

        if q > 1.0 or q < 0:
            raise ValueError("q must be between 0 and 1.")

        if fb > 1.0 or fb < 0:
            raise ValueError("fb must be between 0 and 1.")

        self.fb = fb
        self.q = q

        # don't mess with original
        mj = self.mj.copy()
        Mj = self.Mj.copy()

        # loop through the MS mass bins
        for i in range(self.nms):
            if self.verbose:
                print()
                print(f"current mj: {mj[i]:.3f}")

            # get mass of companion
            companion_mass = mj[i] * self.q
            if self.verbose:
                print(f"{companion_mass = :.3f}")

            # if the companion is smaller than the lightest MS bin, just skip it
            if companion_mass < np.min(mj):
                if self.verbose:
                    print(
                        f"companion mass {companion_mass:.3f} smaller than {np.min(mj):.3f}, skipping"
                    )
                pass
            else:

                # find closest bin to companion mass
                companion_idx = np.argmin(np.abs(mj[: self.nms] - companion_mass))
                if self.verbose:
                    print(f"closest {companion_idx = }")
                # here change the mass of the companion to the mass of the closest bin (TODO: do we want this?)
                companion_mass = mj[companion_idx]
                if self.verbose:
                    print(f"closest {companion_mass = :.3f}")

                # mass of new bin
                binary_mj = mj[i] + companion_mass
                if self.verbose:
                    print(f"new mass: {binary_mj:.3f} ")

                # get total mass of new binary bin, will be (fb * mj) + (fb * companion mass bin)
                primary_Mj = self.fb * Mj[i]
                companion_Mj = self.fb * Mj[companion_idx]

                binary_Mj = primary_Mj + companion_Mj

                # add in new binary mean mass bin
                mj = np.append(mj, binary_mj)

                # add total mass to new binary bin
                Mj = np.append(Mj, binary_Mj)

                # remove mass from both primary, companion bins
                Mj[i] -= primary_Mj
                Mj[companion_idx] -= companion_Mj

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
