# Literature Review Planning

I'll probably have two overall sections, one discussing observations work and one discussing
modelling work. For each paper have the author, title and short summary or important points.

## Observations

- ### Ji+2015 - Binary Frequencies in a Sample of Globular Clusters. II. Sample Analysis and Comparison To Models

    Looks for binaries in Hubble data for clusters in Harris catalogue. Gives estimate on $fb$ for
    $q>0.5$ and for some way of estimating the smaller $q$ values too. Also gives a minimum observed
    value for $q$ which is nice. Further, gives results just within half-light radius and within the
    core. Looks at the radial distribution of binaries in clusters, this is something it might be
    nice to look at in our models.

- ### Fisher+2005 - What a local sample of spectroscopic binaries can tell us about the field binary population

    Looks at binaries (mostly spectroscopic I think) in the solar neighbourhood and gives some
    binary parameter distributions. We use their data for the solar $q$ value distribution in our
    modelling.

- ### Sollima+2007 - The fraction of binary systems in the core of 13 low-density Galactic globular clusters

    This one looks at binaries in low-density (young) clusters, generally finds much higher binary
    fractions, supporting the idea that high-density cluster loose many binaries from dynamical
    interactions. The numbers from this one probably aren't too useful for our clusters but nice to
    have for full overview of clusters.

- ### Milone+2012 - The ACS survey of Galactic globular clusters: XII. Photometric binaries along the main sequence

    Really nice discussion of photometric binaries, long and in depth with nice figures. It also
    gives more support to assuming a flat $q$ distribution in old, dense clusters. Probably will use
    this as the main support for discussing the photometric detection of binaries in globular
    clusters (which makes up most of the detections due to spectroscopy bring more costly).

- ### Giesers+2019 - A stellar census in globular clusters with MUSE: Binaries in NGC 3201

    Uses MUSE to look at binaries in NGC 3201, findings in line with what we expect mostly, binary
    fraction of $\sim6.76\%$ and a radial increase in binaries as you go to the core due to mass
    segregation. Here they are actually able to fit Keplerian orbits to $95$ stars with `The JOKER`
    which gives some insight into their properties.

- ### Lucatello+2015 - The incidence of binaries in globular cluster stellar populations

    Uses VLT FLAMES data to look for binaries in clusters, obviously observing far fewer stars than
    the MUSE data. Overall like $2\%$ binary fraction, but they mostly look at the differences in
    multiple populations.

- ### Albrow+2001 - The Frequency of Binary Stars in the Core of 47 Tucanae

    They use differential photometry to get time-series photometry for many stars and then look for
    variability in these stars, fairly uncommon method I think, would be good to mention it. They
    find core binary fraction of $\sim 13\%$ which seems pretty in line for other studies.

- ### Raghavan+2010 - A survey of stellar families: Multiplicity of solar-type stars

    Look at solar type stars, gives some info on their multiplicity and other parameters, maybe we
    can use something form a discussion of mass ratio distributions?

- ### Davis+2008 - Deep advanced camera for surveys imaging in the globular cluster NGC 6397: The binary fraction

    Uses main-sequence photometry to find binaries in NGC 6397, doesn't look like it has anything
    more than Milone+2012.

## Modelling

- ### Ivanova+2005 - The evolution of binary fractions in globular clusters

    Monte Carlo modelling of GCs with `FewBody` handling the binary interactions. Finds that present
    day binary fractions should be very low, upper bound of $5-10\%$. Goes over how binaries are
    ionized in GCs through dynamical interactions and also stellar evolution. Also gives some
    estimates for WDs which could be useful in the future.

- ### Sollima+2008 - The evolution of the binary population in globular clusters: A full analytical computation

    Presents an analytic model for binary evolution in clusters. Goes through all the possible
    interactions and how they're modelled. Supports low binary fractions in dense clusters and also
    looks at the radial distribution of binaries in clusters. Has a nice blurb about why equal-mass
    binaries are more common.

- ### Hurley+2007 - The Core Binary Fractions of Star Clusters from Realistic Simulations

    Uses N-body models to look at Binary fractions. This one finds that binary fractions basically
    remain unchanged from primordial, Sollima+2008 discusses why these numbers don't match with
    Ivanova, mostly due to different initial conditions and assumptions about the period
    distribution (if we truncate the period distribution at a higher value we will have more long
    period binaries, so the ionization process will be more efficient.). Initial conditions from
    Hurley are closer to the low-density clusters discussed in Sollima+2007 instead of the
    higher-density clusters in Ivanova+2005.

- ### Rodriguez+2021 - Modelling Dense Star Clusters in the Milky Way and Beyond with the Cluster Monte Carlo Code

    CMC paper. Can use for any info needed about monte-carlo models.

- ### Giersz+2013, Hypki+2013, MOCCA Code for Star Cluster Simulation

    Some of the MOCCA papers, should at least mention these when discussing modelling.

- ### Fregeau+2009 - EVOLUTION OF THE BINARY FRACTION IN DENSE STELLAR SYSTEMS

    Use a Monte Carlo code to look at the evolution of binary fractions in GCs. Discusses the
    differences between Hurley and Ivanova like Sollima+2008 does.

## Other

- ### Duchene+2013 - Stellar Multiplicity

    General review of multiple stars, might be useful for random stuff here and there.
