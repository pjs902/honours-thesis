from traceback import print_tb
import emcee
import matplotlib.pyplot as plt
import seaborn
import numpy as np
from tqdm import tqdm
import fitter
import seaborn as sns


# Some plotting config
sns.set(
    context="notebook",
    # style="ticks",
    style="darkgrid",
    font="Times New Roman",
    font_scale=1.75,
)


plt.rcParams.update({"text.usetex": True})


sns.color_palette("mako", as_cmap=True)
plt.rcParams["figure.figsize"] = (12, 10)
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
# plt.rcParams["xtick.top"] = True
# plt.rcParams["ytick.right"] = True
# plt.rcParams["xtick.bottom"] = True
# plt.rcParams["ytick.left"] = True
plt.rcParams["mathtext.fontset"] = "cm"


labels = [
    r"$W_{0}$",
    r"$M/10^6 M_{\odot}$",
    r"$r_h / pc$",
    r"$ log \ r_a / pc$",
    r"$g$",
    r"$\delta$",
    r"$s^2$",
    r"$F$",
    r"$\alpha_1$",
    r"$\alpha_2$",
    r"$\alpha_3$",
    r"$BH_{ret}$",
    r"$d$",
]


def main():

    models = np.load("models.npy", allow_pickle=True)

    # extract the BH info
    BHs = []
    N_BHs = []
    mj_BHs = []
    mj_BHs_weights = []
    for model in tqdm(models):
        total_BH_mass = np.sum(model.BH_Mj.value)
        BHs.append(total_BH_mass)
        total_BH_N = np.sum(model.BH_Nj)
        N_BHs.append(total_BH_N)

        mj_BHs.append(model.BH_mj.value)
        mj_BHs_weights.append(model.BH_Mj)

    # do the plots
    plot_BHs(BHs)
    plot_N_BHs(N_BHs)

    plot_mjs_BHs(mj_BHs, mj_BHs_weights)


def plot_BHs(BHs):
    plt.figure()
    plt.xlabel(r"Mass in BHs [$M_{\odot}$]")
    plt.ylabel("Probability density")

    sns.histplot(x=BHs, kde=False, stat="density")

    vals = np.percentile(BHs, [16, 50, 84])
    q = np.diff(vals)
    # plt.axvline(vals[1], ls="-")
    # plt.axvline(vals[0], ls="--")
    # plt.axvline(vals[2], ls="--")

    txt = "BH Mass = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(vals[1], q[0], q[1])
    print(
        r"BH Mass [Msol]"
        + " ="
        + txt.split("=")[1].split("_")[0]
        + " (+"
        + txt.split("^")[1].split("{")[1].split("}")[0]
        + " "
        + txt.split("_")[1].split("^")[0].split("{")[1].split("}")[0]
        + ")"
    )
    print(
        f"BH mass percentiles:\n 95: {np.percentile(BHs,95)} \n 99: {np.percentile(BHs,99)}"
    )

    # plt.axvline(np.percentile(BHs, 95), ls="--")

    plt.tight_layout()
    plt.plot()
    plt.tight_layout()
    plt.savefig("BH_contours.png", dpi=300)


def plot_mjs_BHs(mj_BHs, mj_BHs_weights):

    plt.figure()
    plt.xlabel(r"Black Hole Mass")
    plt.ylabel("Frequency")
    mj_BHs = np.hstack(mj_BHs).flatten()
    mj_BHs_weights = np.hstack(mj_BHs_weights).flatten()

    sns.histplot(
        x=mj_BHs, weights=mj_BHs_weights, kde=False, stat="frequency", discrete=True
    )

    plt.tight_layout()
    plt.plot()
    plt.tight_layout()
    plt.savefig("BH_mjs_contours.png", dpi=300)


def plot_N_BHs(N_BHs):
    plt.figure()
    plt.xlabel(r"Number of BHs")
    plt.ylabel("Probability density")

    sns.histplot(x=N_BHs, kde=False, stat="density")

    print(
        f"BH Number percentiles:\n 95: {np.percentile(N_BHs,95)} \n 99: {np.percentile(N_BHs,99)}"
    )
    vals = np.percentile(N_BHs, [16, 50, 84])
    q = np.diff(vals)
    txt = "BH Number = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(vals[1], q[0], q[1])
    print(
        r"BH Number "
        + " ="
        + txt.split("=")[1].split("_")[0]
        + " (+"
        + txt.split("^")[1].split("{")[1].split("}")[0]
        + " "
        + txt.split("_")[1].split("^")[0].split("{")[1].split("}")[0]
        + ")"
    )

    plt.plot()
    plt.tight_layout()
    plt.savefig("BH_N_contours.png", dpi=300)


if __name__ == "__main__":
    main()
