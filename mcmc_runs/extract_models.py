from multiprocessing import Pool
import emcee
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import corner
from tqdm import tqdm
import fitter
from tqdm.contrib.concurrent import process_map

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


# Config
DISCARD = 1500
N_MODELS = 1024
FB = 0.02
N_CORES = 4
filename = "NGC0104_sampler.hdf"
obs = fitter.Observations("NGC0104")


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

    # Read file
    reader = emcee.backends.HDFBackend(filename, read_only=True)
    try:
        tau = reader.get_autocorr_time()
        print(tau)
    except:
        pass

    ndim = reader.shape[1]

    flat_samples = reader.get_chain(discard=DISCARD, thin=1, flat=True)

    # Plot Corner plot
    print("Plotting corner plot")
    fig = corner.corner(
        flat_samples, labels=labels, use_math_text=True, show_titles=True
    )
    plt.plot()
    fig.tight_layout()
    plt.savefig("corner.png", dpi=300)

    # plot chains
    print("Plotting chains")
    fig, axes = plt.subplots(ndim, figsize=(10, 20), sharex=True)
    samples = reader.get_chain()
    for i in range(len(labels)):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.plot()
    plt.savefig("walkers.png", dpi=300)

    # compute the models in parallel

    print("Computing models")
    models = process_map(
        get_model,
        tqdm(flat_samples[-N_MODELS:]),
        max_workers=N_CORES,
    )

    print("Saving models")
    np.save("models", np.array(models, dtype=object))


def get_model(theta):

    # compute GCfit Model
    model = fitter.Model(theta=theta, observations=obs, binary_fraction=FB)

    return model


if __name__ == "__main__":
    main()
