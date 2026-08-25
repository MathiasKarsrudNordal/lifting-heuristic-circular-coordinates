"""Computing circular coordinates of the circle"""

import matplotlib.pyplot as plt
import seaborn as sns
from icecream.icecream import ic
from ripser import ripser
from scipy.spatial.distance import pdist, squareform

from lifting_cocycles.circular_coordinates import compute_circular_coordinates
from lifting_cocycles.datasets import circle
from lifting_cocycles.utils import plot_dgms

sns.set_theme("talk", "ticks", palette="tab10")

SEED = 0
PRIME = 47
N = 500


def main():

    X = circle(N, seed=SEED, noise_std=0.1)

    dist = squareform(pdist(X[:, :2], metric="euclidean"))
    rips = ripser(dist, distance_matrix=True, coeff=PRIME, do_cocycles=False)
    dgms = rips["dgms"]

    result = compute_circular_coordinates(
        X[:, :2],
        prime=PRIME,
        lift_scalar=5,
    )

    angles = result.angles

    fig, axes = plt.subplots(1, 3)

    scatter_true = axes[0].scatter(
        *X[:, :2].T,
        c=X[:, 2],
        cmap="winter",
        s=5,
        alpha=0.8,
    )
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_aspect("equal")

    scatter_infer = axes[1].scatter(
        *X[:, :2].T,
        c=angles,
        cmap="winter",
        s=5,
        alpha=0.8,
    )
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].set_aspect("equal")

    scatter_corr = axes[2].scatter(
        X[:, 2].T,
        angles,
        s=5,
        alpha=0.8,
    )
    axes[2].set_xticks([])
    axes[2].set_yticks([])
    axes[2].set_aspect("equal")

    plt.tight_layout()
    plt.show()

    # plotting persistent diagram
    _, ax = plt.subplots()

    _ = plot_dgms(ax, dgms)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
