"""generating trefoil knot"""

import matplotlib.pyplot as plt
import seaborn as sns

from lifting_cocycles.datasets import trefoil_knot

sns.set_theme("talk", "ticks", palette="tab10")


def main():
    seed = 0
    N = 2500

    X = trefoil_knot(N, seed=seed, noise_std=0.1)

    # 2d plot
    fig = plt.figure()
    ax = fig.add_subplot()

    scatter = ax.scatter(
        *X[:, :2].T,
        c=X[:, 3],
        cmap="winter",
        s=5,
        alpha=0.8,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")

    cbar = fig.colorbar(
        scatter,
        ax=ax,
        fraction=0.05,
        pad=0.04,
    )
    cbar.set_ticks([])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
