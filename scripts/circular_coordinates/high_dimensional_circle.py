"""Reproducing figure 6 in the paper"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import Normalize
from sklearn.decomposition import PCA

from lifting_cocycles.circular_coordinates import (
    compute_circular_coordinates_from_cocycle,
    compute_persistent_cocycle,
)
from lifting_cocycles.datasets import high_dim_circle
from lifting_cocycles.utils import plot_dgms

sns.set_theme("talk", "ticks", palette="tab10")

outdir = Path("figs")

SEED = 0
PRIME = 47
N = 1000
LIFT_SCALARS = (1, 3)


def main():
    t, X = high_dim_circle(N, 300)

    persistent_cocycle = compute_persistent_cocycle(
        X,
        prime=PRIME,
        persistence_threshold=20.0,
    )
    results = [
        compute_circular_coordinates_from_cocycle(
            persistent_cocycle,
            lift_scalar=scalar,
            return_scaled_coordinate=True,
        )
        for scalar in LIFT_SCALARS
    ]

    pca = PCA()
    X_pca = pca.fit_transform(X)

    # Persistence diagram on the left, centered beside the 2x2 grid.
    fig = plt.figure(figsize=(12, 7), layout="constrained")
    outer_grid = fig.add_gridspec(1, 2, width_ratios=(1, 2))
    diagram_grid = outer_grid[0, 0].subgridspec(
        3,
        1,
        height_ratios=(0.75, 2.5, 0.75),
    )
    coordinate_grid = outer_grid[0, 1].subgridspec(2, 2)

    diagram_ax = fig.add_subplot(diagram_grid[1, 0])
    axes = np.array(
        [
            [fig.add_subplot(coordinate_grid[row, column]) for column in range(2)]
            for row in range(2)
        ]
    )

    angle_norm = Normalize(0, 2 * np.pi)

    for i, result in enumerate(results):
        scatter = axes[i, 0].scatter(
            X_pca[:, 0],
            X_pca[:, 1],
            s=20,
            c=result.angles,
            cmap="winter",
            norm=angle_norm,
        )
        axes[i, 0].set(xticks=[], yticks=[], aspect="equal", alpha=0.75)

        axes[i, 1].scatter(
            t,
            result.angles,
            s=20,
            c=result.angles,
            cmap="winter",
            norm=angle_norm,
        )
        axes[i, 1].set(xticks=[], yticks=[], aspect="equal", alpha=0.75)

    colorbar = fig.colorbar(
        scatter, ax=axes, ticks=[0, np.pi, 2 * np.pi], aspect=40, pad=0.04
    )
    colorbar.ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])

    plot_dgms(diagram_ax, persistent_cocycle.diagrams)
    diagram_ax.scatter(
        persistent_cocycle.birth,
        persistent_cocycle.death,
        s=300,
        facecolors="none",
        edgecolors="black",
        linewidths=2,
        zorder=3,
    )
    diagram_ax.set(xticks=[], yticks=[])
    plt.savefig(outdir / Path("high-dim-circle.pdf"))
    plt.show()


if __name__ == "__main__":
    main()
