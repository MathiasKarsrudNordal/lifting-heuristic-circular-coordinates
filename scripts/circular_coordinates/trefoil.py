"""Reproducing figure 2 in the paper"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from lifting_cocycles.circular_coordinates import compute_circular_coordinates
from lifting_cocycles.datasets import torus_knot
from lifting_cocycles.utils import plot_dgms

sns.set_theme("talk", "ticks", palette="tab10")

outdir = Path("figs")

SEED = 0
PRIME = 47
N = 2500
NOISE_VARIANCE = np.sqrt(0.1)
LIFT_SCALARS = [1, 2]


def main():

    X = torus_knot(
        N,
        seed=SEED,
        p=2,
        q=3,
        R=5,
        r=2,
        noise_std=NOISE_VARIANCE,
    )

    results = [
        compute_circular_coordinates(
            X[:, :3],
            prime=PRIME,
            coordinate_fraction=0.75,
            lift_scalar=scalar,
        )
        for scalar in LIFT_SCALARS
    ]

    fig, axes = plt.subplots(1, len(results), figsize=(12, 6))

    for ax, scalar, result in zip(
        axes,
        LIFT_SCALARS,
        results,
        strict=True,
    ):
        ax.scatter(
            *X[:, :2].T,
            c=result.angles,
            cmap="winter",
            s=4,
            alpha=0.8,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(outdir / Path("trefoil-scaling.pdf"), format="pdf")
    plt.show()

    # plotting persistent diagram
    fig, ax = plt.subplots()

    ax = plot_dgms(ax, results[0].diagrams)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
