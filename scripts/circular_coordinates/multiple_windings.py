"""Reproducing figure 4 in the paper"""

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from lifting_cocycles.circular_coordinates import compute_circular_coordinates
from lifting_cocycles.datasets import circle

sns.set_theme("talk", "ticks", palette="tab10")

outdir = Path("figs")

SEED = 0
PRIME = 47
N = 1500
LIFT_SCALARS = (1, 2, 5, 10)


def main():

    X = circle(N, seed=SEED, noise_std=0.05)

    results = [
        compute_circular_coordinates(
            X[:, :2],
            prime=PRIME,
            lift_scalar=scalar,
        )
        for scalar in LIFT_SCALARS
    ]

    fig, axes = plt.subplots(1, len(LIFT_SCALARS), figsize=(18, 5))

    for ax, result in zip(
        axes,
        results,
        strict=True,
    ):
        ax.scatter(
            *X[:, :2].T,
            c=result.angles,
            cmap="winter",
            s=2.5,
            alpha=0.8,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(outdir / Path("multiple-windings-circle.pdf"), format="pdf")
    plt.show()


if __name__ == "__main__":
    main()
