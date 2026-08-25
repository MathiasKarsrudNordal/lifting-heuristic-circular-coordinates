"""Plot the proportion of non-liftable lines over finite fields."""

from itertools import product
from math import isqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from numpy.random import Generator

sns.set_theme("talk", "ticks", palette="dark")

outdir = Path("figs")

DIMENSION = 6
K = 3
MAX_PRIME = 300
SAMPLE_SIZE = 1_000_000
SEED = 0


def primes_between(first: int, last: int) -> list[int]:
    """Return all primes in the closed interval [first, last]."""

    return [
        number
        for number in range(max(2, first), last + 1)
        if all(number % divisor for divisor in range(2, isqrt(number) + 1))
    ]


def allowed_residues(prime: int, k: int) -> list[int]:
    """Residues whose centered lift has size at most floor((p - 1) / k)."""

    bound = (prime - 1) // k
    return [*range(bound + 1), *range(prime - bound, prime)]


def normalize_line(vector: tuple[int, ...], prime: int) -> tuple[int, ...]:
    """Choose the representative whose first nonzero coordinate is one."""

    pivot = next(coordinate for coordinate in vector if coordinate != 0)
    inverse = pow(pivot, -1, prime)
    return tuple(inverse * coordinate % prime for coordinate in vector)


def exact_non_liftable_proportion(dimension: int, k: int, prime: int) -> float:
    """Count liftable projective lines without enumerating all projective lines."""

    allowed = allowed_residues(prime, k)
    liftable_lines = {
        normalize_line(vector, prime)
        for vector in product(allowed, repeat=dimension)
        if any(vector)
    }
    total_lines = (prime**dimension - 1) // (prime - 1)
    return 1.0 - len(liftable_lines) / total_lines


def sample_non_liftable_proportion(
    dimension: int,
    k: int,
    prime: int,
    sample_size: int,
    rng: Generator,
) -> float:
    """Estimate the proportion using uniformly sampled projective lines."""

    # Uniform nonzero vectors induce uniform projective lines because every line
    # has exactly p - 1 nonzero representatives.
    lines = rng.integers(0, prime, size=(sample_size, dimension), dtype=np.int64)
    zero_rows = np.all(lines == 0, axis=1)
    while np.any(zero_rows):
        lines[zero_rows] = rng.integers(
            0,
            prime,
            size=(int(zero_rows.sum()), dimension),
        )
        zero_rows = np.all(lines == 0, axis=1)

    # Normalize each line. Its first nonzero coordinate is then one, so a
    # liftable representative can only be obtained with an allowed scalar.
    pivot_columns = np.argmax(lines != 0, axis=1)
    pivots = lines[np.arange(sample_size), pivot_columns]
    inverses = np.array([0, *(pow(x, -1, prime) for x in range(1, prime))])
    unresolved = lines * inverses[pivots, None] % prime

    bound = (prime - 1) // k
    for scalar in allowed_residues(prime, k)[1:]:
        scaled = scalar * unresolved % prime
        centered_size = np.minimum(scaled, prime - scaled)
        is_liftable = np.all(centered_size <= bound, axis=1)
        unresolved = unresolved[~is_liftable]
        if len(unresolved) == 0:
            break

    return len(unresolved) / sample_size


def non_liftable_proportion(
    dimension: int,
    k: int,
    prime: int,
    sample_size: int,
    rng: Generator,
) -> float:
    """Use exact counting when feasible and sampling otherwise."""

    total_lines = (prime**dimension - 1) // (prime - 1)
    if total_lines < sample_size:
        return exact_non_liftable_proportion(dimension, k, prime)
    return sample_non_liftable_proportion(
        dimension,
        k,
        prime,
        sample_size,
        rng,
    )


def main() -> None:
    primes = primes_between(3, MAX_PRIME)
    rng = np.random.default_rng(SEED)
    proportions = []

    for prime in primes:
        proportion = non_liftable_proportion(
            DIMENSION,
            K,
            prime,
            SAMPLE_SIZE,
            rng,
        )
        proportions.append(proportion)
        print(f"p = {prime:3d}: {proportion:.6f}")

    x = [2, *primes]
    y = [0.0, *proportions]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(
        x,
        y,
        step="pre",
        alpha=0.8,
        edgecolor="darkslategray",
        linewidth=1.5,
    )
    ax.set(
        xlabel="Prime $p$",
        xticks=[3, 23, 47, 79, 109, 149, 191],
        xlim=(0, 200),
        ylim=(0, 1.05),
    )
    ax.grid(axis="y", alpha=0.3)
    sns.despine()
    fig.tight_layout()
    plt.savefig(outdir / Path("non-liftable-lines.pdf"), format="pdf")
    plt.show()


if __name__ == "__main__":
    main()
