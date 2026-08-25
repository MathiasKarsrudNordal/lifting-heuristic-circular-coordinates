"""Datasets"""

import numpy as np


def circle(N, *, seed, noise_std=0.01):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, N)
    noise = noise_std * rng.normal(size=(N, 2))

    x = np.cos(t)
    y = np.sin(t)

    coords = np.stack([x, y], axis=-1) + noise
    return np.column_stack((coords, t))


def high_dim_circle(n_points, embedding_dim):
    if embedding_dim < 2:
        raise ValueError("embedding_dim must be at least 2")

    t = np.linspace(0, 2 * np.pi, num=n_points, endpoint=False)
    embedding = np.zeros((n_points, embedding_dim))

    for i, frequency in enumerate(range(1, embedding_dim // 2 + 1)):
        embedding[:, 2 * i] = np.cos(frequency * t)
        embedding[:, 2 * i + 1] = np.sin(frequency * t)

    if embedding_dim % 2:
        frequency = embedding_dim // 2 + 1
        embedding[:, -1] = np.cos(frequency * t)

    return t, embedding


def trefoil_knot(N, *, seed, noise_std=0.01):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, N)
    noise = noise_std * rng.normal(size=(N, 3))

    x = np.sin(t) + 2 * np.sin(2 * t)
    y = np.cos(t) - 2 * np.cos(2 * t)
    z = -np.sin(3 * t)

    coords = np.stack([x, y, z], axis=-1) + noise
    return np.column_stack((coords, t))


def torus_knot(
    N,
    *,
    seed,
    p=2,
    q=3,
    R=2.0,
    r=1.0,
    noise_std=0,
):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, N)

    radial_coordinate = R + r * np.cos(q * t)
    x = radial_coordinate * np.cos(p * t)
    y = radial_coordinate * np.sin(p * t)
    z = r * np.sin(q * t)

    coords = np.stack([x, y, z], axis=-1)
    coords += noise_std * rng.normal(size=(N, 3))
    return np.column_stack((coords, t))
