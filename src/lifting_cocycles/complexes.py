"""Minimal simplicial-complex utilities for circular coordinates.

Ripser constructs the Vietoris--Rips filtration internally.  Once a scale has
been selected, circular-coordinate smoothing only needs its one-skeleton and
the degree-zero coboundary operator

    delta_0: C^0(X; R) -> C^1(X; R).

Edges are oriented from their smaller to their larger vertex index throughout
this module.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import csr_array


@dataclass(frozen=True, slots=True)
class Rips1Skeleton:
    """The oriented edges and degree-zero coboundary at one Rips scale."""

    n_vertices: int
    threshold: float
    edges: NDArray[np.int64]
    delta0: csr_array

    @property
    def n_edges(self) -> int:
        return self.edges.shape[0]


def cocycle_violations(
    skeleton: Rips1Skeleton,
    cochain: ArrayLike,
) -> NDArray[np.int64]:
    """Return nonzero integer coboundaries on triangles in the skeleton."""
    values = np.asarray(cochain, dtype=np.int64)
    if values.shape != (skeleton.n_edges,):
        raise ValueError("cochain must contain one value per oriented edge")

    forward_neighbors = [set() for _ in range(skeleton.n_vertices)]
    edge_values = {}
    for (start, end), value in zip(skeleton.edges, values, strict=True):
        start, end = int(start), int(end)
        forward_neighbors[start].add(end)
        edge_values[start, end] = int(value)

    violations = []
    for first in range(skeleton.n_vertices):
        for second in forward_neighbors[first]:
            common_neighbors = (
                forward_neighbors[first] & forward_neighbors[second]
            )
            for third in common_neighbors:
                coboundary = (
                    edge_values[second, third]
                    - edge_values[first, third]
                    + edge_values[first, second]
                )
                if coboundary != 0:
                    violations.append(coboundary)

    return np.asarray(violations, dtype=np.int64)


def rips_1_skeleton(
    distance_matrix: ArrayLike,
    threshold: float,
) -> Rips1Skeleton:
    """Construct the one-skeleton of a Vietoris--Rips complex.

    An edge ``(u, v)`` with ``u < v`` is present when its distance is at most
    ``threshold``.  The corresponding row of ``delta0`` evaluates a vertex
    cochain ``f`` as ``f(v) - f(u)``.
    """

    distances = np.asarray(distance_matrix, dtype=float)
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("distance_matrix must be a square matrix")
    if distances.shape[0] == 0:
        raise ValueError("distance_matrix must contain at least one vertex")
    if np.isnan(distances).any():
        raise ValueError("distance_matrix must not contain NaN values")
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    if not np.allclose(distances, distances.T):
        raise ValueError("distance_matrix must be symmetric")

    n_vertices = distances.shape[0]
    starts, ends = np.triu_indices(n_vertices, k=1)
    present = distances[starts, ends] <= threshold
    edges = np.column_stack((starts[present], ends[present])).astype(
        np.int64,
        copy=False,
    )

    n_edges = edges.shape[0]
    rows = np.repeat(np.arange(n_edges, dtype=np.int64), 2)
    columns = edges.reshape(-1)
    values = np.tile(np.array([-1.0, 1.0]), n_edges)
    delta0 = csr_array(
        (values, (rows, columns)),
        shape=(n_edges, n_vertices),
    )

    return Rips1Skeleton(
        n_vertices=n_vertices,
        threshold=float(threshold),
        edges=edges,
        delta0=delta0,
    )


def sparse_cocycle_to_vector(
    cocycle: ArrayLike,
    edges: NDArray[np.int64],
    modulus: int,
    *,
    scalar: int = 1,
) -> NDArray[np.int64]:
    """Put a sparse Ripser 1-cocycle into the canonical edge ordering.

    Ripser stores a nonzero term as ``(u, v, coefficient)`` and may orient the
    edge in either direction.  Terms supported on edges outside the selected
    Rips scale are discarded, matching restriction of the cocycle to that
    complex.
    """

    if modulus <= 1:
        raise ValueError("modulus must be greater than one")
    scalar %= modulus
    if scalar == 0:
        raise ValueError("scalar must be nonzero modulo the modulus")

    sparse_cocycle = np.asarray(cocycle)
    if sparse_cocycle.ndim != 2 or sparse_cocycle.shape[1] != 3:
        raise ValueError("cocycle must have shape (n_terms, 3)")

    edge_to_index = {
        (int(start), int(end)): index
        for index, (start, end) in enumerate(edges)
    }
    vector = np.zeros(edges.shape[0], dtype=np.int64)

    for first, second, coefficient in sparse_cocycle:
        first = int(first)
        second = int(second)
        if first == second:
            raise ValueError("a 1-cocycle term cannot be supported on a vertex")

        if first < second:
            edge = (first, second)
            orientation = 1
        else:
            edge = (second, first)
            orientation = -1

        index = edge_to_index.get(edge)
        if index is not None:
            vector[index] = (
                vector[index]
                + orientation * scalar * int(coefficient)
            ) % modulus

    return vector
