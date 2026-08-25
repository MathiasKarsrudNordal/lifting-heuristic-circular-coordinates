"""Algorithm 1 specialized to m=1 from Section 4.1"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from icecream.icecream import ic
from numpy.typing import ArrayLike, NDArray
from scipy.sparse.linalg import lsmr

from lifting_cocycles.circular_coordinates import compute_circular_coordinates
from lifting_cocycles.complexes import Rips1Skeleton, cocycle_violations
from lifting_cocycles.datasets import torus_knot
from lifting_cocycles.utils import plot_dgms

sns.set_theme("talk", "ticks", palette="tab10")

outdir = Path("figs")

SEED = 0
PRIME = 47
N = 2500
NOISE_STD = np.sqrt(0.1)


@dataclass(frozen=True, slots=True)
class SpanningForest:
    """A rooted spanning forest in the skeleton's edge ordering."""

    parent: NDArray[np.int64]
    parent_edge: NDArray[np.int64]
    depth: NDArray[np.int64]
    traversal_order: NDArray[np.int64]
    tree_edges: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class PrimitiveReductionResult:
    """Exact data produced by the degree-one version of Algorithm 1."""

    primitive_cocycle: NDArray[np.int64]
    winding_number: int
    cycle: NDArray[np.int64]
    original_pairing: int
    final_pairing: int
    coboundary_potential: NDArray[np.int64]
    division_primes: tuple[int, ...]


def _spanning_forest(skeleton: Rips1Skeleton) -> SpanningForest:
    """Construct and root a spanning forest using union-find."""

    n_vertices = skeleton.n_vertices
    disjoint_set = np.arange(n_vertices, dtype=np.int64)
    ranks = np.zeros(n_vertices, dtype=np.int64)
    tree_edges = np.zeros(skeleton.n_edges, dtype=bool)
    tree_neighbors: list[list[tuple[int, int]]] = [[] for _ in range(n_vertices)]

    def find(vertex: int) -> int:
        root = vertex
        while int(disjoint_set[root]) != root:
            root = int(disjoint_set[root])
        while vertex != root:
            next_vertex = int(disjoint_set[vertex])
            disjoint_set[vertex] = root
            vertex = next_vertex
        return root

    for edge_index, (start, end) in enumerate(skeleton.edges):
        start, end = int(start), int(end)
        start_root = find(start)
        end_root = find(end)
        if start_root == end_root:
            continue

        if ranks[start_root] < ranks[end_root]:
            start_root, end_root = end_root, start_root
        disjoint_set[end_root] = start_root
        if ranks[start_root] == ranks[end_root]:
            ranks[start_root] += 1

        tree_edges[edge_index] = True
        tree_neighbors[start].append((end, edge_index))
        tree_neighbors[end].append((start, edge_index))

    parent = np.full(n_vertices, -1, dtype=np.int64)
    parent_edge = np.full(n_vertices, -1, dtype=np.int64)
    depth = np.zeros(n_vertices, dtype=np.int64)
    traversal_order: list[int] = []

    for root in range(n_vertices):
        if parent[root] != -1:
            continue

        parent[root] = root
        stack = [root]
        while stack:
            vertex = stack.pop()
            traversal_order.append(vertex)
            for neighbor, edge_index in tree_neighbors[vertex]:
                if parent[neighbor] != -1:
                    continue
                parent[neighbor] = vertex
                parent_edge[neighbor] = edge_index
                depth[neighbor] = depth[vertex] + 1
                stack.append(neighbor)

    return SpanningForest(
        parent=parent,
        parent_edge=parent_edge,
        depth=depth,
        traversal_order=np.asarray(traversal_order, dtype=np.int64),
        tree_edges=tree_edges,
    )


def _coboundary(
    potential: ArrayLike,
    skeleton: Rips1Skeleton,
) -> NDArray[np.int64]:
    """Evaluate the integral degree-zero coboundary edge by edge."""

    values = np.asarray(potential, dtype=np.int64)
    if values.shape != (skeleton.n_vertices,):
        raise ValueError("potential must contain one value per vertex")

    starts = skeleton.edges[:, 0]
    ends = skeleton.edges[:, 1]
    return values[ends] - values[starts]


def _boundary(
    chain: ArrayLike,
    skeleton: Rips1Skeleton,
) -> NDArray[np.int64]:
    """Evaluate the integral boundary of an edge chain."""

    coefficients = np.asarray(chain, dtype=np.int64)
    if coefficients.shape != (skeleton.n_edges,):
        raise ValueError("chain must contain one value per edge")

    boundary = np.zeros(skeleton.n_vertices, dtype=np.int64)
    np.add.at(boundary, skeleton.edges[:, 0], -coefficients)
    np.add.at(boundary, skeleton.edges[:, 1], coefficients)
    return boundary


def _tree_potential(
    cochain: ArrayLike,
    skeleton: Rips1Skeleton,
    forest: SpanningForest,
    *,
    modulus: int | None = None,
) -> NDArray[np.int64]:
    """Integrate a 1-cochain along the rooted spanning forest."""

    coefficients = np.asarray(cochain, dtype=np.int64)
    if coefficients.shape != (skeleton.n_edges,):
        raise ValueError("cochain must contain one value per edge")

    potential = np.zeros(skeleton.n_vertices, dtype=np.int64)
    for vertex_value in forest.traversal_order:
        vertex = int(vertex_value)
        parent = int(forest.parent[vertex])
        if parent == vertex:
            continue

        edge_index = int(forest.parent_edge[vertex])
        start, _ = map(int, skeleton.edges[edge_index])
        increment = int(coefficients[edge_index])
        if parent != start:
            increment = -increment

        potential[vertex] = potential[parent] + increment
        if modulus is not None:
            potential[vertex] %= modulus

    return potential


def _pairing(cochain: ArrayLike, chain: ArrayLike) -> int:
    """Compute the Kronecker pairing without NumPy integer overflow."""

    cochain_values = np.asarray(cochain)
    chain_values = np.asarray(chain)
    if cochain_values.shape != chain_values.shape:
        raise ValueError("cochain and chain must have the same shape")
    return sum(
        int(coefficient) * int(multiplicity)
        for coefficient, multiplicity in zip(
            cochain_values,
            chain_values,
            strict=True,
        )
    )


def _fundamental_cycle(
    cochain: ArrayLike,
    skeleton: Rips1Skeleton,
    forest: SpanningForest,
) -> tuple[NDArray[np.int64], int]:
    """Find a fundamental cycle with nonzero pairing with ``cochain``."""

    coefficients = np.asarray(cochain, dtype=np.int64)
    tree_potential = _tree_potential(coefficients, skeleton, forest)
    defects = coefficients - _coboundary(tree_potential, skeleton)
    candidates = np.flatnonzero((~forest.tree_edges) & (defects != 0))
    if candidates.size == 0:
        raise ValueError("the cocycle has zero pairing with every graph cycle")

    edge_index = min(candidates, key=lambda index: abs(int(defects[index])))
    start, end = map(int, skeleton.edges[edge_index])
    cycle = np.zeros(skeleton.n_edges, dtype=np.int64)
    cycle[edge_index] = 1

    # Complete the oriented non-tree edge start -> end by following the
    # unique tree path end -> start.
    left = end
    right = start
    while forest.depth[left] > forest.depth[right]:
        parent = int(forest.parent[left])
        parent_edge = int(forest.parent_edge[left])
        cycle[parent_edge] += 1 if left < parent else -1
        left = parent
    while forest.depth[right] > forest.depth[left]:
        parent = int(forest.parent[right])
        parent_edge = int(forest.parent_edge[right])
        cycle[parent_edge] -= 1 if right < parent else -1
        right = parent
    while left != right:
        left_parent = int(forest.parent[left])
        left_edge = int(forest.parent_edge[left])
        cycle[left_edge] += 1 if left < left_parent else -1
        left = left_parent

        right_parent = int(forest.parent[right])
        right_edge = int(forest.parent_edge[right])
        cycle[right_edge] -= 1 if right < right_parent else -1
        right = right_parent

    if np.any(_boundary(cycle, skeleton)):
        raise AssertionError("constructed fundamental chain is not a cycle")

    pairing = _pairing(coefficients, cycle)
    if pairing != int(defects[edge_index]):
        raise AssertionError("fundamental-cycle pairing has the wrong sign")
    if pairing < 0:
        cycle = -cycle
        pairing = -pairing
    return cycle, pairing


def _modular_coboundary_potential(
    cochain: ArrayLike,
    modulus: int,
    skeleton: Rips1Skeleton,
    forest: SpanningForest,
) -> NDArray[np.int64] | None:
    """Solve ``delta f = cochain`` modulo ``modulus``, if possible."""

    if modulus < 2:
        raise ValueError("modulus must be at least two")
    coefficients = np.asarray(cochain, dtype=np.int64)
    potential = _tree_potential(
        coefficients,
        skeleton,
        forest,
        modulus=modulus,
    )
    defects = _coboundary(potential, skeleton) - coefficients
    if np.any(defects % modulus):
        return None
    return potential


def _centered_lift(
    coefficients: ArrayLike,
    modulus: int,
) -> NDArray[np.int64]:
    """Choose the integer representatives closest to zero."""

    residues = np.asarray(coefficients, dtype=np.int64) % modulus
    midpoint = modulus // 2
    return ((residues + midpoint) % modulus - midpoint).astype(
        np.int64,
        copy=False,
    )


def _prime_divisors(value: int) -> tuple[int, ...]:
    """Return the distinct prime divisors of a nonzero integer."""

    remaining = abs(int(value))
    if remaining == 0:
        raise ValueError("zero has no finite set of prime divisors")

    divisors: list[int] = []
    if remaining % 2 == 0:
        divisors.append(2)
        while remaining % 2 == 0:
            remaining //= 2

    candidate = 3
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            divisors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate += 2
    if remaining > 1:
        divisors.append(remaining)
    return tuple(divisors)


def reduce_to_primitive_cocycle(
    cocycle: ArrayLike,
    skeleton: Rips1Skeleton,
    *,
    cycle: ArrayLike | None = None,
) -> PrimitiveReductionResult:
    """Reduce an integer 1-cocycle using Algorithm 1 of the manuscript."""

    original = np.asarray(cocycle, dtype=np.int64)
    if original.shape != (skeleton.n_edges,):
        raise ValueError("cocycle must contain one value per oriented edge")
    if cocycle_violations(skeleton, original).size:
        raise ValueError("input is not an integer 1-cocycle")

    forest = _spanning_forest(skeleton)
    if cycle is None:
        beta, original_pairing = _fundamental_cycle(
            original,
            skeleton,
            forest,
        )
    else:
        beta = np.asarray(cycle, dtype=np.int64)
        if beta.shape != (skeleton.n_edges,):
            raise ValueError("cycle must contain one value per oriented edge")
        if np.any(_boundary(beta, skeleton)):
            raise ValueError("the supplied edge chain is not an integer cycle")
        original_pairing = _pairing(original, beta)
        if original_pairing == 0:
            raise ValueError("the cocycle-cycle pairing must be nonzero")
        if original_pairing < 0:
            beta = -beta
            original_pairing = -original_pairing

    current = original.copy()
    winding_number = 1
    correction = np.zeros(skeleton.n_vertices, dtype=np.int64)
    division_primes: list[int] = []
    candidate_primes = _prime_divisors(original_pairing)

    for prime in candidate_primes:
        while True:
            modular_potential = _modular_coboundary_potential(
                current,
                prime,
                skeleton,
                forest,
            )
            if modular_potential is None:
                break

            potential = _centered_lift(modular_potential, prime)
            numerator = current - _coboundary(potential, skeleton)
            if np.any(numerator % prime):
                raise AssertionError("division numerator is not divisible by prime")
            quotient = numerator // prime

            if cocycle_violations(skeleton, quotient).size:
                raise AssertionError("division did not produce an integer cocycle")
            if _pairing(current, beta) != prime * _pairing(quotient, beta):
                raise AssertionError("Kronecker pairing did not divide by prime")

            correction += winding_number * potential
            winding_number *= prime
            current = quotient
            division_primes.append(prime)

            reconstructed = winding_number * current + _coboundary(correction, skeleton)
            if not np.array_equal(original, reconstructed):
                raise AssertionError("global cochain decomposition was not preserved")

    for prime in candidate_primes:
        if (
            _modular_coboundary_potential(
                current,
                prime,
                skeleton,
                forest,
            )
            is not None
        ):
            raise AssertionError("a candidate prime remains removable")

    return PrimitiveReductionResult(
        primitive_cocycle=current,
        winding_number=winding_number,
        cycle=beta,
        original_pairing=original_pairing,
        final_pairing=_pairing(current, beta),
        coboundary_potential=correction,
        division_primes=tuple(division_primes),
    )


def _circular_angles(
    cocycle: ArrayLike,
    skeleton: Rips1Skeleton,
) -> NDArray[np.float64]:
    """Apply the existing least-squares smoothing step to one cocycle."""

    least_squares = lsmr(
        skeleton.delta0,
        np.asarray(cocycle, dtype=float),
        atol=1e-10,
        btol=1e-10,
    )
    if int(least_squares[1]) == 7:
        raise RuntimeError("LSMR reached its iteration limit before converging")
    return np.mod(least_squares[0], 1.0) * (2.0 * np.pi)


def main():
    X = torus_knot(
        N,
        seed=SEED,
        p=2,
        q=3,
        R=5,
        r=2,
        noise_std=NOISE_STD,
    )

    result = compute_circular_coordinates(
        X[:, :3],
        prime=PRIME,
        coordinate_fraction=0.75,
        lift_scalar=None,
    )
    reduction = reduce_to_primitive_cocycle(
        result.integer_cocycle,
        result.skeleton,
    )
    primitive_angles = _circular_angles(
        reduction.primitive_cocycle,
        result.skeleton,
    )
    ic(
        result.scalar,
        reduction.original_pairing,
        reduction.final_pairing,
        reduction.winding_number,
        reduction.division_primes,
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].scatter(
        *X[:, :2].T,
        c=result.angles,
        cmap="winter",
        s=4,
        alpha=0.8,
    )
    axes[0].set(
        # title=f"Winding number {reduction.winding_number}",
        xticks=[],
        yticks=[],
        aspect="equal",
    )
    axes[1].scatter(
        *X[:, :2].T,
        c=primitive_angles,
        cmap="winter",
        s=4,
        alpha=0.8,
    )
    axes[1].set(
        # title="Primitive cocycle",
        xticks=[],
        yticks=[],
        aspect="equal",
    )

    plt.tight_layout()
    plt.savefig(outdir / Path("winding-reduction.pdf"), format="pdf")
    plt.show()

    # plotting persistent diagram
    fig, ax = plt.subplots()

    ax = plot_dgms(ax, result.diagrams)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
