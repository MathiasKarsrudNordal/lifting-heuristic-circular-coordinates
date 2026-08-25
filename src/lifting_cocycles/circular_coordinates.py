"""Compute circular coordinates from a persistent 1-cocycle."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from ripser import ripser
from scipy.sparse.linalg import lsmr
from scipy.spatial.distance import pdist, squareform

from .complexes import (
    Rips1Skeleton,
    cocycle_violations,
    rips_1_skeleton,
    sparse_cocycle_to_vector,
)


@dataclass(frozen=True, slots=True)
class CircularCoordinatesResult:
    """The coordinates and intermediate data used to construct them."""

    angles: NDArray[np.float64]
    integer_cocycle: NDArray[np.int64]
    finite_field_cocycle: NDArray[np.int64]
    skeleton: Rips1Skeleton
    diagrams: tuple[NDArray[np.float64], ...]
    cocycle_index: int
    birth: float
    death: float
    scalar: int
    prime: int
    lsmr_stop_code: int
    lsmr_iterations: int
    lsmr_residual_norm: float


@dataclass(frozen=True, slots=True)
class PersistentCocycle:
    """A selected Ripser cocycle and the data needed to smooth it."""

    distance_matrix: NDArray[np.float64]
    sparse_cocycle: NDArray[np.int64]
    diagrams: tuple[NDArray[np.float64], ...]
    cocycle_index: int
    birth: float
    death: float
    prime: int
    persistence_threshold: float


def centered_lift(
    coefficients: ArrayLike,
    prime: int,
) -> NDArray[np.int64]:
    """Lift to the integer closest to zero."""

    residues = np.asarray(coefficients, dtype=np.int64) % prime
    midpoint = prime // 2
    return ((residues + midpoint) % prime - midpoint).astype(
        np.int64,
        copy=False,
    )


def compute_circular_coordinates(
    points: ArrayLike,
    *,
    prime: int = 47,
    persistence_threshold: float = np.inf,
    coordinate_threshold: float | None = None,
    coordinate_fraction: float = 0.25,
    cocycle_index: int | None = None,
    lift_scalar: int | None = 1,
    return_scaled_coordinate: bool = True,
) -> CircularCoordinatesResult:
    """Compute one circular coordinate using persistent cohomology.

    Parameters
    ----------
    points:
        Array of shape ``(n_points, ambient_dimension)``.
    prime:
        Odd prime used by Ripser.
    persistence_threshold:
        Maximum edge length passed to Ripser.  A finite value is useful for
        large point clouds.
    coordinate_threshold:
        Rips scale on which to construct the coordinate.  If omitted, choose
        a scale inside the selected class's persistence interval.
    coordinate_fraction:
        Position in the persistence interval used when
        ``coordinate_threshold`` is omitted.  The default 0.25 avoids the
        dense, nearly filled-in complex close to the death scale.
    cocycle_index:
        Index in Ripser's one-dimensional persistence diagram.  By default,
        select the class with greatest persistence.
    lift_scalar:
        Nonzero finite-field scalar applied before centered lifting.  The
        default uses the representative returned by Ripser without scaling.

    Notes
    -----
    The centered lift is used directly.  This implementation does not check
    whether it is an integer 1-cocycle or correct it when it is not.
    """

    persistent_cocycle = compute_persistent_cocycle(
        points,
        prime=prime,
        persistence_threshold=persistence_threshold,
        cocycle_index=cocycle_index,
    )
    return compute_circular_coordinates_from_cocycle(
        persistent_cocycle,
        coordinate_threshold=coordinate_threshold,
        coordinate_fraction=coordinate_fraction,
        lift_scalar=lift_scalar,
        return_scaled_coordinate=return_scaled_coordinate,
    )


def compute_persistent_cocycle(
    points: ArrayLike,
    *,
    prime: int = 47,
    persistence_threshold: float = np.inf,
    cocycle_index: int | None = None,
) -> PersistentCocycle:
    """Run Ripser and select one persistent one-dimensional cocycle."""

    point_cloud = np.asarray(points, dtype=float)
    if point_cloud.ndim != 2:
        raise ValueError("points must have shape (n_points, ambient_dimension)")
    if point_cloud.shape[0] < 3:
        raise ValueError("at least three points are required")
    if not np.isfinite(point_cloud).all():
        raise ValueError("points must contain only finite values")
    if not _is_odd_prime(prime):
        raise ValueError("prime must be an odd prime")
    if persistence_threshold <= 0:
        raise ValueError("persistence_threshold must be positive")

    distance_matrix = squareform(pdist(point_cloud))
    ripser_result = ripser(
        distance_matrix,
        distance_matrix=True,
        maxdim=1,
        coeff=prime,
        do_cocycles=True,
        thresh=persistence_threshold,
    )

    diagrams = tuple(
        np.asarray(persistence_diagram, dtype=float)
        for persistence_diagram in ripser_result["dgms"]
    )
    diagram = diagrams[1]
    if diagram.size == 0:
        raise ValueError("Ripser found no one-dimensional persistent class")

    selected_index = _select_cocycle_index(
        diagram,
        cocycle_index,
        persistence_threshold,
    )
    birth, death = map(float, diagram[selected_index])

    return PersistentCocycle(
        distance_matrix=distance_matrix,
        sparse_cocycle=np.asarray(
            ripser_result["cocycles"][1][selected_index],
            dtype=np.int64,
        ),
        diagrams=diagrams,
        cocycle_index=selected_index,
        birth=birth,
        death=death,
        prime=prime,
        persistence_threshold=float(persistence_threshold),
    )


def compute_circular_coordinates_from_cocycle(
    persistent_cocycle: PersistentCocycle,
    *,
    coordinate_threshold: float | None = None,
    coordinate_fraction: float = 0.25,
    lift_scalar: int | None = 1,
    return_scaled_coordinate: bool = True,
) -> CircularCoordinatesResult:
    """Construct circular coordinates from a previously computed cocycle."""

    birth = persistent_cocycle.birth
    death = persistent_cocycle.death
    prime = persistent_cocycle.prime
    scale = _coordinate_scale(
        birth,
        death,
        coordinate_threshold,
        coordinate_fraction,
        persistent_cocycle.persistence_threshold,
    )

    skeleton = rips_1_skeleton(persistent_cocycle.distance_matrix, scale)
    sparse_cocycle = persistent_cocycle.sparse_cocycle

    # no provided scaler; search over all scaleres until succesfull lift is guaranteed
    if lift_scalar is None:
        for scalar in range(1, prime):
            finite_field_cocycle = sparse_cocycle_to_vector(
                sparse_cocycle,
                skeleton.edges,
                prime,
                scalar=scalar,
            )
            scalar_integer_cocycle = centered_lift(finite_field_cocycle, prime)
            integer_cocycle = scalar_integer_cocycle

            # scale in Z to obtain lirft of original cocycle
            if not return_scaled_coordinate:
                inverse_scalar = pow(scalar, -1, prime)
                if inverse_scalar > prime // 2:
                    inverse_scalar -= prime
                integer_cocycle = inverse_scalar * scalar_integer_cocycle

            if cocycle_violations(skeleton, integer_cocycle).size == 0:
                break
        else:
            raise ValueError("no scalar gives an integer-cocycle centered lift")

    # provided scaler
    else:
        scalar = int(lift_scalar) % prime

        if scalar != 0:
            finite_field_cocycle = sparse_cocycle_to_vector(
                sparse_cocycle,
                skeleton.edges,
                prime,
                scalar=scalar,
            )
            if not np.any(finite_field_cocycle):
                raise ValueError(
                    "the selected Ripser cocycle has empty support at the coordinate threshold"
                )
            integer_cocycle = centered_lift(finite_field_cocycle, prime)

            if not return_scaled_coordinate:
                inverse_scalar = pow(scalar, -1, prime)
                if inverse_scalar > prime // 2:
                    inverse_scalar -= prime
                integer_cocycle = inverse_scalar * integer_cocycle

        elif scalar == 0:
            raise ValueError("lift_scalar must be nonzero modulo prime")

    least_squares = lsmr(
        skeleton.delta0,
        integer_cocycle.astype(float),
        atol=1e-10,
        btol=1e-10,
    )
    vertex_potential = least_squares[0]
    stop_code = int(least_squares[1])
    if stop_code == 7:
        raise RuntimeError("LSMR reached its iteration limit before converging")

    angles = np.mod(vertex_potential, 1.0) * (2.0 * np.pi)
    return CircularCoordinatesResult(
        angles=angles,
        integer_cocycle=integer_cocycle,
        finite_field_cocycle=finite_field_cocycle,
        skeleton=skeleton,
        diagrams=persistent_cocycle.diagrams,
        cocycle_index=persistent_cocycle.cocycle_index,
        birth=birth,
        death=death,
        scalar=scalar,
        prime=prime,
        lsmr_stop_code=stop_code,
        lsmr_iterations=int(least_squares[2]),
        lsmr_residual_norm=float(least_squares[3]),
    )


def _select_cocycle_index(
    diagram: NDArray[np.float64],
    requested_index: int | None,
    persistence_threshold: float,
) -> int:
    if requested_index is not None:
        if requested_index < 0 or requested_index >= diagram.shape[0]:
            raise IndexError("cocycle_index is outside the persistence diagram")
        return requested_index

    effective_deaths = diagram[:, 1].copy()
    infinite = np.isinf(effective_deaths)
    if infinite.any() and np.isfinite(persistence_threshold):
        effective_deaths[infinite] = persistence_threshold
    persistence = effective_deaths - diagram[:, 0]
    return int(np.argmax(persistence))


def _coordinate_scale(
    birth: float,
    death: float,
    requested_scale: float | None,
    coordinate_fraction: float,
    persistence_threshold: float,
) -> float:
    if requested_scale is None:
        if not 0.0 < coordinate_fraction < 1.0:
            raise ValueError("coordinate_fraction must lie strictly between 0 and 1")

        if np.isfinite(death):
            interval_end = death
        elif np.isfinite(persistence_threshold):
            interval_end = float(persistence_threshold)
        else:
            raise ValueError(
                "coordinate_threshold is required when the selected class "
                "has infinite death and persistence_threshold is infinite"
            )
        scale = birth + coordinate_fraction * (interval_end - birth)
    else:
        scale = float(requested_scale)

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("coordinate_threshold must be finite and positive")
    if scale <= birth:
        raise ValueError(
            f"coordinate_threshold={scale} must exceed class birth={birth}"
        )
    if np.isfinite(death) and scale >= death:
        raise ValueError(
            f"coordinate_threshold={scale} must be below class death={death}"
        )
    if np.isfinite(persistence_threshold) and scale > persistence_threshold:
        raise ValueError("coordinate_threshold cannot exceed persistence_threshold")
    return scale


def _is_odd_prime(value: int) -> bool:
    if not isinstance(value, (int, np.integer)) or value < 3 or value % 2 == 0:
        return False
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True
