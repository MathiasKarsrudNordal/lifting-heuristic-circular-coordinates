# Lifting cocycles and circular coordinates

This project contains a minimal Python implementation of circular coordinates
based on representative persistent cocycles returned by Ripser.

```python
from lifting_cocycles.circular_coordinates import compute_circular_coordinates
from lifting_cocycles.datasets import circle

data = circle(500, seed=0, noise_std=0.02)
result = compute_circular_coordinates(data[:, :2], prime=47)

angles = result.angles  # values in [0, 2*pi)
```

The implementation constructs the selected Vietoris--Rips one-skeleton and
its sparse degree-zero coboundary operator. It directly applies the centered
integer lift to Ripser's finite-field cocycle, then smooths it with SciPy's
sparse LSMR solver. It currently neither detects nor corrects a failed lift.

Run the tests with:

```bash
uv run python -m unittest discover -s tests
```
