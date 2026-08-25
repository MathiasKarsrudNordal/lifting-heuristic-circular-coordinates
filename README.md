# Lifting Cocycles: From Heuristic to Theory

The circular coordinates algorithm extracts circle-valued features from data by
lifting persistent cocycles from a finite field to the integers. This repository
accompanies our study of when the commonly used coefficientwise lift succeeds,
how finite-field rescaling can make it succeed, and how an integral cocycle can
be reduced to winding number 1. We also extend the lifting principles to
homology cycles.

## Scaling before lifting

[![Circular coordinates on a trefoil before and after rescaling the finite-field cocycle](figs/trefoil-scaling.png)](figs/trefoil-scaling.pdf)


#### Reproduction of figures

**Fig. 2:**
Run `scripts/circular_coordinates/trefoil.py`.

**Fig. 3:**
Run `scripts/non_liftable_lines.py`.

**Fig. 4:**
Run `scripts/circular_coordinates/multiple_windings.py`.

**Fig. 5:**
Run `scripts/winding_number_reduction.py`.

**Fig. 6:**
Run `scripts/circular_coordinates/high_dimensional_circle.py`.
