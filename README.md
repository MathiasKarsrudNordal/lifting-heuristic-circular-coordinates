# Lifting Cocycles: From Heuristic to Theory

The circular coordinates algorithm extracts circle-valued features from data by
lifting persistent cocycles from a finite field to the integers. This repository
accompanies our study of when the commonly used coefficientwise lift succeeds,
how finite-field rescaling can make it succeed, and how an integral cocycle can
be reduced to winding number 1. We also extend the lifting principles to
homology cycles.

## Scaling before lifting

[![Circular coordinates on a trefoil before and after rescaling the finite-field cocycle](figs/trefoil-scaling.png)](figs/trefoil-scaling.pdf)

*For $p = 47$, directly lifting the cocycle returned by Ripser fails to produce
the expected circular coordinates (left). Scaling the cocycle by 2 in the
finite field before lifting yields an integral cocycle and recovers the circular
coordinate (right). Click the figure for the vector PDF.*
