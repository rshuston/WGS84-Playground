# WGS84-Playground

A workspace to experiment with and compare WGS84 coordinate conversion algorithms.

Only one geodetic-to-ECEF conversion function is coded, since the conversion is
straightforward.

Four ECEF-to-geodetic conversion functions are coded:
- Iterative method - typically converges in less than 5 iterations for terrestrial locations
- Bowring's method - where a single iteration is sufficiently accurate for most applications
- Heikkenin's method - optimized to compute fewest terms, but requires cube root
- Olson's method - accurate as Heikkinen, less computations, doesn't needccube root, and 63% quicker

The main program performs timing and accuracy comparisons.

The unit tests exercise significant earth positions and some hand-picked positions.
