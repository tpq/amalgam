## amalgam 0.0.4
---------------------
* New scalar product methods
    * New `adot` computes Aitchison scalar product for two vectors
    * New `amat` computes Aitchison scalar product for two matrices
    * New unit tests

## amalgam 0.0.3
---------------------
* Update `plot.amalgam`
    * Add argument to select which amalgams to show in ternary plot

## amalgam 0.0.2
---------------------
* Move `zCompositions` to Suggests

## amalgam 0.0.1
---------------------
* New `int` methods
    * `int.fromNbit` converts a binary vector into integers
* New `weight` methods
    * `weight.Nto1` builds a many-to-one amalgamation matrix
    * `weight.NtoN` builds a many-to-many amalgamation matrix
* New `objective` methods
    * `objective.keepDist` seeks to preserve Aitchison distance
    * `objective.maxRDA` seeks to maximize constraint of RDA
* New `amalgam` wrapper
    * This function runs a genetic algorithm for the objective
    * `print.amalgam` shows the object structure
    * `plot.amalgam` makes 4-panel figure
