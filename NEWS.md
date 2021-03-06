## amalgam 0.1.8
---------------------
* New zero-handling procedure
    * If `shrink = TRUE`, bin frequencies are estimated by James-Stein-type shrinkage
* New objective functions
    * Add `asSLR` support to `objective.keepWADIST` by comparing with Euclidean distance
    * `objective.maxDSKL` maximizes the RDA of a PCoA of the SKL divergences

## amalgam 0.1.7
---------------------
* New objective function
    * `objective.keepSKL` preserves the symmetric Kullback-Leibler (SKL) divergences

## amalgam 0.1.6
---------------------
* New objective function
    * `objective.keepWADIST` preserves the weighted Aitchison distances

## amalgam 0.1.5
---------------------
* Update information objectives
    * `objective.keepEntropy` now recloses the amalgamation (as it should)
    * `objective.maxDEntropy` now recloses the amalgamation (as it should)

## amalgam 0.1.4
---------------------
* New backend information functions
    * `shannon` function calculates Shannon's index
    * `uniformity` function calculates a normalized Shannon's index
* New objective functions
    * `objective.keepEntropy` preserves the relative entropy of the amalgams
    * `objective.maxDEntropy` maximizes the between-group difference in amalgam entropy
* Update `amalgam` function
    * `prepareArgs` now forces a re-closure of the data

## amalgam 0.1.3
---------------------
* Switch `objective.maxRDA` and `objective.maxRDA2` (again)

## amalgam 0.1.2
---------------------
* Switch `objective.maxRDA` and `objective.maxRDA2`

## amalgam 0.1.1
---------------------
* Update `objective.maxRDA` and `objective.maxRDA2`
    * Now sums inertia across all RD axes (as it should)

## amalgam 0.1.0
---------------------
* Update `objective.maxRDA2` method
    * Like `objective.maxRDA` but makes "z" the community matrix

## amalgam 0.0.9
---------------------
* Update `objective.maxRDA` method
    * Handle rare case when constrained axes explain no variance
    * Add check for "z" argument

## amalgam 0.0.8
---------------------
* Have `prepareArgs` coerce "z" as a `data.frame`
* Update `objective.maxRDA` method
    * Now save the intermediate data to debug

## amalgam 0.0.7
---------------------
* Update `objective.maxRDA` method
    * Add `tryCatch` to handle the Lapack routine 'dgesdd' error
    * This error occurs when the SVD fails to converge

## amalgam 0.0.6
---------------------
* New `as.slr` function turns components into a set of summed log-ratios
* Update `amalgam`
    * New `asSLR` argument passed to objective functions
    * Objective functions now support SLRs
    * Update `plot` method for SLRs
    
## amalgam 0.0.5
---------------------
* Update `amalgam`
    * Zeros now replaced with `zCompositions::cmultRepl`
    * Correctly coerce input as matrix

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
