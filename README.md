
<!-- README.md is generated from README.Rmd. Please edit that file -->
Quick start
-----------

Welcome to the `amalgam` GitHub page!

It is difficult to visualize high-dimensional compositional data. One approach involves amalgamating the parts into a new set of composite features called "amalgams". This package uses genetic algorithms to find an optimal set of amalgams that satisfies an objective function. Although amalgams themselves are low-dimensional compositions, they are not sub-compositions. Rather, amalgamation acts like a kind of non-linear dimension reduction that can facilitate data understanding. Note that inferences made on amalgamations may not agree with inferences made on the full composition.

``` r
library(devtools)
devtools::install_github("tpq/amalgam")
library(amalgam)
?amalgam
```

Data-driven amalgamation
------------------------

Before using `amalgam`, the analyst should decide what kind of amalgamation they want.

First, the analyst should decide *how* to amalgamate. The simplest amalgamation is "N-to-1" (handled by the `weight.Nto1` function), where each component only contributes to one amalgam. A more complex amalgamation is "N-to-N" (handled by the `weight.NtoN` function), where each component may contribute to multiple amalgams. The former would produce amalgams that are easier to interpret, but the latter may fit the objective better.

Second, the analyst should decide on the objective. One objective function is to preserve the Aitchison distance between the sample amalgams (see `objective.keepDist`). This is an unsupervised dimension reduction method akin to a PCoA. Another objective function is to maximize the variance explained by a constraining matrix, `z` (see `objective.maxRDA`). This is a supervised dimension reduction method akin to a non-linear RDA.

Once this is decided, running the algorithm is easy! Use the `n.amalgams` argument to choose how many amalgams you want. Use the `maxiter` argument to choose how long the algorithm runs.

``` r
set.seed(3220)
sampleData <- randAcomp(10, 15)
out <- amalgam(sampleData, n.amalgams = 3, z = NULL,
               objective = objective.keepDist,
               weight = weight.Nto1)
print(out)
#> AMALGAM OBJECT
#> 
#> Dimension of original data:
#> [1] 10 15
#> 
#> Dimension of amalgamated data:
#> [1] 10  3
#> 
#> Size of amalgams:
#> [1] 3 4 3
#> 
#> Size of contributions:
#> 
#>  0  1 
#>  5 10 
#> 
#> Fitness:
#> [1] 0.8347408
#> 
#> Use plot() method!
```

Visualization
-------------

``` r
fakeClusters <- cutree(hclust(dist(sampleData)), 3)
plot(out, col = fakeClusters, center = TRUE, scale = TRUE)
```

![](README-unnamed-chunk-5-1.png)

This package only offers a few objectives right now, but we will add more later.
