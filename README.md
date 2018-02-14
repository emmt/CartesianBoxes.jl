# Simple and efficient bounding-boxes in Julia

This module implements bounding-boxes which are aliases to Cartesian ranges
and methods to easily construct bounding-boxes and help writing fast
multi-dimensional algorithms.

A bounding-box is defined by `BoundingBox{N}` and `BoundingBox(args...)`
are shortcuts for `CartesianRange{CartesianIndex{N}}` and
`CartesianRange(args...)`.  In addition to the usual ways to build a
Cartesian range, the following methods are provided:

```julia
BoundingBox(A) -> BoundingBox(indices(A))
BoundingBox((imin,jmin,...), (imax,jmax,...))
BoundingBox(imin:imax, jmin:jmax, ...)
```

where `A` is an array and `imin`, `imax`, *etc.* are integers.  These
definitions do not interfere with those of `CartesianRange`.

The call:

```julia
intersection(R, S)
```

yields the bounding-box (*i.e.* the Cartesian range) which is the
intersection of the two Cartesian ranges `R` and `S`.  This method is
similar to `intersect(R,S) = R ∩ S` which yields an array of Cartesian
indices and is **much** slower (and hence less useful).

The call:

```julia
isnonemptypartof(R, S)
```

yields whether region `R` is nonempty and a valid part of `S`.  If `R` is a
Cartesian range, then `S` can be an array to check whether `R` is a valid
nonempty region of interest of `S`.

The call:

```julia
boundingbox(A [, B])
```

yields the bounding-box of non-zero values in array `A`.  If optional
argumen `B` is provided, then the bounding-box of the sub-region `B` of `A`
is returned.


## Restrictions

The algorithm for find the bounding-box of non-zero values is pretty simple
and scales as `O(length(A))`.


## Installation

`BoundingBoxes.jl` is not yet an
[official Julia package](https://pkg.julialang.org/) so you have to clone the
repository to install the package:

```julia
Pkg.clone("https://github.com/emmt/BoundingBoxes.jl.git")
```

There is nothing to build so no needs to call `Pkg.build("BoundingBoxes")`.

Later, it is sufficient to do:

```julia
Pkg.update("BoundingBoxes")
```

to pull the latest version.
