# Simple and efficient Cartesian boxes in Julia

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/CartesianBoxes.jl.svg?branch=master)](https://travis-ci.org/emmt/CartesianBoxes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/CartesianBoxes.jl?branch=master)](https://ci.appveyor.com/project/emmt/CartesianBoxes-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/emmt/CartesianBoxes.jl/badge.svg?branch=master)](https://coveralls.io/github/emmt/CartesianBoxes.jl?branch=master)
[![codecov.io](http://codecov.io/github/emmt/CartesianBoxes.jl/coverage.svg?branch=master)](http://codecov.io/github/emmt/CartesianBoxes.jl?branch=master)

This module implements rectangular regions of `N`-dimensional indices which are
instances of `CartesianBox{N}`.  For multi-dimensional loops, these instances
are as fast as `CartesianIndices{N}` or as `CartesianRange{CardinalIndex{N}}`
(depending on your version of Julia).  They can thus be used as a *fast* and
*portable* replacement (see [*Rationale*](#rationale) below) of these different
representations of rectangular multi-dimensional regions.

A bounding-box is created by calling the constructor `CartesianBox(args...)`
with one of the following forms:

```julia
CartesianBox(A)
CartesianBox(axes(A))
CartesianBox(CartesianIndex(imin,jmin,...), CartesianIndex(imax,jmax,...))
CartesianBox((imin:imax, jmin:jmax, ...))
CartesianBox((dim1, dim2, ...))
CartesianBox(R)
```

where `A` is an array, `imin`, `imax`, ... are integers, `dim1`, `dim2`,
... are integer dimension lenghts, `R` is an instance of `CartesianIndices` or,
on Julia ≤ 0.6, an instance of `CartesianRange`.  Note that ranges like `k:k`
can be abbreviated by just specifying `k`.  However beware that a tuple of
integers is interpreted as a list of dimensions similar to `(1:dim1, 1:dim2,
...)`.  There is the same ambiguity in the constructors of `CartesianIndices`
and of `CartesianRange`.

An instance of `CartesianBox` can be used in a loop as follows:

```julia
B = CartesianBox(...)
for i in B
   A[i] = ...
end
```

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.

The call:

```julia
intersection(R, S)
```

yields the bounding-box (as a `CartesianBox`) which is the intersection of the
two regions defined by `R` and `S`.  This method is similar to `intersect(R,S)
= R ∩ S` which yields an array of Cartesian indices and is **much** slower (and
hence less useful).

The call:

```julia
isnonemptypartof(R, S)
```

yields whether region `R` is nonempty and a valid part of `S`.  If `R` is a
Cartesian range, then `S` can be an array to check whether `R` is a valid
nonempty region of interest of `S`.

The call:

```julia
boundingbox([pred,] A [, B])
```

yields the bounding-box of values in array `A` for which the predicate function
`pred` is true.  If the predicate function `pred` is omitted, the result is the
bounding-box of non-zero values in array `A` or of the `true` values in `A` if
its elements are of type `Bool`.  Optional argument `B` is to only consider a
sub-region `B` of `A` (`B` can be a `CartesianBox`, a `CartesianIndices`, a
`CartesianRange` or a tuple of integer unit ranges).


## Restrictions

The algorithm for finding the bounding-box of valid values is pretty simple
and scales as `O(length(A))`.


## Installation

`CartesianBoxes.jl` is not yet an
[official Julia package](https://pkg.julialang.org/) so you have to clone the
repository to install the package:

```julia
Pkg.clone("https://github.com/emmt/CartesianBoxes.jl.git")
```

There is nothing to build so no needs to call `Pkg.build("CartesianBoxes")`.

Later, it is sufficient to do:

```julia
Pkg.update("CartesianBoxes")
```

to pull the latest version.


## Rationale

In the past, rectangular regions of `N`-dimensional indices were defined by
instances of `CartesianRange{CardinalIndex{N}}` in Julia and have a number of
related methods which make coding [multi-dimensional
algorithms](https://julialang.org/blog/2016/02/iteration) not only *easy* but
also *very efficient*.  Recent Julia versions (≥ 0.7) introduced a change in
the representation of such sets of multi-dimensional indices which are now
called
[`CartesianIndices{N}`](https://github.com/JuliaLang/julia/issues/20974).
There have been a few changes in the API but, in general, it is sufficient to
replace `CartesianRange{CardinalIndex{N}}` by `CartesianIndices{N}` in the
code.  For backward compatibility, [`using
Compat`](https://github.com/JuliaLang/Compat.jl) let you use
`CartesianIndices{N}` with Julia ≤ 0.6.  However, while the performances have
been maintained or even improved with `CartesianIndices{N}` in Julia ≥ 0.7
compared to `CartesianRange{CardinalIndex{N}}` in Julia ≤ 0.6, this is not true
if you are using `CartesianIndices{N}` in Julia 0.6 via the
[Compat](https://github.com/JuliaLang/Compat.jl) package.  For instance, I
measured (with [BenchmarkTools](http://github.com/JuliaCI/BenchmarkTools.jl))
slowdowns worse than a factor of 30 for simple additions of arrays.  My guess
is that this is because [Compat](https://github.com/JuliaLang/Compat.jl) does
not extend `simd_outer_range()`, `simd_inner_length()` nor `simd_index()`
methods for `CartesianIndices{N}`.

Another issue is that, I was not quite satisfied with the API and wanted to add
functionalities but did not want to interfere with those implemented by
`CartesianIndices` or `CartesianRange`.
