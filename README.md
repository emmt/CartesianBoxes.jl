# Flexible and efficicient multi-dimensional index boxes in Julia

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/CartesianBoxes.jl.svg?branch=master)](https://travis-ci.org/emmt/CartesianBoxes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/CartesianBoxes.jl?branch=master)](https://ci.appveyor.com/project/emmt/CartesianBoxes-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/emmt/CartesianBoxes.jl/badge.svg?branch=master)](https://coveralls.io/github/emmt/CartesianBoxes.jl?branch=master)
[![codecov.io](http://codecov.io/github/emmt/CartesianBoxes.jl/coverage.svg?branch=master)](http://codecov.io/github/emmt/CartesianBoxes.jl?branch=master)

This module implements `CartesianBox{N}` to define rectangular regions of
`N`-dimensional indices in Julia arrays.  Cartesian boxes are similar to
`CartesianIndices` or, for Julia version ≤ 0.6, to `CartesianRange` but, being
a different type, they can be used to specifically extend methods without
introducing unexpected behaviors in other Julia modules.

For multi-dimensional loops, instances of `CartesianBox{N}` are as fast as
`CartesianIndices{N}` or as `CartesianRange{CardinalIndex{N}}`.  They can thus
be used as a *fast* and *portable* replacement (see [*Rationale*](#rationale)
below) of these different representations of rectangular multi-dimensional
regions.


## Usage

### Construction

A Cartesian box is created by calling the constructor `CartesianBox(args...)`
with a variety of arguments.  For instance:

```julia
CartesianBox(A)
```

yields the Cartesian box which contains all indices of array `A`.  An arbitrary
region whose first an last multi-dimensional indices are `(imin,jmin,...)` and
`(imax,jmax,...)` can be defined by one of the following calls:

```julia
CartesianBox(CartesianIndex(imin,jmin,...), CartesianIndex(imax,jmax,...))
CartesianBox((imin:imax, jmin:jmax, ...))
```

hence

```julia
CartesianBox(axes(A))
```

also defines a Cartesian box encompassing all indices of array `A`.  For normal
arrays (with 1-based indices), it is sufficient to provide the dimensions of
the array:

```julia
CartesianBox(size(A))
CartesianBox((dim1, dim2, ...))
```

Note that index ranges like `k:k` can be abbreviated by just specifying `k`.
However, a tuple of integers, say `(n1,n2,...)`, is interpreted as a list of
dimensions, as if `(1:n1, 1:n2, ...)` as been specified.  There is the same
ambiguity in the constructors of `CartesianIndices` and of `CartesianRange`.

Finally, it is possible to convert an instance, say `R`, of `CartesianIndices`
or an instance of `CartesianRange` into a `CartesianBox`:

```julia
B = CartesianBox(R)
```

The reverse operation is also possible, `CartesianIndices(B)` and
`CartesianRange(B)` work as expected.


### Fast iterations

An instance of `CartesianBox` can be used in a loop as follows:

```julia
B = CartesianBox(...)
for i in B
   A[i] = ...
end
```

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.  To benefit from faster loops
you may suppress bound checking and activate
[SIMD](https://fr.wikipedia.org/wiki/Single_instruction_multiple_data)
vectorization:

```julia
B = CartesianBox(...)
@inbounds @simd for i in B
   A[i] = ...
end
```

### Indexation and views

You may extract the region defined by a Cartesian box `B` from an array `A`
by calling:

```julia
C = A[B]
```

Setting values is also possible with

```julia
A[B] = C
```

where `C` is an array of same dimensions as the region defined by `B`.  To fill
the region `B` of array `A` with a scalar `x`, just do:

```julia
fill!(A, B, x) -> A
```

A *view* can be created by:

```julia
V = view(A, B)
```

which yields a sub-array `V` sharing its elements with `A` in the region
defined by `B`.


### Methods

The call:

```julia
intersection(R, S)
```

yields the Cartesian box which is the intersection of the Cartesian regions
defined by `R` and `S`.  This method is similar to `intersect(R,S) = R ∩ S`
which yields an array of Cartesian indices and is **much** slower (and hence
less useful).

The call:

```julia
isnonemptypartof(R, S)
```

yields whether the region defined by `R` is nonempty and a valid part of the
region defined by `S` or of the contents of `S` if it is an array.  If this
method returns `false`, you may call:

```julia
isempty(R)
```

to check whether `R` was empty.

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

* The algorithm for finding the bounding-box of valid values is pretty simple
  and scales as `O(length(A))`.

* There is no way to define an *empty* Cartesian box when `N=0`.


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

Another motivation for this module, was that I wanted to add some
functionalities in a such a way that is does not interfere with how
`CartesianIndices` or `CartesianRange` are used by others.
