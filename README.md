# Flexible and efficicient multi-dimensional index boxes in Julia

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/CartesianBoxes.jl.svg?branch=master)](https://travis-ci.org/emmt/CartesianBoxes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/CartesianBoxes.jl?branch=master)](https://ci.appveyor.com/project/emmt/CartesianBoxes-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/emmt/CartesianBoxes.jl/badge.svg?branch=master)](https://coveralls.io/github/emmt/CartesianBoxes.jl?branch=master)
[![codecov.io](http://codecov.io/github/emmt/CartesianBoxes.jl/coverage.svg?branch=master)](http://codecov.io/github/emmt/CartesianBoxes.jl?branch=master)

This module implements `CartesianBox{N}` to define rectangular regions of
`N`-dimensional indices in Julia arrays.  Cartesian boxes are similar to
`CartesianIndices` but, being a different type, they can be used to
specifically extend methods without introducing unexpected behaviors in other
Julia modules.

For multi-dimensional loops, instances of `CartesianBox{N}` are as fast as
`CartesianIndices{N}`.  They can thus be used as a *fast* and *portable*
replacement (see [*Rationale*](#rationale) below) of these different
representations of rectangular multi-dimensional regions.


## Usage

### Construction

A Cartesian box is created by calling the constructor `CartesianBox(args...)`
with a variety of arguments.  For instance:

```julia
CartesianBox(A)
```

yields the Cartesian box which contains all indices of array `A`.  An arbitrary
region whose first an last multi-dimensional indices are `(imin,jmin,...)` and
`(imax,jmax,...)` can be defined by one of:

```julia
CartesianBox(CartesianIndex(imin,jmin,...), CartesianIndex(imax,jmax,...))
CartesianBox((imin:imax, jmin:jmax, ...))
CartesianBox((imin,jmin,...), (imax,jmax,...))
```

hence

```julia
CartesianBox(axes(A))
```

also defines a Cartesian box containing all indices of array `A`.  For normal
arrays (with 1-based indices), it is sufficient to provide the dimensions of
the array:

```julia
CartesianBox(size(A))
CartesianBox((dim1, dim2, ...))
```

This is however not recommended, `CartesianBox(axes(A))` or, for short,
`CartesianBox(A)` are more likely to be coorect for any kind of array `A`.

It is also possible to convert an instance `R` of `CartesianIndices` into a
`CartesianBox` by calling the constructor:

```julia
B = CartesianBox(R)
```

The reverse operation is also possible and is lossless as shown by:

```julia
CartesianIndices(B) === R
```

which is always true.

To retrieve the `N`-tuple of ranges that constitute a Cartesian box `B`, call
`Tuple(B)`.  This is not the same as `axes(B)` which yields the ranges to index
`B` itself.


```julia
B = CartesianBox(2:7, 3:5)
Tuple(B) -> (2:7, 3:5)
axes(B) -> (Base.OneTo(6),Base.OneTo(3))
```


### Fast (and safe) iterations

An instance of `CartesianBox` can be used in a loop as follows:

```julia
for i in CartesianBox(...)
   A[i] = ...
end
```

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.  To benefit from faster loops
you may suppress bound checking and activate
[SIMD](https://fr.wikipedia.org/wiki/Single_instruction_multiple_data)
vectorization:

```julia
@inbounds @simd for i in CartesianBox(...)
   A[i] = ...
end
```

When at least one of `A` or `B` is a Cartesian box, the expression `A ∩ B`, or
equivalently `intersect(A,B)`, yields the Cartesian box contining all indices
in `A` and `B`.  This may be used to write safe loops like:

```julia
A = ...               # some array
B = CartesianBox(...) # some region of interest
@inbounds for i in B ∩ A
    A[i] = ...
end
```

to operate on the indices of the Cartesian box `B` that are valid for `A`.

When at least one of `A` or `B` is a Cartesian box, the expression `A ⊆ B`, or
equivalently `intersect(A,B)`, yields whether all Cartesian indices defined by
`A` are also indices of `B`.  This may be used as:

```julia
A = ...               # some array
B = CartesianBox(...) # some region of interest
if B ⊆ A
    @inbounds for i in B
        A[i] = ...
    end
end
```

to only access the indices of the Cartesian box `B` if they are all valid for
`A`.


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

### Shifting a CartesianBox

Operations `B + I` and `B - I` can be used to shift `B`, an instance of
`CartesianBox{N}`, by offset `I` specified as an instance of
`CartesianIndex{N}` or as an `N`-tuple of integers.


### Exported or extended methods

The call:

```julia
intersect(CartesianBox, A, B)
```

yields the Cartesian box given by the intersection of the Cartesian regions
defined by `A` and `B`.  In this context, a Cartesian region can specified by a
Cartesian box, a list of integer valued ranges, a list of dimensions, or an
instance of `CartesianIndices`.  This method is equivalent to `intersect(A,B)`,
or `A ∩ B` for short, when at least one of `A` or `B` is a Cartesian box.

The call:

```julia
isnonemptypartof(A, B)
```

yields whether the region defined by `A` is nonempty and a valid part of the
region defined by `B` or of the contents of `B` if `B` is an array.  This is
equivalent to:

```julia
!isempty(CartesianBox(A)) && (CartesianBox(A) ⊆ CartesianBox(B))
```

except that `A` may not be an array.

The call:

```julia
boundingbox([pred,] A [, B])
```

yields the bounding-box of values in array `A` for which the predicate function
`pred` is true.  If the predicate function `pred` is omitted, the result is the
bounding-box of non-zero values in array `A` or of the `true` values in `A` if
its elements are of type `Bool`.  Optional argument `B` is to only consider a
sub-region `B` of `A` (`B` can be a `CartesianBox`, a `CartesianIndices`, or a
tuple of integer valued unit ranges).


## Restrictions

* The algorithm for finding the bounding-box of valid values is pretty simple
  and scales as `O(length(A))`.

* There is no way to define an *empty* Cartesian box when `N=0`.


## Installation

The easiest way to install `CartesianBoxes` is via Julia registry
[`EmmtRegistry`](https://github.com/emmt/EmmtRegistry):

```julia
using Pkg
pkg"registry add General" # if you have not yet any registries
pkg"registry add https://github.com/emmt/EmmtRegistry"
pkg"add CartesianBoxes"
```


## Rationale

For pre-0.7 Julia versions, rectangular regions of `N`-dimensional indices were
defined by instances of `CartesianRange{CardinalIndex{N}}` in Julia and have a
number of related methods which make coding [multi-dimensional
algorithms](https://julialang.org/blog/2016/02/iteration) not only *easy* but
also *very efficient*.  More recent Julia versions (≥ 0.7) introduced a change
in the representation of such sets of multi-dimensional indices which are now
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
