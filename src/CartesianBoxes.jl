#
# CartesianBoxes.jl -
#
# Extends CartesianIndices (or CartesianRange on Julia ≤ 0.6).
#
#-------------------------------------------------------------------------------
#
# This file is part of the `CartesianBoxes.jl` package which is licensed under
# the MIT "Expat" License.
#
# Copyright (c) 2017-2018 Éric Thiébaut.
#

__precompile__(true)

module CartesianBoxes

export
    CartesianBox,
    boundingbox,
    intersection,
    isnonemptypartof,
    isnonzero

using Base: tail
using Base.IteratorsMD: inc

import Base: eachindex, isempty, first, last, ndims, eltype, length, size,
    view, getindex, setindex!, fill!
import Base: simd_outer_range, simd_inner_length, simd_index

# Deal with compatibility issues.
using Compat
@static if isdefined(Base, :CartesianIndices)
    import Base: CartesianIndices
else
    import Base: CartesianRange
    import Compat: CartesianIndices
end
@static if isdefined(Base, :axes)
    import Base: axes
else
    import Base: indices
    const axes = indices
end

"""

`CartesianBox{N}` defines a rectangular region of `N`-dimensional indices and
can be constructed by:

```julia
CartesianBox(A)
CartesianBox(axes(A))
CartesianBox(CartesianIndex(imin,jmin,...), CartesianIndex(imax,jmax,...))
CartesianBox((imin:imax, jmin:jmax, ...))
```

where `A` is an array (to define a region consisting in all the indices in the
array), `imin`, `imax`, etc. are integers (to define a region from
`(imin,jmin,...)` to `(imax,jmax,...)`.

An instance of `CartesianBox` can also be constructed from an instance of
`CartesianIndices` and, for Julia version ≤ 0.6, from an instance of
`CartesianRange`.

An instance of `CartesianBox` can be used in a loop as follows:

```julia
B = CartesianBox(...)
for i in B
   A[i] = ...
end
```

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.

See also: [`boundingbox`](@ref), [`CartesianIndices`](@ref),
[`CartesianRange`](@ref), [`CartesianIndex`](@ref), [`intersection`](@ref).

"""
struct CartesianBox{N} <: AbstractArray{CartesianIndex{N},N}
    first::CartesianIndex{N}
    last::CartesianIndex{N}
end
CartesianBox(B::CartesianBox) = B
CartesianBox(A::AbstractArray) = CartesianBox(axes(A))
CartesianBox(R::CartesianIndices) = CartesianBox(R.indices)
@static if ! isdefined(Base, :CartesianIndices)
    CartesianBox(R::CartesianRange) = CartesianBox(R.start, R.stop)
end
CartesianBox(dims::NTuple{N,Integer}) where N =
    CartesianBox(one(CartesianIndex{N}), CartesianIndex(dims))
CartesianBox(rngs::NTuple{N,AbstractUnitRange{<:Integer}}) where N =
    CartesianBox(CartesianIndex(map(first, rngs)),
                 CartesianIndex(map(last,  rngs)))
CartesianBox(first::NTuple{N,Integer}, last::NTuple{N,Integer}) where {N} =
    CartesianRange(CartesianIndex(first), CartesianIndex(last))
CartesianBox(::Tuple{}) =
    CartesianBox{0}(CartesianIndex{0}(()), CartesianIndex{0}(()))

@inline first(B::CartesianBox) = B.first
@inline last(B::CartesianBox) = B.last
@inline isempty(B::CartesianBox) = any(map(>, first(B).I, last(B).I))
ndims(::CartesianBox{N}) where {N} = N
eltype(::CartesianBox{N}) where {N} = CartesianIndex{N}
length(B::CartesianBox) = prod(size(B))
size(B::CartesianBox) = map(_dimensionlength, first(B).I, last(B).I)
size(B::CartesianBox, d) =
    _dimensionlength(first(B).I[d], last(B).I[d])
axes(B::CartesianBox) = map((i,j) -> i:j, first(B).I, last(B).I)
axes(B::CartesianBox, d) = (first(B).I[d]:last(B).I[d])

_dimensionlength(first::Int, last::Int) = max(0, last - first + 1)


view(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
    view(A, map((i,j) -> i:j, first(B).I, last(B).I)...)

function getindex(A::AbstractArray{T,N}, B::CartesianBox{N}) where {T,N}
    empty = isempty(B)
    empty || isnonemptypartof(B, A) ||
        throw(BoundsError("sub-region is not part of array"))
    C = Array{T}(undef, size(B))
    if ! empty
        if any(i -> i != 1, first(B).I)
            off = CartesianIndex(map(i -> i - 1, first(B).I))
            @inbounds @simd for i in B
                C[i - off] = A[i]
            end
        else
            @inbounds @simd for i in B
                C[i] = A[i]
            end
        end
    end
    return C
end

setindex!(A::AbstractArray{<:Any,N}, x, B::CartesianBox{N}) where {N} =
    A[map((i,j) -> i:j, first(B).I, last(B).I)...] = x

fill!(A::AbstractArray{T,N}, B::CartesianBox{N}, x) where {T,N} =
    fill!(A, B, convert(T, x))
function fill!(A::AbstractArray{T,N},
               B::CartesianBox{N}, x::T) where {T,N}
    if ! isempty(B)
        isnonemptypartof(B, A) ||
            throw(BoundsError("sub-region is not part of array"))
        @inbounds @simd for i in B
            A[i] = x
        end
    end
    return A
end


"""

`CartesianBoxable{N}` is a union of types (other than `CartesianBox{N}`) which
can be automatically converted into a `CartesianBox{N}`.

*Remarks:* Although the constructor `CartesianBox` can be also applied to any
instance of `AbstractArray`, this abstract type does not belong to the union
`CartesianBoxable` as it is considered that such a conversion cannot be
automatic.

"""
CartesianBoxable
@static if isdefined(Base, :CartesianIndices)
    # Union of argument types accepted by CartesianBox outer constructors.
    const CartesianBoxable{N} = Union{CartesianIndices{N},NTuple{N,AbstractUnitRange{<:Integer}},NTuple{N,Integer}}
else
    # Idem but with CartesianRange.
    const CartesianBoxable{N} = Union{CartesianIndices{N},NTuple{N,AbstractUnitRange{<:Integer}},NTuple{N,Integer},CartesianRange{CartesianIndex{N}}}
end


# Extend eachindex() method and provide converter to CartesianIndices (and
# CartesianRange for Julia ≤ 0.6).
eachindex(::IndexCartesian, B::CartesianBox) = B
CartesianIndices(B::CartesianBox) =
    CartesianIndices(map((i,j) -> i:j, first(B).I, last(B).I))
@static if ! isdefined(Base, :CartesianIndices)
    CartesianRange(B::CartesianBox) = CartesianRange(first(B), last(B))
end

# Make CartesianBox iterable.
@static if isdefined(Base, :iterate)
    import Base: iterate
    @inline iterate(iter::CartesianBox) =
        isempty(iter) ? nothing : (first(iter), first(iter))
    @inline function iterate(iter::CartesianBox, state)
        next = CartesianIndex(inc(state.I, first(iter).I, last(iter).I))
        next.I[end] > last(iter).I[end] ? nothing : (next, next)
    end
else
    import Base: start, done, next
    @inline start(iter::CartesianBox) =
        isempty(iter) ? last(iter) + 1 : first(iter)
    @inline next(iter::CartesianBox{N}, state) where {N} =
        state, CartesianIndex{N}(inc(state.I, first(iter).I, last(iter).I))
    @inline done(iter::CartesianBox{N}, state) where {N} =
        state.I[end] > last(iter).I[end]
end

# 0-d cartesian ranges are special-cased to iterate once and only once
@static if isdefined(Base, :iterate)
    iterate(iter::CartesianBox{0}, done=false) =
        done ? nothing : (CartesianIndex(), true)
else
    start(iter::CartesianBox{0}) = false
    next(iter::CartesianBox{0}, state) = first(iter), true
    done(iter::CartesianBox{0}, state) = state
end

# Extend methods for fast SIMD iterations.
simd_outer_range(iter::CartesianBox{0}) = iter
simd_outer_range(iter::CartesianBox) =
    CartesianBox(CartesianIndex(tail(first(iter).I)),
                 CartesianIndex(tail(last(iter).I)))
simd_inner_length(iter::CartesianBox{0}, ::CartesianIndex) = 1
simd_inner_length(iter::CartesianBox, I::CartesianIndex) =
    last(iter).I[1] - first(iter).I[1] + 1
simd_index(iter::CartesianBox{0}, ::CartesianIndex, I1::Int) = first(iter)
@inline simd_index(iter::CartesianBox, Ilast::CartesianIndex, I1::Int) =
    CartesianIndex((I1 + first(iter)[1], Ilast.I...))

"""

```julia
intersection(R, S)
```

yields the Cartesian box which is the intersection of the Cartesian regions
defined by `R` and `S`.  This method is similar to `intersect(R,S) = R ∩ S`
which yields an array of Cartesian indices and is **much** slower (and less
useful).

See also: [`CartesianBox`](@ref), [`isnonemptypartof`](@ref).

"""
intersection(R::CartesianBox{N}, S::CartesianBox{N}) where {N} =
    CartesianBox(max(first(R), first(S)), min(last(R), last(S)))
intersection(R::CartesianBox{N}, S::CartesianBoxable{N}) where {N} =
    intersection(R, CartesianBox(S))
intersection(R::CartesianBoxable{N}, S::CartesianBox{N}) where {N} =
    intersection(CartesianBox(R), S)
intersection(R::CartesianBoxable{N}, S::CartesianBoxable{N}) where {N} =
    intersection(CartesianBox(R), CartesianBox(S))

"""
```julia
isnonemptypartof(R, S)
```

yields whether the region defined by `R` is nonempty and a valid part of the
region defined by `S` or of the contents of `S` if it is an array.  If this
method returns `false`, you may call `isempty(R)` to check whether `R` was
empty.

See also: [`CartesianBox`](@ref), [`intersection`](@ref).

"""
function isnonemptypartof(R::Union{CartesianBoxable{N},CartesianBox{N}},
                          S::Union{CartesianBoxable{N},
                                   AbstractArray{<:Any,N}}) where {N}
    isnonemptypartof(R, CartesianBox(S))
end

isnonemptypartof(R::CartesianBoxable{N}, S::CartesianBox{N}) where {N} =
    isnonemptypartof(CartesianBox(R), S)

isnonemptypartof(R::CartesianBox{N}, S::CartesianBox{N}) where {N} =
    first(S) ≤ first(R) ≤ last(R) ≤ last(S)


"""
```julia
boundingbox([pred,] A [, B])
```

yields the bounding-box of the elements in array `A` for which the predicate
function `pred` is true.  If the predicate function `pred` is omitted, the
result is the bounding-box of non-zero values in array `A` or of the `true`
values in `A` if its elements are of type `Bool`.  Optional argument `B` is to
only consider a sub-region `B` of `A` (`B` can be a `CartesianBox`, a
`CartesianIndices`, a `CartesianRange` or a tuple of integer valued unit
ranges).

FIXME: The algorithm is pretty silly for now and could be made faster than
       `O(length(A))`.

See also: [`CartesianBox`](@ref), [`intersection`](@ref), [`isnonzero`](@ref).

"""
boundingbox(A::AbstractArray{Bool,N}) where {N} = boundingbox(identity, A)
boundingbox(A::AbstractArray{<:Any,N}) where {N} = boundingbox(isnonzero, A)
boundingbox(A::AbstractArray{Bool,N}, B::CartesianBox{N}) where {N} =
    boundingbox(identity, A, B)
boundingbox(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
    boundingbox(isnonzero, A, B)
boundingbox(A::AbstractArray{<:Any,N}, B::CartesianBoxable{N}) where {N} =
    boundingbox(A, CartesianBox(B))
function boundingbox(pred::Function, A::AbstractArray{<:Any,N},
                     B::CartesianBoxable{N}) where {N}
    boundingbox(pred, A, CartesianBox(B))
end

function boundingbox(pred::Function,
                     A::AbstractArray{T,N}) where {T,N}
    inds = axes(A)
    Imin = CartesianIndex(map(r -> last(r) + 1, inds))
    Imax = CartesianIndex(map(r -> first(r) - 1, inds))
    @inbounds for I in CartesianBox(A)
        if pred(A[I])
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return CartesianBox(Imin, Imax)
end

function boundingbox(pred::Function,
                     A::AbstractArray{T,N},
                     B::CartesianBox{N}) where {T,N}
    R = intersection(CartesianBox(A), B)
    inds = axes(R)
    Imin = CartesianIndex(map(r -> last(r) + 1, inds))
    Imax = CartesianIndex(map(r -> first(r) - 1, inds))
    @inbounds for I in R
        if pred(A[I])
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return CartesianBox(Imin, Imax)
end

"""

```julia
isnonzero(x)
```

yields whether `x` is non-zero, the negation of `iszero(x)`.

See also [`iszero`](@ref).

"""
isnonzero(x) = !iszero(x)

end # module
