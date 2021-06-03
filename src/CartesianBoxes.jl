#
# CartesianBoxes.jl -
#
# Extends CartesianIndices.
#
#-------------------------------------------------------------------------------
#
# This file is part of the `CartesianBoxes.jl` package which is licensed under
# the MIT "Expat" License.
#
# Copyright (c) 2017-2021 Éric Thiébaut.
#

__precompile__(true)

module CartesianBoxes

export
    CartesianBox,
    boundingbox,
    intersection,
    isnonemptypartof,
    isnonzero,
    ranges

using Base: tail, @propagate_inbounds

import Base: CartesianIndices, IndexStyle, axes, eachindex, isempty,
    iterate, first, last, ndims, eltype, length, size, view,
    getindex, setindex!, fill!, show
import Base: simd_outer_range, simd_inner_length, simd_index

# Deal with compatibility issues.
using Compat

"""

`CartesianBox{N}` defines a rectangular region of `N`-dimensional indices and
can be constructed by:

    CartesianBox(A)
    CartesianBox(axes(A))
    CartesianBox(CartesianIndex(imin,jmin,...), CartesianIndex(imax,jmax,...))
    CartesianBox((imin,jmin,...), (imax,jmax,...))
    CartesianBox((imin:imax, jmin:jmax, ...))

where `A` is an array (to define a region consisting in all the indices in the
array), `imin`, `imax`, etc. are integers (to define a region from
`(imin,jmin,...)` to `(imax,jmax,...)`.

An instance of `CartesianBox` can also be constructed from an instance of
`CartesianIndices` and conversely:

    B = CartesianBox(...)
    R = CartesianIndices(B)
    CartesianBox(R) === B # yields true

An instance of `CartesianBox` can be used in a loop as follows:

    B = CartesianBox(...)
    for i in B
       A[i] = ...
    end

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.

See also: [`boundingbox`](@ref), `CartesianIndices`, `CartesianIndex`,
[`intersection`](@ref).

"""
struct CartesianBox{N,R<:CartesianIndices{N}} <: AbstractArray{CartesianIndex{N},N}
    inds::R
end
CartesianBox(B::CartesianBox) = B
CartesianBox(A::AbstractArray) = CartesianBox(axes(A))
CartesianBox(inds::Tuple{Vararg{Union{<:Integer,AbstractUnitRange{<:Integer}}}}) =
        CartesianBox(CartesianIndices(inds))
#CartesianBox(::Tuple{}) = CartesianBox(CartesianIndices(()))
CartesianBox(first::CartesianIndex{N}, last::CartesianIndex{N}) where {N} =
    CartesianBox(first.I, last.I)
CartesianBox(first::NTuple{N,Integer}, last::NTuple{N,Integer}) where {N} =
    CartesianBox(map((i,j) -> i:j, first, last))

CartesianIndices(B::CartesianBox) = B.inds

"""
    ranges(B) -> inds

yields the `N`-tuple of index ranges in the Cartesian box `B` (an instance of
[`CartesianBox{N}`](@ref)).

"""
ranges(B::CartesianBox) = CartesianIndices(B).indices

first(B::CartesianBox) = first(CartesianIndices(B))
last(B::CartesianBox) = last(CartesianIndices(B))
size(B::CartesianBox) = size(CartesianIndices(B))
axes(B::CartesianBox) = axes(CartesianIndices(B))

isempty(B::CartesianBox) = _isempty(first(B), last(B))
@inline _isempty(first::CartesianIndex{N}, last::CartesianIndex{N}) where {N} =
    any(map(isless, last.I, first.I))

show(io::IO, ::MIME"text/plain", B::CartesianBox{N}) where {N} = begin
    print(io, "CartesianBox{",N,"}((")
    if N >= 1
        I = first(B).I
        J = last(B).I
        print(io, I[1], ":", J[1])
        for k in 2:N
            print(io, ", ", I[k], ":", J[k])
        end
    end
    print(io, "))")
end

view(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
    view(A, ranges(B)...)

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
    A[ranges(B)...] = x

fill!(A::AbstractArray{T,N}, B::CartesianBox{N}, x) where {T,N} =
    fill!(A, B, convert(T, x)::T)
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

IndexStyle(::Type{<:CartesianBox{N,R}}) where {N,R} = IndexStyle(R)
@inline @propagate_inbounds getindex(B::CartesianBox, I...) =
    getindex(CartesianIndices(B), I...)

"""

`CartesianBoxable{N}` is a union of types (other than `CartesianBox{N}`) which
can be automatically converted into a `CartesianBox{N}`.

!!! note
    Although the constructor `CartesianBox` can be also applied to any instance
    of `AbstractArray`, this abstract type does not belong to the union
    `CartesianBoxable` as it is considered that such a conversion cannot be
    automatic.

"""
const CartesianBoxable{N} = Union{CartesianIndices{N},
                                  NTuple{N,AbstractUnitRange{<:Integer}},
                                  NTuple{N,Integer}}


# Extend eachindex() method.
eachindex(::IndexCartesian, B::CartesianBox) = B

# Make CartesianBox iterable.
iterate(iter::CartesianBox) = iterate(CartesianIndices(iter))
iterate(iter::CartesianBox, state) = iterate(CartesianIndices(iter), state)

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
    intersection(R, S)

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
    isnonemptypartof(R, S)

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
    boundingbox([pred,] A [, B])

yields the bounding-box of the elements in array `A` for which the predicate
function `pred` is true.  If the predicate function `pred` is omitted, the
result is the bounding-box of non-zero values in array `A` or of the `true`
values in `A` if its elements are of type `Bool`.  Optional argument `B` is to
only consider a sub-region `B` of `A` (`B` can be a `CartesianBox`, a
`CartesianIndices`, or a tuple of integer valued unit ranges).

!!! warning
    The algorithm is pretty silly for now and could be made faster than
    `O(length(A))`.

See also: [`CartesianBox`](@ref), [`intersection`](@ref), [`isnonzero`](@ref).

"""
boundingbox(A::AbstractArray{Bool}) = boundingbox(identity, A)
boundingbox(A::AbstractArray) = boundingbox(isnonzero, A)
boundingbox(A::AbstractArray{Bool,N}, B::CartesianBox{N}) where {N} =
    boundingbox(identity, A, B)
boundingbox(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
    boundingbox(isnonzero, A, B)
boundingbox(A::AbstractArray{<:Any,N}, B::CartesianBoxable{N}) where {N} =
    boundingbox(A, CartesianBox(B))
boundingbox(pred, A::AbstractArray{<:Any,N}, B::CartesianBoxable{N}) where {N} =
    boundingbox(pred, A, CartesianBox(B))

function boundingbox(pred,
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
    if _isempty(Imin, Imax)
        Imin = oneunit(CartesianIndex{N})
        Imax = zero(CartesianIndex{N})
    end
    return CartesianBox(Imin, Imax)
end

function boundingbox(pred,
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
    if _isempty(Imin, Imax)
        Imin = oneunit(CartesianIndex{N})
        Imax = zero(CartesianIndex{N})
    end
    return CartesianBox(Imin, Imax)
end

"""
    isnonzero(x)

yields whether `x` is non-zero, the negation of `iszero(x)`.

"""
isnonzero(x) = !iszero(x)

end # module
