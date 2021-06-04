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
    intersect, issubset, getindex, setindex!, fill!, show
import Base: simd_outer_range, simd_inner_length, simd_index

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

`CartesianBox(size(A))` yields a Cartesian region of the same size as `A` but
whose indices all start at 1.  In most cases, `CartesianBox(axes(A))` is more
likely to do the right thing.  In fact, `CartesianBox(A)` is equivalent to
`CartesianBox(axes(A))`.

An instance of `CartesianBox` can also be constructed from an instance of
`CartesianIndices` and conversely.  The conversion is lossless in the sense of
the following example:

    B = CartesianBox(...)
    R = CartesianIndices(B)
    CartesianBox(R) === B # always true

An instance of `CartesianBox` can be efficiently used in a loop as follows:

    for i in CartesianBox(...)
       A[i] = ...
    end

where `i` will be set to a `CartesianIndex` with all the multi-dimensional
indices of the rectangular region defined by `B`.

When at least one of `A` or `B` is a Cartesian box, the expression `A ∩ B`, or
equivalently `intersect(A,B)`, yields the Cartesian box contining all indices
in `A` and `B`.  This may be used to write safe loops like:

    A = ...               # some array
    B = CartesianBox(...) # some region of interest
    @inbounds for i in B ∩ A
        A[i] = ...
    end

to operate on the indices of the Cartesian box `B` that are valid for `A`.

When at least one of `A` or `B` is a Cartesian box, the expression `A ⊆ B`, or
equivalently `intersect(A,B)`, yields whether all Cartesian indices defined by
`A` are also indices of `B`.  This may be used as:

    A = ...               # some array
    B = CartesianBox(...) # some region of interest
    if B ⊆ A
        @inbounds for i in B
            A[i] = ...
        end
    end

to only access the indices of the Cartesian box `B` if they are all valid for
`A`.

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

isempty(B::CartesianBox) = isempty_(first(B), last(B))
@inline isempty_(first::CartesianIndex{N}, last::CartesianIndex{N}) where {N} =
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
    intersection(A, B)

yields the Cartesian box given by the intersection of the Cartesian regions
defined by `A` and `B`.  In this context, a Cartesian region can specified by a
Cartesian box, a list of integer valued ranges, a list of dimensions, or an
instance of `CartesianIndices`.

This method is equivalent to `intersect(A,B)`, or `A ∩ B` for short, when at
least one of `A` or `B` is a Cartesian box.

See also: [`CartesianBox`](@ref), [`isnonemptypartof`](@ref).

"""
intersection(A::CartesianBox{N}, B::CartesianBox{N}) where {N} = A ∩ B
intersection(A::CartesianBox{N}, B::CartesianBoxable{N}) where {N} =
    A ∩ CartesianBox(B)
intersection(A::CartesianBoxable{N}, B::CartesianBox{N}) where {N} =
    CartesianBox(A) ∩ B
intersection(A::CartesianBoxable{N}, B::CartesianBoxable{N}) where {N} =
    CartesianBox(A) ∩ CartesianBox(B)

# Override ∩(A,B) and ⊆(A,B) when at least one of A or B is a Cartesian box,
# the other being boxable or an abstract array.
for f in (:intersect, :issubset)
    @eval begin
        $f(A::CartesianBox{N}, B::CartesianBoxable{N}) where {N} =
            $f(A, CartesianBox(B))
        $f(A::CartesianBox{N}, B::AbstractArray{<:Any,N}) where {N} =
            $f(A, CartesianBox(B))
        $f(A::CartesianBoxable{N}, B::CartesianBox{N}) where {N} =
            $f(CartesianBox(A), B)
        $f(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
            $f(CartesianBox(A), B)
    end
end
@inline intersect(A::CartesianBox{N}, B::CartesianBox{N}) where {N} =
    CartesianBox(max(first(A), first(B)), min(last(A), last(B)))
@inline issubset(A::CartesianBox{N}, B::CartesianBox{N}) where {N} =
    isempty(A) || isnonemptypartof(A,B)

"""
    isnonemptypartof(A, B)

yields whether the region defined by `A` is nonempty and a valid part of the
region defined by `B` or of the contents of `B` if it is an array.  If this
method returns `false`, you may call `isempty(A)` to check whether `A` is
empty.

When at least one of `A` or `B` is a Cartesian box, the expression `A ⊆ B` or
`issubset(A,B)` is equivalent to:

    isempty(A) || isnonemptypartof(A, B)

See also: [`CartesianBox`](@ref), [`intersection`](@ref).

"""
function isnonemptypartof(A::Union{CartesianBoxable{N},CartesianBox{N}},
                          B::Union{CartesianBoxable{N},
                                   AbstractArray{<:Any,N}}) where {N}
    isnonemptypartof(A, CartesianBox(B))
end

isnonemptypartof(A::CartesianBoxable{N}, B::CartesianBox{N}) where {N} =
    isnonemptypartof(CartesianBox(A), B)

@inline isnonemptypartof(A::CartesianBox{N}, B::CartesianBox{N}) where {N} =
    first(B) ≤ first(A) ≤ last(A) ≤ last(B)


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
    Imin = typemax_(CartesianIndex{N})
    Imax = typemin_(CartesianIndex{N})
    @inbounds for I in CartesianBox(A)
        if pred(A[I])
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    if isempty_(Imin, Imax)
        Imin = one_(CartesianIndex{N})
        Imax = zero(CartesianIndex{N})
    end
    return CartesianBox(Imin, Imax)
end

function boundingbox(pred,
                     A::AbstractArray{T,N},
                     B::CartesianBox{N}) where {T,N}
    Imin = typemax_(CartesianIndex{N})
    Imax = typemin_(CartesianIndex{N})
    @inbounds for I in CartesianBox(A) ∩ B
        if pred(A[I])
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    if isempty_(Imin, Imax)
        Imin = one_(CartesianIndex{N})
        Imax = zero(CartesianIndex{N})
    end
    return CartesianBox(Imin, Imax)
end

@inline typemin_(::Type{<:CartesianIndex{N}}) where {N} =
    CartesianIndex(ntuple(i -> typemin(Int), Val(N)))
@inline typemax_(::Type{<:CartesianIndex{N}}) where {N} =
    CartesianIndex(ntuple(i -> typemax(Int), Val(N)))

one_(I::CartesianIndex) = one_(typeof(I))
one_(T::Type{<:CartesianIndex}) =
    @static if VERSION < v"1.1.0-rc1"
        one(T)
    else
        oneunit(T)
    end

"""
    isnonzero(x)

yields whether `x` is non-zero, the negation of `iszero(x)`.

"""
isnonzero(x) = !iszero(x)

end # module
