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
# Copyright (c) 2017-2022 Éric Thiébaut.
#

__precompile__(true)

"""
`CartesianBoxes` extends `CartesianIndices`.
"""
module CartesianBoxes

export
    CartesianBox,
    boundingbox,
    isnonemptypartof,
    isnonzero

using Base: tail, @propagate_inbounds

import Base:
    intersect, issubset,
    simd_outer_range, simd_inner_length, simd_index

# Starting with Julia-1.6.0-beta1, Cartesian indices can have non-unit step.
const IndexRange{I} = @static if VERSION < v"1.6.0-beta1"
    AbstractUnitRange{I}
else
    OrdinalRange{I,I}
end

"""

`CartesianBox{N}` defines a rectangular region of `N`-dimensional indices and
can be constructed by:

    CartesianBox(A)
    CartesianBox(axes(A))
    CartesianBox(CartesianIndex(istart,jstart,...), CartesianIndex(istop,jstop,...))
    CartesianBox((istart,jstart,...), (istop,jstop,...))
    CartesianBox((istart:[istep:]istop, jstart:[jstep:]jstop, ...))

where `A` is an array (to define a region consisting in all the indices in the
array), `istart`, `istop`, etc. are integers (to define a region from
`(istart,jstart,...)` to `(istop,jstop,...)`.  Starting with Julia 1.6, a
non-unit step `istep`, `jstep`, etc. may also be specified.

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
[`intersect`](@ref).

"""
struct CartesianBox{N,I<:NTuple{N,IndexRange{Int}}} <: AbstractArray{CartesianIndex{N},N}
    indices::I
end

"""
    indices(A) -> inds

yields an `N`-tuple of index ranges if `A` is an instance of
[`CartesianBox{N}`](@ref) or of `CartesianIndices{N}`; yields an `N`-tuple of
indices if `A` is an instance of `CartesianIndex{N}` or an `N`-tuple of
integers.

"""
indices(R::CartesianBox) = getfield(R, :indices)
indices(R::CartesianIndices) = getfield(R, :indices)
indices(I::CartesianIndex) = Tuple(I)
indices(I::NTuple{N,Integer}) where {N} = I

Base.Tuple(R::CartesianBox) = indices(R)

# Fast conversion between CartesianIndices and CartesianBox.
Base.CartesianIndices(R::CartesianBox{N,I}) where {N,I} =
    CartesianIndices{N,I}(indices(R))
CartesianBox(R::CartesianIndices{N,I}) where {N,I} =
    CartesianBox{N,I}(indices(R))

# Other constructors.
CartesianBox(B::CartesianBox) = B
CartesianBox(A::AbstractArray) = CartesianBox(axes(A))
CartesianBox(inds::Tuple{Vararg{Union{<:Integer,IndexRange{<:Integer}}}}) =
        CartesianBox(CartesianIndices(inds))
CartesianBox(first::CartesianIndex{N}, last::CartesianIndex{N}) where {N} =
    CartesianBox(indices(first), indices(last))
CartesianBox(first::NTuple{N,Integer}, last::NTuple{N,Integer}) where {N} =
    CartesianBox(map((i,j) -> i:j, first, last))

@deprecate ranges(B::CartesianBox) indices(B) false

Base.first(B::CartesianBox) = first(CartesianIndices(B))
Base.last(B::CartesianBox) = last(CartesianIndices(B))
Base.size(B::CartesianBox) = size(CartesianIndices(B))
Base.axes(B::CartesianBox) = axes(CartesianIndices(B))

isempty(B::CartesianBox) = isempty_(first(B), last(B))
@inline isempty_(first::CartesianIndex{N}, last::CartesianIndex{N}) where {N} =
    any(map(isless, indices(last), indices(first)))

Base.show(io::IO, ::MIME"text/plain", B::CartesianBox) = show(io, B)
Base.show(io::IO, B::CartesianBox) = begin
    print(io, "CartesianBox(")
    print(io, indices(B))
    print(io, ")")
end

Base.view(A::AbstractArray{<:Any,N}, B::CartesianBox{N}) where {N} =
    view(A, indices(B)...)

function Base.getindex(A::AbstractArray{T,N}, B::CartesianBox{N}) where {T,N}
    empty = isempty(B)
    empty || isnonemptypartof(B, A) ||
        throw(BoundsError("sub-region is not part of array"))
    C = Array{T}(undef, size(B))
    if ! empty
        Ifirst = first(B)
        if any(i -> i != 1, indices(Ifirst))
            off = CartesianIndex(map(i -> i - 1, indices(Ifirst)))
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

Base.setindex!(A::AbstractArray{<:Any,N}, x, B::CartesianBox{N}) where {N} =
    A[indices(B)...] = x

Base.fill!(A::AbstractArray{T,N}, B::CartesianBox{N}, x) where {T,N} =
    fill!(A, B, convert(T, x)::T)
function Base.fill!(A::AbstractArray{T,N},
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

Base.IndexStyle(::Type{<:CartesianBox{N,R}}) where {N,R} = IndexStyle(R)
@inline @propagate_inbounds Base.getindex(B::CartesianBox, I...) =
    getindex(CartesianIndices(B), I...)

"""
    CartesianBoxable{N}

is a union of types (other than `CartesianBox{N}`) which can be automatically
converted into a `CartesianBox{N}`.

!!! note
    Although the constructor `CartesianBox` can be also applied to any instance
    of `AbstractArray`, this abstract type does not belong to the union
    `CartesianBoxable` as it is considered that such a conversion cannot be
    automatic.

"""
const CartesianBoxable{N} = Union{CartesianIndices{N},
                                  NTuple{N,IndexRange{<:Integer}},
                                  NTuple{N,Integer}}

# Extend eachindex() method.
Base.eachindex(::IndexCartesian, B::CartesianBox) = B

# Make CartesianBox iterable.
Base.iterate(iter::CartesianBox) = iterate(CartesianIndices(iter))
Base.iterate(iter::CartesianBox, state) = iterate(CartesianIndices(iter), state)

# Extend methods for fast SIMD iterations.
simd_outer_range(iter::CartesianBox{0}) = iter
simd_outer_range(iter::CartesianBox) = CartesianBox(tail(indices(iter)))
simd_inner_length(iter::CartesianBox{0}, ::CartesianIndex) = 1
simd_inner_length(iter::CartesianBox, I::CartesianIndex) = length(indices(iter)[1])
simd_index(iter::CartesianBox{0}, ::CartesianIndex, I1::Int) = first(iter)
@inline @propagate_inbounds function simd_index(iter::CartesianBox,
                                                Ilast::CartesianIndex, I1::Int)
    # For maximum portability, delegate work to do to the embedded
    # CartesianIndices.
    simd_index(CartesianIndices(iter), Ilast, I1)
end

# Increment a range.
incr(rng::AbstractUnitRange, adj::Number) =
    (first(rng) + adj):(last(rng) + adj)
incr(rng::AbstractRange, adj::Number) =
    (first(rng) + adj):step(rng):(last(rng) + adj)

# Decrement a range.
decr(rng::AbstractUnitRange, adj::Number) =
    (first(rng) - adj):(last(rng) - adj)
decr(rng::AbstractRange, adj::Number) =
    (first(rng) - adj):step(rng):(last(rng) - adj)

# Shifting of a CartesianBox.
function Base.:(+)(R::CartesianBox{N},
                   I::Union{NTuple{N,Integer},CartesianIndex{N}}) where {N}
    CartesianBox(map((r,i) -> incr(r, i), indices(R), indices(I)))
end
function Base.:(-)(R::CartesianBox{N},
                   I::Union{NTuple{N,Integer},CartesianIndex{N}}) where {N}
    CartesianBox(map((r,i) -> decr(r, i), indices(R), indices(I)))
end

"""
    intersect(CartesianBox, A, B)

yields the Cartesian box given by the intersection of the Cartesian regions
defined by `A` and `B`.  In this context, a Cartesian region can specified by a
Cartesian box, a list of integer valued ranges, a list of dimensions, or an
instance of `CartesianIndices`.

This method is equivalent to `intersect(A,B)`, or `A ∩ B` for short, when at
least one of `A` or `B` is a Cartesian box.

See also: [`CartesianBox`](@ref), [`isnonemptypartof`](@ref).

""" intersect

intersect(::Type{CartesianBox}, A::CartesianBox{N}, B::CartesianBox{N}) where {N} = A ∩ B
intersect(::Type{CartesianBox}, A::CartesianBox{N}, B::CartesianBoxable{N}) where {N} =
    A ∩ CartesianBox(B)
intersect(::Type{CartesianBox}, A::CartesianBoxable{N}, B::CartesianBox{N}) where {N} =
    CartesianBox(A) ∩ B
intersect(::Type{CartesianBox}, A::CartesianBoxable{N}, B::CartesianBoxable{N}) where {N} =
    CartesianBox(A) ∩ CartesianBox(B)

@deprecate intersection(A, B) intersect(CartesianBox, A, B) false

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

See also: [`CartesianBox`](@ref), [`intersect`](@ref).

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

See also: [`CartesianBox`](@ref), [`intersect`](@ref), [`isnonzero`](@ref).

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

# Provide equivalent methods for `typemin`, `typemax`, and `one` applied
# to a Cartesian index with different names to avoid type piracy.
@inline typemin_(::Type{<:CartesianIndex{N}}) where {N} =
    CartesianIndex(ntuple(i -> typemin(Int), Val(N)))
@inline typemax_(::Type{<:CartesianIndex{N}}) where {N} =
    CartesianIndex(ntuple(i -> typemax(Int), Val(N)))

one_(I::CartesianIndex) = one_(typeof(I))
@static if VERSION < v"1.1.0-rc1"
    one_(T::Type{<:CartesianIndex}) = one(T)
else
    one_(T::Type{<:CartesianIndex}) = oneunit(T)
end

"""
    isnonzero(x)

yields whether `x` is non-zero, the negation of `iszero(x)`.

"""
isnonzero(x) = !iszero(x)

end # module
