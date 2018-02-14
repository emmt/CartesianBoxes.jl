#
# BoundingBoxes.jl -
#
# Extends CartesianRange.
#
#-------------------------------------------------------------------------------
#
# Copyright (c) 2017 Éric Thiébaut.
#
# All rights reserved.
#

module BoundingBoxes

export
    BoundingBox,
    boundingbox,
    intersection,
    isnonemptypartof

"""

`BoundingBox{N}` and `BoundingBox(args...)` are shortcuts for
`CartesianRange{CartesianIndex{N}}` and `CartesianRange(args...)`.  In addition
to the usual ways to build a Cartesian range, the following methods are
provided:

    BoundingBox(A) -> BoundingBox(size(A))
    BoundingBox((imin,jmin,...), (imax,jmax,...))
    BoundingBox(imin:imax, jmin:jmax, ...)

where `A` is an array and `imin`, `imax`, etc. are integers.  These definitions
do not interfere with those of `CartesianRange`.

See also: `boundingbox`[@ref], `CartesianRange`[@ref], `CartesianIndex`[@ref],
 `intersection`[@ref].

"""
const BoundingBox{N} = CartesianRange{CartesianIndex{N}}
BoundingBox(args...) = CartesianRange(args...)
BoundingBox(A::AbstractArray) = CartesianRange(size(A))
BoundingBox(start::NTuple{N,Integer}, stop::NTuple{N,Integer}) where {N} =
    CartesianRange(CartesianIndex(start), CartesianIndex(stop))
BoundingBox(rngs::AbstractUnitRange...) =
    CartesianRange(rngs)

"""

    intersection(R, S)

yields the Cartesian range which is the intersection of the two Cartesian
ranges `R` and `S`.  This method is similar to `intersect(R,S) = R ∩ S` which
yields an array of Cartesian indices and is **much** slower (and less useful).

See also: `BoundingBox`[@ref], `isnonemptypartof`[@ref].

"""
intersection(R::BoundingBox{N}, S::BoundingBox{N}) where {N} =
    BoundingBox(max(first(R), first(S)), min(last(R), last(S)))

"""
    isnonemptypartof(R, S)

yields whether region `R` is nonempty and a valid part of `S`.  If `R` is a
Cartesian range, then `S` can be an array to check whether `R` is a valid
nonempty region of interest of `S`.

See also: `BoundingBox`[@ref], `intersection`[@ref].

"""
isnonemptypartof(R, S) = false

function isnonemptypartof(R::BoundingBox{N},
                          A::AbstractArray{T,N}) where {T,N}
    isnonemptypartof(R, CartesianRange(size(A)))
end

function isnonemptypartof(R::BoundingBox{N},
                          S::BoundingBox{N}) where {N}
    first(S) ≤ first(R) ≤ last(R) ≤ last(S)
end


"""
    boundingbox(A)

yields the bounding-box of non-zero values in array `A`.  To only consider a
sub-region `B` of `A`, call:

    boundingbox(A, B)

FIXME: The algorithm is pretty silly for now and could be faster than
       `O(length(A))`.

See also: `BoundingBox`[@ref], `intersection`[@ref].

"""
boundingbox(A::AbstractArray) = _boundingbox(A, BoundingBox(A))

boundingbox(A::AbstractArray{T,N}, B::BoundingBox{N}) where {T,N} =
    _boundingbox(A, intersection(BoundingBox(A), B))

function _boundingbox(A::AbstractArray{T,N},
                      B::BoundingBox{N}) where {T,N}
    Imin = CartesianIndex(ntuple(i -> typemax(Int), Val{N}))
    Imax = CartesianIndex(ntuple(i -> typemin(Int), Val{N}))
    for I in B
        if A[i] != zero(T)
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return BoundingBox(Imin, Imax)
end

function _boundingbox(A::AbstractArray{Bool,N},
                      B::BoundingBox{N}) where {N}
    Imin = CartesianIndex(ntuple(i -> typemax(Int), Val{N}))
    Imax = CartesianIndex(ntuple(i -> typemin(Int), Val{N}))
    for I in B
        if A[i]
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return BoundingBox(Imin, Imax)
end

end # module
