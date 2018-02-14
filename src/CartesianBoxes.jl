#
# CartesianBoxes.jl -
#
# Extends CartesianRange.
#
#-------------------------------------------------------------------------------
#
# This file is part of the `CartesianBoxes.jl` package which is licensed under
# the MIT "Expat" License.
#
# Copyright (c) 2017 Éric Thiébaut.
#

__precompile__(true)

module CartesianBoxes

export
    CartesianBox,
    boundingbox,
    intersection,
    isnonemptypartof

"""

`CartesianBox{N}` and `CartesianBox(args...)` are shortcuts for
`CartesianRange{CartesianIndex{N}}` and `CartesianRange(args...)`.  In addition
to the usual ways to build a Cartesian range, the following methods are
provided:

    CartesianBox(A) -> CartesianBox(indices(A))
    CartesianBox((imin,jmin,...), (imax,jmax,...))
    CartesianBox(imin:imax, jmin:jmax, ...)

where `A` is an array and `imin`, `imax`, etc. are integers.  These definitions
do not interfere with those of `CartesianRange`.

See also: `boundingbox`[@ref], `CartesianRange`[@ref], `CartesianIndex`[@ref],
 `intersection`[@ref].

"""
const CartesianBox{N} = CartesianRange{CartesianIndex{N}}
CartesianBox(args...) = CartesianRange(args...)
CartesianBox(A::AbstractArray) = CartesianRange(indices(A))
CartesianBox(start::NTuple{N,Integer}, stop::NTuple{N,Integer}) where {N} =
    CartesianRange(CartesianIndex(start), CartesianIndex(stop))
CartesianBox(rngs::AbstractUnitRange...) =
    CartesianRange(rngs)

"""

    intersection(R, S)

yields the Cartesian range which is the intersection of the two Cartesian
ranges `R` and `S`.  This method is similar to `intersect(R,S) = R ∩ S` which
yields an array of Cartesian indices and is **much** slower (and less useful).

See also: `CartesianBox`[@ref], `isnonemptypartof`[@ref].

"""
intersection(R::CartesianBox{N}, S::CartesianBox{N}) where {N} =
    CartesianBox(max(first(R), first(S)), min(last(R), last(S)))

"""
    isnonemptypartof(R, S)

yields whether region `R` is nonempty and a valid part of `S`.  If `R` is a
Cartesian range, then `S` can be an array to check whether `R` is a valid
nonempty region of interest of `S`.

See also: `CartesianBox`[@ref], `intersection`[@ref].

"""
isnonemptypartof(R, S) = false

function isnonemptypartof(R::CartesianBox{N},
                          A::AbstractArray{T,N}) where {T,N}
    isnonemptypartof(R, CartesianBox(A))
end

function isnonemptypartof(R::CartesianBox{N},
                          S::CartesianBox{N}) where {N}
    first(S) ≤ first(R) ≤ last(R) ≤ last(S)
end


"""
    boundingbox(A)

yields the bounding-box of non-zero values in array `A`.  To only consider a
sub-region `B` of `A`, call:

    boundingbox(A, B)

FIXME: The algorithm is pretty silly for now and could be faster than
       `O(length(A))`.

See also: `CartesianBox`[@ref], `intersection`[@ref].

"""
boundingbox(A::AbstractArray) = _boundingbox(A, CartesianBox(A))

boundingbox(A::AbstractArray{T,N}, B::CartesianBox{N}) where {T,N} =
    _boundingbox(A, intersection(CartesianBox(A), B))

function _boundingbox(A::AbstractArray{T,N},
                      B::CartesianBox{N}) where {T,N}
    Imin = CartesianIndex(ntuple(i -> typemax(Int), Val{N}))
    Imax = CartesianIndex(ntuple(i -> typemin(Int), Val{N}))
    for I in B
        if A[i] != zero(T)
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return CartesianBox(Imin, Imax)
end

function _boundingbox(A::AbstractArray{Bool,N},
                      B::CartesianBox{N}) where {N}
    Imin = CartesianIndex(ntuple(i -> typemax(Int), Val{N}))
    Imax = CartesianIndex(ntuple(i -> typemin(Int), Val{N}))
    for I in B
        if A[i]
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return CartesianBox(Imin, Imax)
end

end # module
