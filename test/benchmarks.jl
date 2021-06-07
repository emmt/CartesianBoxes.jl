#isdefined(Main, :CartesianBoxes) ||
#    include(joinpath("..", "src", "CartesianBoxes.jl"))
#push!(LOAD_PATH, normpath(joinpath(dirname(@__FILE__), "../src")))

module CartesianBoxesBenchmarks

using CartesianBoxes
using BenchmarkTools

# Deal with compatibility issues.
using Compat
@static if isdefined(Base, :CartesianIndices)
    import Base: CartesianIndices
    const CARTESIAN_REGIONS = (CartesianIndices, CartesianBox)
else
    import Base: CartesianRange
    import Compat: CartesianIndices
    const CARTESIAN_REGIONS = (CartesianIndices, CartesianRange, CartesianBox)
end
const CartesianRegions = Union{CARTESIAN_REGIONS...}
@static if isdefined(Base, :axes)
    import Base: axes
else
    import Base: indices
    const axes = indices
end

function add!(dst::AbstractArray{T,N},
              A::AbstractArray{T,N},
              B::AbstractArray{T,N}) where {T,N}
    inds = axes(dst)
    @assert axes(A) == inds
    @assert axes(B) == inds
    @inbounds @simd for i in eachindex(dst, A, B)
        dst[i] = A[i] + B[i]
    end
    return dst
end

function add!(::typeof(eachindex),
              dst::AbstractArray{T,N},
              A::AbstractArray{T,N},
              B::AbstractArray{T,N}) where {T,N}
    inds = axes(dst)
    @assert axes(A) == inds
    @assert axes(B) == inds
    @inbounds @simd for i in eachindex(dst, A, B)
        dst[i] = A[i] + B[i]
    end
    return dst
end

function add!(::Type{R},
              dst::AbstractArray{T,N},
              A::AbstractArray{T,N},
              B::AbstractArray{T,N}) where {R<:CartesianRegions,T,N}
    inds = axes(dst)
    @assert axes(A) == inds
    @assert axes(B) == inds
    @inbounds @simd for i in R(inds)
        dst[i] = A[i] + B[i]
    end
    return dst
end

function inner(A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    @assert axes(A) == axes(B)
    s = zero(T)
    @inbounds @simd for i in eachindex(A, B)
        s += A[i]*B[i]
    end
    return s
end

function inner(::typeof(eachindex),
               A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    @assert axes(A) == axes(B)
    s = zero(T)
    @inbounds @simd for i in eachindex(A, B)
        s += A[i]*B[i]
    end
    return s
end

function inner(::Type{R},
               A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {R<:CartesianRegions,
                                             T<:AbstractFloat,N}
    inds = axes(A)
    @assert axes(B) == inds
    s = zero(T)
    @inbounds @simd for i in R(inds)
        s += A[i]*B[i]
    end
    return s
end

T = Float32
dims = (20,30,40)
A = rand(T,dims...)
B = rand(T,dims...)
C1 = Array{T}(undef, dims)
C2 = Array{T}(undef, dims)
C3 = Array{T}(undef, dims)

brief(::typeof(eachindex)) = "eachindex";
brief(::Type{CartesianBox}) = "CartesianBox";
@static if isdefined(Base, :CartesianRange)
    brief(::Type{CartesianRange}) = "CartesianRange";
end
@static if isdefined(Base, :CartesianIndices)
    brief(::Type{CartesianIndices}) = "CartesianIndices";
end

const VERBOSE = false
add!(C1, A, B)

println("Addition of arrays of size $dims ($(prod(dims)) operations)")
for R in (eachindex, CARTESIAN_REGIONS...)
    flg = (R == last(CARTESIAN_REGIONS))
    pfx1 = (flg ? " └─ " : " ├─ ")
    pfx2 = (flg ? "    " : " │  ")
    println(pfx1, "Computations with `$(brief(R))`:")

    fill!(C2, 0)
    add!(R, C2, A, B)
    println(pfx2, " ├─ diff. ∈ ", extrema(C1 - C2))

    fill!(C3, 0)
    if VERBOSE
        t = @benchmark add!($R, $C3, $A, $B)
        show(stdout, MIME"text/plain"(), t)
    else
        print(pfx2, " └─ @btime:")
        @btime add!($R, $C3, $A, $B)
    end
end

w = inner(A, B)
println("Inner of arrays of size $dims ($(2*prod(dims)) operations)")
for R in (eachindex, CARTESIAN_REGIONS...)
    flg = (R == last(CARTESIAN_REGIONS))
    pfx1 = (flg ? " └─ " : " ├─ ")
    pfx2 = (flg ? "    " : " │  ")
    println(pfx1, "Computations with  `$(brief(R))`:")

    wp = inner(R, A, B)
    println(pfx2, " ├─ result = $wp ($w)")

    if VERBOSE
        t = @benchmark inner($R, $A, $B)
        show(stdout, MIME"text/plain"(), t)
    else
        print(pfx2, " └─ @btime:")
        @btime inner($R, $A, $B)
    end
end

#show(STDOUT, MIME"text/plain"(), @benchmark $p(dat, a, img, b))

end # module
