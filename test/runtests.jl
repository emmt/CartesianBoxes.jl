
if ! isdefined(Base, :CartesianBoxes)
    include(joinpath("..", "src", "CartesianBoxes.jl"))
end

module CartesianBoxesTests

using Compat
using Compat.Test
using CartesianBoxes

# Deal with compatibility issues.
@static if isdefined(Base, :CartesianIndices)
    import Base: CartesianIndices
    const CARTESIAN_REGIONS = (CartesianIndices, CartesianBox)
    const CHECK_CARTESIAN_RANGES = false
else
    import Base: CartesianRange
    import Compat: CartesianIndices
    const CARTESIAN_REGIONS = (CartesianIndices, CartesianBox, CartesianRange)
    const CHECK_CARTESIAN_RANGES = true
end
const CartesianRegions = Union{CARTESIAN_REGIONS...}
@static if isdefined(Base, :axes)
    import Base: axes
else
    import Base: indices
    const axes = indices
end

@inline _add(a,b) = a + b
@inline _sub(a,b) = a - b
@inline _mul(a,b) = a * b
@inline _div(a,b) = a / b

for op in (:_add, :_sub, :_mul, :_div)
    slow = Symbol(:slow,op,:!)
    fast = Symbol(:fast,op,:!)
    @eval begin
        function $slow(dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            for i in eachindex(dst, A, B)
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end

        function $slow(::Type{R},
                       dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {R<:CartesianRegions,T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            for i in R(inds)
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end

        function $slow(R::CartesianRegions,
                       dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            @assert isnonemptypartof(CartesianBox(R), A)
            for i in R
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end

        # Implementation using what is provided by Julia base (probably linear
        # indexing for regular arrays).
        function $fast(dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            @inbounds @simd for i in eachindex(dst, A, B)
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end

        # Implementation using a specific iterator.
        function $fast(::Type{R},
                       dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {R<:CartesianRegions,T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            @inbounds @simd for i in R(inds)
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end

        # Operation over a sub-region.
        function $fast(R::CartesianRegions,
                       dst::AbstractArray{T,N},
                       A::AbstractArray{T,N},
                       B::AbstractArray{T,N}) where {T,N}
            inds = axes(dst)
            @assert axes(A) == inds
            @assert axes(B) == inds
            @assert isnonemptypartof(CartesianBox(R), A)
            @inbounds @simd for i in R
                dst[i] = $op(A[i], B[i])
            end
            return dst
        end
    end
end

maxabsdif(A, B) = maximum(abs.(A .- B))

function stupidcount(iter)
    n = 0
    for i in iter
        n += 1
    end
    return n
end

SIZES = (345, (21,22), (11,12,13), (5,6,7,8))
TYPES = (Float64, Float32)

@testset "CartesianBoxes" begin
    @testset "Basic operations" for dims in SIZES, T in TYPES
        A = Array{T}(undef, dims...)
        r = CartesianIndices(A)
        b = CartesianBox(A)
        @test ndims(b) == ndims(r)
        @test eltype(b) == eltype(r)
        @test length(b) == stupidcount(b)
        @test length(b) == length(r)
        @test size(b) == size(r)
        @test axes(b) == axes(r)
        for d in 1:ndims(b)
            @test size(b,d) == size(r,d)
            @test axes(b,d) == axes(r,d)
        end

        @test CartesianBox(r) == b
        @test CartesianIndices(b) == r
        @static if CHECK_CARTESIAN_RANGES
            s = CartesianRange(axes(A))
            @test CartesianBox(s) == b
            @test CartesianRange(b) == s
        end
    end

    @testset "Array indexing" for dims in SIZES, T in TYPES
        A = randn(T, dims...)
        B = randn(T, dims...)
        C1 = Array{T}(undef, dims...)
        C2 = Array{T}(undef, dims...)
        C3 = Array{T}(undef, dims...)

        @testset "Operation $oper" for oper in (slow_add!, slow_sub!,
                                                slow_mul!)
            oper(C1, A, B)
            oper(CartesianBox, fill!(C2, 0), A, B)
            @test maxabsdif(C1, C2) == 0
            if ndims(A) ≥ 1
                # Take a sub-region.
                S = CartesianBox(map(n -> (min(2,n) : max(n-1,1)), size(A)))
                oper(CartesianBox(S),     fill!(C2, 0), A, B)
                oper(CartesianIndices(S), fill!(C3, 0), A, B)
                @test maxabsdif(C2, C3) == 0
            end
        end

        @testset "Operation $oper" for oper in (fast_add!, fast_sub!,
                                                fast_mul!)
            oper(C1, A, B)
            oper(CartesianBox, fill!(C2, 0), A, B)
            @test maxabsdif(C1, C2) ≤ 4*eps(T)
            if ndims(A) ≥ 1
                # Take a sub-region.
                S = CartesianBox(map(n -> (min(2,n) : max(n-1,1)), size(A)))
                oper(CartesianBox(S),     fill!(C2, 0), A, B)
                oper(CartesianIndices(S), fill!(C3, 0), A, B)
                @test maxabsdif(C2, C3) == 0
            end
         end
    end
end

end # module
