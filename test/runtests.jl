module TestingCartesianBoxes

using Test
using CartesianBoxes
using CartesianBoxes: indices, ranges

# Deal with compatibility issues.
import Base: CartesianIndices, axes
const CartesianRegions = Union{CartesianIndices, CartesianBox}

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

# This method is to exercise the loop and check the computed length.
function stupidcount(iter)
    n = 0
    for i in iter
        n += 1
    end
    return n
end

SIZES = ((), (45,), (21,22), (11,12,13), (5,6,7,8))
TYPES = (Float64, Float32)

@testset "CartesianBoxes" begin
    @testset "Basic operations" for dims in SIZES, T in TYPES
        A = Array{T}(undef, dims)
        r = CartesianIndices(A)
        b = CartesianBox(A)
        @test ndims(b) == ndims(r)
        @test eltype(b) === eltype(r)
        @test length(b) == stupidcount(b)
        @test length(b) == length(r)
        @test size(b) == size(r)
        @test axes(b) == axes(r)
        for d in 1:ndims(b)
            @test size(b,d) == size(r,d)
            @test axes(b,d) == axes(r,d)
        end
        @test first(b) === first(r)
        @test last(b) === last(r)
        buf = IOBuffer();
        show(buf, MIME"text/plain"(), b);
        @test String(take!(buf)) == repr(b)

        @test CartesianBox(r) == b
        @test CartesianIndices(b) == r
        @test CartesianIndices(CartesianBox(r)) === r
        @test CartesianBox(CartesianIndices(b)) === b

        I0 = CartesianIndex(ntuple(i -> 0, ndims(b)))
        I1 = CartesianIndex(ntuple(i -> i, ndims(b)))
        @test indices(I1) === I1.I
        @test indices(b) == axes(b)
        @test (b + I0) == b
        @test (b - I0) == b
        @test (b + I0.I) == b
        @test (b - I0.I) == b
        for f in (+, -)
            @test f(b,I1) == CartesianBox(map((r,i) -> f(first(r),i):f(last(r),i),
                                              axes(b), indices(I1)))
        end
    end

    @testset "Array indexing" for dims in SIZES, T in TYPES
        A = randn(T, dims)
        B = randn(T, dims)
        C1 = Array{T}(undef, dims)
        C2 = Array{T}(undef, dims)
        C3 = Array{T}(undef, dims)

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

    # FIXME: for now, we skip tests with ndims = 0
    @testset "Array functions" for dims in SIZES[2:end], T in (TYPES..., Bool)
        A = rand([zero(T), one(T)], dims)
        B = boundingbox(A)
        fill!(A, B, 0)
        C = boundingbox(A)
        @test isempty(C)
        @test length(C) == 0
        B = CartesianBox(map(r -> first(r)+1:last(r)-1, axes(A)))
        C = CartesianBox(map(r -> first(r)+2:last(r)-2, axes(A)))
        fill!(A, zero(T))
        fill!(A, B, one(T))
        @test boundingbox(A) === B
        A[B] .= zero(T)
        @test isempty(boundingbox(A))
        if T !== Bool
            X = rand(T, size(C))
            fill!(A, zero(T))
            fill!(A, B, typemax(T))
            A[C] = X
            @test all(map(isequal, A[C], X))
            @test boundingbox(x -> x != zero(T) && x < typemax(T), A) == C
            @test boundingbox(x -> x != typemax(T), A, B) == C
            inds = CartesianBoxes.indices(B)
            @test (@test_deprecated ranges(B)) === inds
            @test boundingbox(x -> x != typemax(T), A, inds) == C
            V = view(A,C)
            @test all(map(isequal, V, X))
            # Indices starting at 1.
            @. A = randn(T)
            inds = map(r -> first(r):last(r)-1, axes(A))
            B = CartesianBox(inds)
            @test maxabsdif(A[inds...], A[B]) == 0
        end
    end

    @testset "Miscellaneous" begin
        B = CartesianBox((1:3, -7:6))
        @test eachindex(IndexCartesian(), B) == B
        Z = CartesianBox(())
        @test stupidcount(Z) == 1
    end
end

end # module
