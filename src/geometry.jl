## geometry of simulation domain and naive grid generator

export gen_mesh, SimpleLine, SimpleLineGrid, SimpleRect, SimpleRectGrid

abstract type AbstractGeometry{dim} end
abstract type SimpleGeometry{dim} <: AbstractGeometry{dim} end

abstract type AbstractMesh{dim} end
abstract type SimpleMesh{dim} <: AbstractMesh{dim} end

struct SimpleLine{A<:Number, T<:NTuple{3,A}, U<:NTuple{1}} <: SimpleGeometry{1}
    origin::T
    span::U
end

struct SimpleRect{A<:Number, T<:NTuple{3,A}, U<:NTuple{2}} <: SimpleGeometry{2}
    origin::T
    span::U
end

struct SimpleLineGrid{T<:AbstractArray, U<:Number, I<:Integer, S<:AbstractArray} <: SimpleMesh{1}
    ξ::T
    Δξ::U
    nξ::I
    aξ::S
end

struct SimpleRectGrid{T<:AbstractArray, U<:Number, I<:Integer, S<:AbstractArray} <: SimpleMesh{2}
    x::T
    ξ::T
    Δx::U
    Δξ::U
    nx::I
    nξ::I
    ax::S
    aξ::S
end

function gen_mesh(::Val{:SimpleGrid}, geo::SimpleLine, dξ::Number)
    ξ, nξ, aξ = __divide_segment__(Val(:halfspace), promote(geo.span[1], dξ)...)
    return SimpleLineGrid(ξ, convert(eltype(ξ), dξ), nξ, aξ)
end

function gen_mesh(::Val{:SimpleGrid}, geo::SimpleRect, dx::Number, dξ::Number)
    ξ, nξ, aξ = __divide_segment__(Val(:halfspace), promote(geo.span[2], dξ)...)
    x, nx, ax = __divide_segment__(Val(:symmatzero), promote(geo.span[1], dx)...)
    return SimpleRectGrid(x, ξ, convert(eltype(x), dx), convert(eltype(ξ), dξ), nx, nξ, ax, aξ)
end

function __divide_segment__(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = range(zero(T), stop=-x+Δx, step=-Δx) .- Δx/2 |> collect
    xs = __segment_boarder__(xi, Δx)
    return xi, length(xi), xs
end

function __divide_segment__(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx) |> collect
    xs = __segment_boarder__(xi, Δx)
    return xi, length(xi), xs
end

__segment_boarder__(xi::AbstractVector{<:Number}, Δx::Number) = [[xx - Δx/2, xx + Δx/2] for xx in xi]
