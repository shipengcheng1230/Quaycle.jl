module BEM

include(joinpath(@__DIR__, "Fault.jl"))
include(joinpath(@__DIR__, "RateStateFriction.jl"))
include(joinpath(@__DIR__, "dc3d.jl"))

using SimpleTraits
using Parameters
using ToeplitzMatrices
using TensorOperations

using Base.Threads

using .RateStateFriction
using .Fault

export AbstractDifferenceGrids, BEMGrid

"""
# Traits:
Symmetric properties based on `BEMGrid` on 2D fault plane.

- Translational Symmetry:
    Stiffness remain the same if both observation and dislocation patches are moved by the same displacement.
- Reflective (A)Symmetry:
    Stiffness remain the same(opposite) if swapping the along strike position between observation and dislocation patches.
"""
@traitdef IsTranslationalSymmetrical{ftype}
@traitdef IsReflectiveSymmetrical{ftype}
@traitdef IsReflectiveASymmetrical{ftype}

@traitimpl IsTranslationalSymmetrical{NormalFault}
@traitimpl IsTranslationalSymmetrical{ThrustFault}
@traitimpl IsTranslationalSymmetrical{StrikeSlipFault}

@traitimpl IsReflectiveSymmetrical{NormalFault}
@traitimpl IsReflectiveSymmetrical{ThrustFault}
@traitimpl IsReflectiveASymmetrical{StrikeSlipFault}

abstract type AbstractDifferenceGrids{dim, isuniform} end

"Along-strike will be places on x-axis while along-downdip on yz-plane, w.r.t Okada's dc3d coordinates."
struct BEMGrid{N, U<:AbstractVector} <: AbstractDifferenceGrids{N, true}
    x::U # along-strike
    ξ::U # along-downdip
    Δx
    Δξ
    nx
    nξ

    function BEMGrid(x::U, ξ::U, Δx, Δξ, nx, nξ) where {U<:AbstractVector}
        N = isempty(x) ? 1 : 2
        N == 2 ? new{N, U}(x, ξ, Δx, Δξ, nx, nξ) : new{N, U}(x, ξ, nothing, Δξ, nothing, nξ)
    end
end

BEMGrid(ξ, Δξ, nξ) = BEMGrid([], ξ, nothing, Δξ, nothing, nξ)

function discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ) where {ftype}
    ξ, nξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    BEMGrid(ξ, Δξ, nξ)
end

function discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ) where {ftype}
    ξ, nξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    x, nx = _divide_segment(Val(:symmatzero), fa[:x], Δx)
    BEMGrid(x, ξ, Δx, Δξ, nx, nξ)
end

function _divide_segment(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x+Δx, stop=zero(T), step=Δx)) .- Δx/2
    return xi, length(xi)
end

function _divide_segment(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx))
    return xi, length(xi)
end

abstract type AbstractElasticProperties end
"""
Okada's [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) only applies on isotropic materials,
therefore, elastic modulus are constrained to be scalars except for `ρ` and `cs`.
"""
@with_kw struct IsotropicElasticProperties{T} <: AbstractElasticProperties
    λ::T # Lamé's first parameter
    μ::T # shear modulus
    α::T = (λ + μ) / (λ + 2μ) # material constants which equals (λ + μ) / (λ + 2μ)
    ρ = μ ./ cs^2 # density
    cs = sqrt(μ ./ ρ) # shear wave velocity
end

@with_kw struct StaticMaterialProperties{dim, perturb, tensorMD, T} <: AbstractMaterialProperties{dim}
    a::AbstractVecOrMat{T} # contrib from velocity
    b::AbstractVecOrMat{T} # contrib from state
    L::AbstractVecOrMat{T} # critical distance
    k::AbstractArray{T} # stiffness tensor
    σ::AbstractVecOrMat{T} # effective normal stress
    η::AbstractVecOrMat{T} # radiation damping
    vpl::AbstractVecOrMat{T} # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T # ref. frictional coeff
    v0::T # ref. velocity

    function MaterialProperties(a, b, L, k, σ, η, vpl, f0, v0)
        dims = maximum([ndims(x) for x in (a, b, L, σ, η, vpl)])
        # If k contains completely 4 dims, use `@tensor` from TensorOperations
        # If k is of toeplitz matrix, use BLAS implemented in ToeplitzMatrices
        # If k is reduced, i.e. 3 dims, use fft from FFTW
        ndims(k) == 4 ? tensorMD = Val(:tensor) : eltype(k) <: ToeplitzMatrices.AbstractToeplitz ? tensorMD = Val(:toeplitz) : tensorMD = Val{:fft}
        new{dims, Nothing}(a, b, L, k, σ, η, vpl, f0, v0)
    end
end

"""
    shear_traction(::Type{<:PlaneFault}, u12, λ, μ, dip)

Calculate the shear traction on the fault plane w.r.t. fault types.

# Arguments
- `u12::AbstractArray{<:Number, 1}`: the output from dc3d_okada
- `λ::Number`: Lamé's first parameter
- `μ::Number`: shear modulus
- `dip::Number`: plane dip angle

# Reference
- A good reference is at [Displacement & Strain & Stress](https://nnakata.oucreate.com/page/Teaching_files/GEOPHYS130/GEOPHYS130_notes_all.pdf).

"""
shear_traction(ftype::Type{<:AbstractFault}) = NaN

@traitfn @inline function shear_traction(::Type{FT}, u12, λ, μ, dip) where {FT<:PlaneFault; IsReflectiveSymmetrical{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@traitfn @inline function shear_traction(::Type{FT}, u12, λ, μ, dip) where {FT<:PlaneFault; IsReflectiveASymmetrical{FT}}
    # Special case for `dip = 0.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

"""
Calculate relationship between shear traction at one location and displacement at another.

# Arguments
- same as `dc3d_okada`, see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
"""
function stiffness_element(x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T}, nrept::Integer) where {T}
    u = zeros(T, 12)
    @fastmath for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i*lf, aw, disl)
    end
    return u
end

"""
    stiffness_tensor(fa::PlaneFaultDomain, gd::BEMGrid, mp::BEM)

Calculate the reduced stiffness tensor based on the symmetric properties.
"""
function stiffness_tensor(fa::PlaneFaultDomain{1}, gd::BEMGrid{1}, mp::IsotropicElasticProperties)
    cdip = cospi(fa.dip / 180)
    sdip = sinpi(fa.dip / 180)

end

end # module
