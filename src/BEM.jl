module BEM

using Distributed

include(joinpath(@__DIR__, "Fault.jl"))
include(joinpath(@__DIR__, "RateStateFriction.jl"))
include(joinpath(@__DIR__, "dc3d.jl"))

using SimpleTraits
using Parameters
using ToeplitzMatrices
using TensorOperations
using FFTW

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
abstract type BoundaryElementGrid{dim, isuniform} <: AbstractDifferenceGrids{dim, isuniform} end

struct BEMGrid_1D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{1, true}
    ξ::U1 # along-downdip
    Δξ::T
    nξ::I
    aξ::U2 # cell boundary along-downdip (width)
end

"Along-strike will be places on x-axis while along-downdip on yz-plane, w.r.t Okada's dc3d coordinates."
struct BEMGrid_2D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{2, true}
    x::U1 # along-strike
    ξ::U1 # along-downdip
    Δx::T
    Δξ::T
    nx::I
    nξ::I
    ax::U2 # cell boundary along-strike (length)
    aξ::U2 # cell boundary along-downdip (width)
end

function discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ) where {ftype}
    ξ, nξ, aξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    BEMGrid_1D(ξ, Δξ, nξ, aξ)
end

function discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ) where {ftype}
    ξ, nξ, aξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    x, nx, ax = _divide_segment(Val(:symmatzero), fa[:x], Δx)
    BEMGrid_2D(x, ξ, Δx, Δξ, nx, nξ, ax, aξ)
end

function _divide_segment(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(zero(T), stop=-x+Δx, step=-Δx)) .- Δx/2
    xs = cat(xi .- Δx/2, xi .+ Δx/2, dims=2)
    return xi, length(xi), xs
end

function _divide_segment(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx))
    xs = cat(xi .- Δx/2, xi .+ Δx/2, dims=2)
    return xi, length(xi), xs
end

abstract type AbstractElasticProperties end

"""
Okada's [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) only applies on isotropic materials,
    therefore, elastic modulus are constrained to be scalars except for `ρ` and `cs`.
"""
@with_kw struct HomogeneousElasticProperties{T<:Number} <: AbstractElasticProperties
    λ::T # Lamé's first parameter
    μ::T # shear modulus
    cs::T # shear wave velocity
    α::T = (λ + μ) / (λ + 2μ) # material constants which equals (λ + μ) / (λ + 2μ)
end

@with_kw struct StaticPlaneMaterialProperties{dim, U<:AbstractVecOrMat, P<:AbstractArray, T<:Number} <: AbstractMaterialProperties{dim}
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    k::P # stiffness tensor
    σ::U # effective normal stress
    η::U # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    function StaticPlaneMaterialProperties(a::U, b::U, L::U, k::P, σ::U, η::U, vpl::T, f0::T, v0::T) where {U, P, T}
        dims = maximum([ndims(x) for x in (a, b, L, σ, η)])
        new{dims, U, P, T}(a, b, L, k, σ, η, vpl, f0, v0)
    end
end

"""
    shear_traction(::Type{<:PlaneFault}, u12, λ, μ, dip)

Calculate the shear traction on the fault plane w.r.t. fault types.

# Arguments
- `u::AbstractArray{<:Number, 1}`: the output from dc3d_okada
- `λ::Number`: Lamé's first parameter
- `μ::Number`: shear modulus
- `dip::Number`: plane dip angle

# Reference
- A good reference is at [Displacement & Strain & Stress](https://nnakata.oucreate.com/page/Teaching_files/GEOPHYS130/GEOPHYS130_notes_all.pdf).

"""
shear_traction(ftype::Type{<:AbstractFault}) = NaN

@inline function shear_traction(::Type{FT}, u, λ, μ, dip) where {FT<:PlaneFault}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline function shear_traction(::Type{StrikeSlipFault}, u, λ, μ, dip)
    # Special case for `dip = 90.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

"""
Periodic boundary condition for 2D faults.

# Arguments
- same as `dc3d_okada`, see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
"""
function stiffness_periodic_boundary_condition(x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T}, lf::T, nrept::Integer) where {T}
    u = zeros(T, 12)
    @fastmath for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i*lf, aw, disl)
    end
    return u
end

"""
    stiffness_tensor(fa::PlaneFaultDomain, gd::BEMGrid, mp::BEM, ::Vararg)

# Arguments
- lengthexpensionratio: ratio of along-strike length against along downdip for mimicing infinitely extending 1D fault

Calculate the reduced stiffness tensor based on the symmetric properties for 2D case.
Results for 1D fault will be 2D matrix and 3D Array for 2D faults.
"""
function stiffness_tensor(fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1, true}, ep::HomogeneousElasticProperties, ax_ratio=50; ncpus=1) where {ftype, T}
    x, y, z, dep = zero(T), gd.ξ .* cospi(fa.dip / 180), gd.ξ .* sinpi(fa.dip / 180), zero(T)
    disl = applied_unit_dislocation(ftype)
    ax = fa[:ξ] * ax_ratio .* [-one(T), one(T)]
    u12 = zeros(T, 12)

    if ncpus == 1
        k = zeros(T, gd.nξ, gd.nξ)
        @inbounds @threads for j = 1: gd.nξ
            for i = 1: gd.nξ
                u12 .= dc3d_okada(x, y[i], z[i], ep.α, dep, fa.dip, ax, gd.aξ[j,:], disl)
                k[i,j] = shear_traction(ftype, u12, ep.λ, ep.μ, fa.dip)
            end
        end
    else
        k = SharedArray{T}(gd.nξ, gd.nξ)
        @sync @distributed for j = 1: gd.nξ
            for i = 1: gd.nξ
                u12 .= dc3d_okada(x, y[i], z[i], ep.α, dep, fa.dip, ax, gd.aξ[j,:], disl)
                k[i,j] = shear_traction(ftype, u12, ep.λ, ep.μ, fa.dip)
            end
        end
        k = sdata(k)
    end
    return k
end

function stiffness_tensor(fa::PlaneFaultDomain{ftype, 2, T}, gd::BoundaryElementGrid{2, true}, ep::HomogeneousElasticProperties) where {ftype, T}

end

"""
Return the proper style of stiffness tensor for 2D cases in consideration of using proper contraction methods.

    * If k contains completely 4 dims, use `@tensor` from `TensorOperations`
    * If k is of toeplitz matrix, use BLAS implemented in `ToeplitzMatrices`
    * If k is reduced, i.e. 3 dims, use fft from `FFTW`

The latter two utilize the symmetry properties in stiffness tensor of 2D cases by converting contraction to circular convolution.
"""
function full_stiffness_tensor(contraction_method, k)

end

"""
For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault
    if dip angle are constrained within [0, π/2] in fact.

The unit of *unit dislocation* below is the same of `v * t` at set by user so to avoid normalization step.
"""
applied_unit_dislocation(::Type{NormalFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{ThrustFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{StrikeSlipFault}) = [1., 0., 0.]

"Temporal variable in solving ODEs aimed to avoid alloction overheads."
struct tmp_variable{T<:AbstractVecOrMat}
    dμ_dt::T
    dμ_dθ::T
    dμ_dv::T
    ψ1::T
    ψ2::T
end

end # module
