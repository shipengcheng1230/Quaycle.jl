include(joinpath(@__DIR__, "dc3d.jl"))

using SimpleTraits
using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra

export AbstractDifferenceGrids, BoundaryElementGrid
export HomogeneousElasticProperties, PlaneMaterialProperties
export shear_traction, stiffness_tensor
export derivations!, EarthquakeCycleProblem

@traitdef IsOnYZPlane{faulttype}
@traitimpl IsOnYZPlane{faulttype} <- isonyz(faulttype)

isonyz(::Type{NormalFault}) = true
isonyz(::Type{ThrustFault}) = true
isonyz(::Type{StrikeSlipFault}) = false

abstract type AbstractDifferenceGrids{dim, isuniform, T} end
abstract type BoundaryElementGrid{dim, isuniform, T} <: AbstractDifferenceGrids{dim, isuniform, T} end

struct BEMGrid_1D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{1, true, T}
    ξ::U1 # along-downdip
    Δξ::T
    nξ::I
    aξ::U2 # cell boundary along-downdip (width)
end

"Along-strike will be places on x-axis while along-downdip on yz-plane, w.r.t Okada's dc3d coordinates."
struct BEMGrid_2D{U1<:AbstractVector, U2<:AbstractMatrix, T<:Number, I<:Integer} <: BoundaryElementGrid{2, true, T}
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
    α::T = (λ + μ) / (λ + 2μ) # material constants which equals (λ + μ) / (λ + 2μ)
end

@with_kw struct PlaneMaterialProperties{dim, U<:AbstractVecOrMat, P<:AbstractArray, T<:Number} <: AbstractMaterialProperties{dim}
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    k::P # stiffness tensor
    σ::U # effective normal stress
    η::U # radiation damping
    vpl::T # plate rate, unlike pure Rate-State Friction simulation, here is restrained to be constant
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    function PlaneMaterialProperties(a::U, b::U, L::U, k::P, σ::U, η::U, vpl::T, f0::T, v0::T) where {U, P, T}
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

@inline @traitfn function shear_traction(::Type{FT}, u, λ, μ, dip) where {FT<:PlaneFault; IsOnYZPlane{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline @traitfn function shear_traction(::Type{FT}, u, λ, μ, dip) where {FT<:PlaneFault; !IsOnYZPlane{FT}}
    # Special case for `dip = 90.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

"""
Periodic boundary condition for 2D faults.

# Arguments
- same as `dc3d_okada`, see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
- `ax::T`: along-strike fault length
- `nrept::Integer`: number of repetition
- `buffer::T`: length of buffer size for introducing *zero-dislocation* area at along-strike edges of defined fault domain.
"""
function stiffness_periodic_boundary_condition(x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    ax::T, nrept::Integer, buffer::T=ax) where {T<:Number}
    u = zeros(T, 12)
    ax_extend = ax + buffer
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i * ax_extend, aw, disl)
    end
    return u
end

"""
    stiffness_tensor(fa::PlaneFaultDomain, gd::BEMGrid, mp::BEM, ::Vararg)

# Arguments
- ax_ratio: ratio of along-strike length against along downdip for mimicing infinitely extending 1D fault

# Note
- Faults are originated from surface and extends downwards, thus `dep = 0`

Calculate the reduced stiffness tensor based on the symmetric properties for 2D case.
Results for 1D fault will be 2D matrix and 3d matrix for 2D fault.
"""
function stiffness_tensor(fa::PlaneFaultDomain{ftype, 1, T}, gd::BoundaryElementGrid{1, true}, ep::HomogeneousElasticProperties,
    ax_ratio=50) where {ftype<:PlaneFault, T<:Number}

    if nprocs() == 1
        y, z = gd.ξ .* cospi(fa.dip / 180), gd.ξ .* sinpi(fa.dip / 180)
        disl = applied_unit_dislocation(ftype)
        ax = fa[:ξ] * ax_ratio .* [-one(T), one(T)]
        u12 = zeros(T, 12)
        ST = zeros(T, gd.nξ, gd.nξ)
        @threads for j = 1: gd.nξ
            @inbounds @fastmath @simd for i = 1: gd.nξ
                u12 .= dc3d_okada(0., y[i], z[i], ep.α, 0., fa.dip, ax, gd.aξ[j,:], disl)
                ST[i,j] = shear_traction(ftype, u12, ep.λ, ep.μ, fa.dip)
            end
        end
        return ST
    else
        ST = SharedArray{T}(gd.nξ, gd.nξ)
    end
end

"""
For noraml fault, it should of course be [0., -1., 0.]. However, in term of force balance, it is quivalent to thrust fault
    if dip angle are constrained within [0, π/2] in fact.

The unit of *unit dislocation* below is the same of `v * t` at set by user so to avoid normalization step.
"""
applied_unit_dislocation(::Type{NormalFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{ThrustFault}) = [0., 1., 0.]
applied_unit_dislocation(::Type{StrikeSlipFault}) = [1., 0., 0.]

"Temporal variable in solving ODEs aimed to avoid allocation overheads."
@with_kw struct TmpVariable{T<:AbstractVecOrMat}
    dμ_dt::T
    dμ_dθ::T
    dμ_dv::T
    ψ1::T
    ψ2::T
    relv::T
    relv_dft::Union{T, Nothing}
    dμ_dt_dft::Union{T, Nothing}
    plan_1d::Union{Plan, Nothing}
    plan_2d::Union{Plan, Nothing}
end

function create_tmp_var(gd::BoundaryElementGrid{dim, true, T}) where {dim, T}
    if dim == 1
        tvar = TmpVariable([Vector{T}(undef, gd.nξ) for _ in 1: 6]..., fill(nothing, 4)...)
    elseif dim == 2
        FFTW.set_num_threads(Sys.CPU_THREADS)
        x1 = Vector{T}(undef, gd.nx)
        x2 = Matrix{T}(undef, gd.nx, gd.nξ)
        p1 = plan_rfft(x1, flags=FFTW.MEASURE)
        p2 = plan_rfft(x2, 1, flags=FFTW.MEASURE)
        rfft_size = floor(Int, gd.nx/2) + 1
        tvar = TmpVariable(
            [Matrix{T}(undef, gd.nx, gd.nξ) for _ in 1: 6]...,
            [Matrix{T}(undef, rfft_size, gd.nξ) for _ in 1: 2]...,
            p1, p2)
    else
        @warn "Unexpected dim: $dim received."
        return nothing
    end
    return tvar
end

function derivations!(du, u, mp::PlaneMaterialProperties{dim}, tvar::TmpVariable, se, fform) where {dim}
    _ndim = dim + 1
    v = selectdim(u, _ndim, 1)
    θ = selectdim(u, _ndim, 2)
    dv = selectdim(du, _ndim, 1)
    dθ = selectdim(du, _ndim, 2)

    @. tvar.relv = mp.vpl - v
    dθ_dt!(dθ, se, v, θ, mp.L)
    dμ_dt!(mp, tvar)
    dv_dθ_dt!(fform, dv, dθ, v, θ, mp, tvar)
end

@inline function dμ_dt!(mp::PlaneMaterialProperties{1}, tvar::TmpVariable)
    mul!(tvar.dμ_dt, mp.k, tvar.relv)
end

@inline function dμ_dt!(mp::PlaneMaterialProperties{2}, tvar::TmpVariable)

end

@inline dθ_dt!(dθ, se, v, θ, L) =  @. dθ = dθ_dt(se, v, θ, L)

@inline function dv_dθ_dt!(::CForm, dv, dθ, v, θ, mp::PlaneMaterialProperties, tvar::TmpVariable)
    @. tvar.dμ_dθ = mp.σ * mp.b / θ
    @. tvar.dμ_dv = mp.σ * mp.a / v
    @. dv = dv_dt(tvar.dμ_dt, tvar.dμ_dv, tvar.dμ_dθ, dθ, mp.η)
end

@inline function dv_dθ_dt!(::RForm, dv, dθ, v, θ, mp::PlaneMaterialProperties, tvar::TmpVariable, se, fform)
    @. tvar.ψ1 = exp((mp.f0 + mp.b * log(mp.v0 * θ / mp.L)) / mp.a) / 2mp.v0
    @. tvar.ψ2 = mp.σ * tvar.ψ1 / sqrt(1 + (v * tvar.ψ1)^2)
    @. tvar.dμ_dv = mp.a * tvar.ψ2
    @. tvar.dμ_dθ = mp.b / θ * v * tvar.ψ2
    @. dv = dv_dt(tvar.dμ_dt, tvar.dμ_dv, tvar.dμ_dθ, dθ, mp.η)
end

function EarthquakeCycleProblem(gd::BoundaryElementGrid, p::MaterialProperties{dim}, u0, tspan; se=DieterichStateLaw(), fform=CForm()) where {dim}
    (fform == RForm() && p.η ≈ 0.0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    tvar = create_tmp_var(gd)
    f! = (du, u, p, t) -> derivations!(du, u, p, tvar, se, fform)
    ODEProblem(f!, u0, tspan, p)
end
